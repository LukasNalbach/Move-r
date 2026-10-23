/**
 * part of LukasNalbach/Move-r
 *
 * MIT License
 *
 * Copyright (c) Lukas Nalbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef MOVE_R_MOVE_RB_LOCATE_CONTEXT_TPP
#define MOVE_R_MOVE_RB_LOCATE_CONTEXT_TPP

#include <move_rb/locate_context.hpp>

template <move_r_support support>
template <locate_output_t output>
inline uint64_t move_rb<support>::locate_context_t<output>::sampled_occ() const
{
    const auto& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s_kept] = ctx->template const_vars<LEFT>();
    const sample_info_pack_t s = ctx->sample();

    if (s.d() == LEFT) {
        return (s.t() == BEG ? idx->idx_fwd.SA_run_start(s.x()) : idx->idx_fwd.SA_run_end(s.x())) - s.i();
    } else if (s.d() == NO_DIR) [[unlikely]] {
        return s.x(); // an explicit SA-sample stores its value in x
    } else /* if (s.d() == RIGHT) */ {
        uint64_t val = s.t() == BEG ? idx->idx_bwd.SA_run_start(s.x()) : idx->idx_bwd.SA_run_end(s.x());
        return idx->n - val - (ctx->m - s.i() + 1);
    }
}

template <move_r_support support>
template <locate_output_t output>
inline uint64_t move_rb<support>::locate_context_t<output>::compute_center()
    requires(supports_multiple_locate)
{
    if (centered) return c;
    centered = true;

    const auto& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s_kept] = ctx->template const_vars<LEFT>();
    const sample_info_pack_t s = ctx->sample();

    SA_c = sampled_occ();

    if (s.d() == LEFT || (stores_SA_pos && s.d() == NO_DIR)) {
        c = b + s.o();
    } else if constexpr (!stores_SA_pos) {
        /* the position of the sample is not stored; apply Phi^{-1} to SA_c until the LCP value is
           smaller than m, i.e., until e is reached; the decoded occurrences are stored in buf
           (if reports_sa_positions) */
        const auto& M_Phi_inv = idx->idx_fwd.M_Phi_inv();
        uint64_t sa = SA_c;
        uint64_t it = idx->phi_inv_interval(sa);

        if constexpr (reports_sa_positions) { buf.clear(); buf_pos = 0; }
        uint64_t num_steps = 0;

        while (num_steps + 1 < ctx->num_occ() &&
            idx->plcp_at_phi_inv(sa, it) >= ctx->m)
        {
            M_Phi_inv.move(sa, it);
            num_steps++;
            if constexpr (reports_sa_positions) buf.emplace_back(sa);
        }

        c = e - num_steps;
    } else /* if (s.d() == RIGHT) */ {
        /* SA[b] is decoded from the SA-sample at the start of the BWT run containing b, if that reads fewer
           values than the m - i LF-steps from the position of the backward SA-sample (see SA_pos) */
        const auto& fwd = idx->idx_fwd;
        const uint64_t steps = ctx->m - s.i();
        const uint64_t x_b = idx->template mlf_interval<LEFT>(b);
        const uint64_t x_run = fwd.first_interval_of_run(x_b);
        const uint64_t p_run = fwd.M_LF().p(x_run);

        if (b - p_run < steps) {
            SA_c = fwd.SA_run_start(x_run);

            if (p_run < b) {
                auto dec = fwd.SA_rlz().decode_context();
                dec.init_right(p_run + 1);
                dec.set_value(SA_c);
                dec.skip_right(b + 1);
                SA_c = dec.value();
            }

            c = b;
        } else {
            c = idx->SA_pos(s.x(), s.t() == BEG);
            uint64_t c_ = idx->template mlf_interval<LEFT>(c);
            assert(fwd.M_LF().p(c_) <= c && c < fwd.M_LF().p(c_ + 1));

            for (uint64_t j = 0; j < steps; j++) {
                fwd.M_LF().move(c, c_);
            }
        }
    }

    assert(b <= c && c <= e);
    return c;
}

template <move_r_support support>
template <locate_output_t output>
inline uint64_t move_rb<support>::locate_context_t<output>::first_occ()
    requires(supports_multiple_locate)
{
    const auto& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = ctx->template const_vars<LEFT>();

    if constexpr (!stores_SA_pos) {
        // start at the SA-sample; the occurrences after it are decoded with Phi^{-1} (see next_occ)
        SA_c = sampled_occ();
        SA_i = SA_c;
        scan_dir = RIGHT;
        x_phi = idx->phi_inv_interval(SA_i);
    } else {
        // start at c (see compute_center)
        compute_center();
        SA_i = SA_c;
        scan_dir = RIGHT;
        i = c;
        SA_rlz_r = idx->idx_fwd.SA_rlz().decode_context();
        SA_rlz_r.set_value(SA_i);
        if (c < e) SA_rlz_r.init_right(c + 1);
    }

    occ_rem--;
    return SA_i;
}

template <move_r_support support>
template <locate_output_t output>
inline bool move_rb<support>::locate_context_t<output>::is_right_limit() const
    requires(supports_multiple_locate && !stores_SA_pos)
{
    // PLCP[Phi^{-1}(SA_i)] = LCP[i + 1] < m <=> i = e
    return idx->plcp_at_phi_inv(SA_i, x_phi) < ctx->m;
}

template <move_r_support support>
template <locate_output_t output>
inline uint64_t move_rb<support>::locate_context_t<output>::next_occ()
    requires(supports_multiple_locate)
{
    if (scan_dir == NO_DIR) [[unlikely]] {
        return first_occ();
    }

    const auto& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = ctx->template const_vars<LEFT>();

    if constexpr (!stores_SA_pos) {
        // after the last occurrence to the right of the SA-sample, continue to the left of it with Phi
        if (scan_dir == RIGHT && is_right_limit()) [[unlikely]] {
            scan_dir = LEFT;
            SA_i = SA_c;
            x_phi = idx->phi_interval(SA_i);
        }

        if (scan_dir == RIGHT) {
            idx->idx_fwd.M_Phi_inv().move(SA_i, x_phi);
        } else {
            idx->idx_fwd.M_Phi().move(SA_i, x_phi);
        }
    } else {
        if (scan_dir == RIGHT && i == e) [[unlikely]] {
            scan_dir = LEFT;
            i = c;
            SA_rlz_l = idx->idx_fwd.SA_rlz().decode_context();
            SA_rlz_l.set_value(SA_c);
            SA_rlz_l.init_left(c);
        }

        if (scan_dir == RIGHT) {
            SA_rlz_r.next();
            SA_i = SA_rlz_r.value();
            i++;
        } else {
            SA_rlz_l.prev();
            SA_i = SA_rlz_l.value();
            i--;
        }
    }

    occ_rem--;
    return SA_i;
}

template <move_r_support support>
template <locate_output_t output>
template <typename report_fnc_t>
inline void move_rb<support>::locate_context_t<output>::locate(report_fnc_t report)
    requires(supports_multiple_locate)
{
    locate_from_sample(report);
}

template <move_r_support support>
template <locate_output_t output>
template <typename report_fnc_t>
inline void move_rb<support>::locate_context_t<output>::locate_from_sample(
    report_fnc_t report) requires(supports_multiple_locate)
{
    const auto& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = ctx->template const_vars<LEFT>();

    // reports an occurrence (with its position in the suffix array, if reports_sa_positions)
    auto emit = [&](uint64_t position, uint64_t occurrence) {
        if constexpr (reports_sa_positions) report(position, occurrence);
        else report(occurrence);
    };

    /* without positions: report the SA-sample, then the occurrences to the right of it with Phi^{-1} (until
       the LCP value is smaller than m), then the remaining occ_rem occurrences to the left of it with Phi */
    if constexpr (!reports_sa_positions && !stores_SA_pos) {
        if (occ_rem == ctx->num_occ()) {
            SA_c = sampled_occ();
            report(SA_c);
            occ_rem--;
        }

        SA_i = SA_c;
        uint64_t it = idx->phi_inv_interval(SA_i);

        while (occ_rem > 0 && idx->plcp_at_phi_inv(SA_i, it) >= ctx->m) {
            idx->idx_fwd.M_Phi_inv().move(SA_i, it);
            report(SA_i);
            occ_rem--;
        }

        if (occ_rem > 0) {
            SA_i = SA_c;
            it = idx->phi_interval(SA_i);

            while (occ_rem > 0) {
                idx->idx_fwd.M_Phi().move(SA_i, it);
                report(SA_i);
                occ_rem--;
            }
        }

        return;
    }

    // compute c and report SA[c] (see compute_center)
    if (occ_rem == ctx->num_occ()) {
        compute_center();
        emit(c, SA_c);
        occ_rem--;
    }

    if constexpr (!stores_SA_pos) {
        // report SA[b..c-1] from right to left
        if (c > b) {
            SA_i = SA_c;
            i = c;
            x_phi = idx->phi_interval(SA_i);

            while (i > b) {
                idx->idx_fwd.M_Phi().move(SA_i, x_phi);
                i--;
                emit(i, SA_i);
            }
        }

        // report SA[c+1..e] from buf (if compute_center stored it there), else with Phi^{-1}
        bool from_buffer = false;

        if constexpr (reports_sa_positions) {
            while (buf_pos < buf.size()) {
                emit(c + 1 + buf_pos, buf[buf_pos]);
                buf_pos++;
            }

            from_buffer = !buf.empty();
        }

        if (!from_buffer && c < e) {
            SA_i = SA_c;
            i = c;
            x_phi = idx->phi_inv_interval(SA_i);

            while (i < e) {
                idx->idx_fwd.M_Phi_inv().move(SA_i, x_phi);
                i++;
                emit(i, SA_i);
            }
        }
    } else {
        if (c > b) {
            SA_rlz_l = idx->idx_fwd.SA_rlz().decode_context();
            SA_rlz_l.set_value(SA_c);
            SA_rlz_l.init_left(c);
            SA_rlz_l.report_left(b, emit);
        }

        if (c < e) {
            SA_rlz_r = idx->idx_fwd.SA_rlz().decode_context();
            SA_rlz_r.set_value(SA_c);
            SA_rlz_r.init_right(c + 1);
            SA_rlz_r.report_right(e, emit);
        }
    }

    occ_rem = 0;
}

template <move_r_support support>
template <locate_output_t output>
std::vector<uint48_t> move_rb<support>::locate_context_t<output>::locate()
    requires(supports_multiple_locate)
{
    std::vector<uint48_t> Occ;
    Occ.reserve(occ_rem);

    if constexpr (reports_sa_positions) {
        locate([&](uint64_t, uint64_t occ){Occ.emplace_back(occ);});
    } else {
        locate([&](uint64_t occ){Occ.emplace_back(occ);});
    }

    return Occ;
}

#endif // MOVE_R_MOVE_RB_LOCATE_CONTEXT_TPP
