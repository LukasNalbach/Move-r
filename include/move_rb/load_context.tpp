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

#ifndef MOVE_R_MOVE_RB_LOAD_CONTEXT_TPP
#define MOVE_R_MOVE_RB_LOAD_CONTEXT_TPP

#include <data_structures/elias_fano_bit_vector.hpp>

#include <move_rb/move_rb.hpp>

template <move_r_support support>
template <query_support_t query_support>
void move_rb<support>::load_context(
    search_context_t<query_support>& ctx, uint64_t len, uint32_t idx, uint32_t idx_R) const
{
    ctx.idx = this;
    std::tie(ctx.b, ctx.e) = _context_tables.template interval<_forward_order>(len, idx);

    // both intervals have the same width
    ctx.b_R = _context_tables.template interval_beg<_reverse_order>(len, idx_R);
    ctx.e_R = ctx.b_R + (ctx.e - ctx.b);

    /* the tables store an SA-sample only for an entry whose interval contains no SA-sample of the index;
       else the sample is marked with d = BOTH_DIRS and computed on demand (see search_context_t::sample) */
    if constexpr (supports_locate && query_support == LOCATE) {
        uint32_t smpl_idx = _context_tables.sample_rank(len, idx);

        if (smpl_idx != no_ctx_entry) [[unlikely]] {
            const auto& smpl = _context_tables.samples(len);
            ctx.s = {
                direction_t(smpl.template get<_smpl_d>(smpl_idx)),
                sample_t(smpl.template get<_smpl_t>(smpl_idx)),
                smpl.template get<_smpl_o>(smpl_idx),
                smpl.template get<_smpl_i>(smpl_idx),
                smpl.template get<_smpl_x>(smpl_idx) };
        } else {
            ctx.s = { BOTH_DIRS, BEG, 0, 0, 0 };
        }
    }

    // the input intervals containing the limits are unknown (see valid_interval_dirs)
    ctx.b_ = 0;
    ctx.e_ = 0;
    ctx.b_R_ = 0;
    ctx.e_R_ = 0;

    ctx.m = len;
    ctx.err = 0;
    ctx.ctx_entry = idx;
    ctx.ctx_entry_R = idx_R;
    ctx.valid_interval_dirs = NO_DIR;
}

#endif // MOVE_R_MOVE_RB_LOAD_CONTEXT_TPP
