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

#pragma once

#include <rlzsa_opt/rlzsa_opt.hpp>

// ############################# DECODING #############################

template <typename pos_t>
inline void rlzsa_opt<pos_t>::decode_context_t::init_right(pos_t i)
{
    this->i = i;

    // index in SCP_S of the last sampled copy phrase starting before or at i
    pos_t x_scps = rlz->SCP_S().rank_1(i + 1);

    if (x_scps == 0) [[unlikely]] {
        // i lies before the first copy phrase
        s_np = i + 1;
        x_p = i;
        x_lp = i;
        x_cp = 0;
        x_r = rlz->SR(0);
    } else {
        // copy-phrase index of the last copy-phrase starting before or at i
        pos_t x_cp_lcp = (x_scps - 1) * sample_rate_scp;
        // starting position of the last copy-phrase starting before or at i
        pos_t s_lcp = rlz->SCP_S().select_1(x_scps);
        // phrase index of the last copy-phrase starting before or at i
        pos_t x_p_lcp = rlz->PT().select_0(x_cp_lcp + 1);
        // phrase index of the x_cp_lcp+1-th copy phrase
        pos_t x_p_ncp = rlz->PT().select_0(x_cp_lcp + 2);
        // number of literal phrases between the current and the next copy-phrase
        pos_t n_lp = x_p_ncp - x_p_lcp - 1;
        // starting position of the next copy-phrase
        pos_t s_ncp = s_lcp + rlz->CPL(x_cp_lcp) + n_lp;

        // find the last copy-phrase starting before or at i
        while (s_ncp <= i) {
            x_cp_lcp++;
            s_lcp = s_ncp;
            x_p_lcp = x_p_ncp;
            x_p_ncp = rlz->PT().select_0(x_cp_lcp + 2);
            n_lp = x_p_ncp - x_p_lcp - 1;
            s_ncp += rlz->CPL(x_cp_lcp) + n_lp;
        }

        if (i >= s_ncp - n_lp) {
            // there is a literal phrase at position i
            s_np = i + 1;
            x_cp = x_cp_lcp + 1;
            x_r = rlz->SR(x_cp);
            x_p = x_p_lcp + 1 + (i - (s_ncp - n_lp));
        } else {
            // i lies within a copy-phrase
            x_cp = x_cp_lcp;
            x_p = x_p_lcp;
            x_r = rlz->SR(x_cp) + (i - s_lcp);
            s_np = s_ncp - n_lp;
        }

        x_lp = x_p - x_cp;
    }
}

template <typename pos_t>
inline void rlzsa_opt<pos_t>::decode_context_t::init_left(pos_t i)
{
    init_right(i);
    turn_left();
}

template <typename pos_t>
inline void rlzsa_opt<pos_t>::decode_context_t::turn_left()
{
    // adjust the navigation state to advance to the left instead of to the right
    if (rlz->PT(x_p)) {
        s_np--;

        if (x_cp > 0) [[likely]] {
            x_cp--;
            x_r = rlz->SR(x_cp) + rlz->CPL(x_cp) - 1;
        }
    } else {
        s_np -= rlz->CPL(x_cp);
        x_lp--;
    }
}

template <typename pos_t>
inline void rlzsa_opt<pos_t>::decode_context_t::turn_right()
{
    // adjust the navigation state to advance to the right instead of to the left (the exact inverse of
    // turn_left())
    if (rlz->PT(x_p)) {
        s_np++;

        // x_cp is the previous copy-phrase (left-ready); recompute the next copy-phrase (right-ready)
        // directly from x_p, because a plain x_cp++ is ambiguous in the leading literal run (no previous
        // copy-phrase)
        x_cp = rlz->PT().rank_0(x_p);
        x_r = rlz->SR(x_cp);
    } else {
        s_np += rlz->CPL(x_cp);
        x_lp++;
    }
}

template <typename pos_t>
inline void rlzsa_opt<pos_t>::decode_context_t::skip_right(pos_t e)
{
    while (i < e) {
        // decode all copy-phrases before the next literal phrase
        while (i < e && !rlz->PT(x_p)) {
            // decode the x_cp-th copy-phrase
            while (i < s_np && i < e) {
                s += rlz->R(x_r);
                s -= rlz->n;
                i++;
                x_r++;
            }

            if (i < e || i == s_np) [[likely]] {
                x_p++;
                x_cp++;
                x_r = rlz->SR(x_cp);
                s_np += rlz->PT(x_p) ? 1 : rlz->CPL(x_cp);
            }
        }

        // decode all literal phrases before the next copy-phrase
        while (i < e && rlz->PT(x_p)) {
            // decode the x_lp-th literal phrase
            s += rlz->LP(x_lp);
            s -= rlz->n;
            i++;
            x_p++;
            x_lp++;

            if (rlz->PT(x_p)) [[likely]] {
                s_np++;
            } else {
                s_np += rlz->PT(x_p) ? 1 : rlz->CPL(x_cp);
            }
        }
    }
}

template <typename pos_t>
inline void rlzsa_opt<pos_t>::decode_context_t::skip_left(pos_t b)
{
    while (i > b) {
        // skip all copy-phrases after the next literal phrase
        while (i > b && !rlz->PT(x_p)) {
            // skip the x_cp-th copy-phrase
            while (i >= s_np && i > b) {
                s += rlz->n;
                s -= rlz->R(x_r);
                i--;
                x_r--;
            }

            if (i > b || i + 1 == s_np) [[likely]] {
                if (x_cp > 0) [[likely]] {
                    x_cp--;
                    x_r = rlz->SR(x_cp) + rlz->CPL(x_cp) - 1;
                }

                x_p--;
                s_np -= rlz->PT(x_p) ? 1 : rlz->CPL(x_cp);
            }
        }

        // skip all literal phrases after the next copy-phrase
        while (i > b && rlz->PT(x_p)) {
            // skip the x_lp-th literal phrase
            s += rlz->n;
            s -= rlz->LP(x_lp);
            i--;
            x_p--;
            x_lp--;

            if (rlz->PT(x_p)) [[likely]] {
                s_np--;
            } else {
                s_np -= rlz->CPL(x_cp);
            }
        }
    }
}

template <typename pos_t>
inline pos_t rlzsa_opt<pos_t>::decode_context_t::prev()
{
    if (rlz->PT(x_p)) {
        // literal phrase
        s += rlz->n;
        s -= rlz->LP(x_lp);
        i--;
        x_p--;
        x_lp--;

        if (i > 0) [[likely]] {
            if (rlz->PT(x_p)) {
                // the previous phrase is a literal phrase
                s_np--;
            } else {
                // the previous phrase is a copy-phrase
                s_np -= rlz->CPL(x_cp);
            }
        }
    } else {
        // copy-prhase
        s += rlz->n;
        s -= rlz->R(x_r);
        i--;
        x_r--;

        // i lies within the previous phrase
        if (s_np > 0 && i < s_np) [[unlikely]] {
            x_p--;

            if (x_cp > 0) [[likely]] {
                x_cp--;
                x_r = rlz->SR(x_cp) + rlz->CPL(x_cp) - 1;
            }

            if (rlz->PT(x_p)) {
                // the next phrase is a literal phrase
                s_np--;
            } else {
                // the next phrase is a copy-phrase
                s_np -= rlz->CPL(x_cp);
            }
        }
    }

    return s;
}

template <typename pos_t>
inline pos_t rlzsa_opt<pos_t>::decode_context_t::next()
{
    if (rlz->PT(x_p)) {
        // literal phrase
        s += rlz->LP(x_lp);
        s -= rlz->n;
        i++;
        x_p++;
        x_lp++;

        if (i < rlz->n) [[likely]] {
            if (rlz->PT(x_p)) {
                // the next phrase is a literal phrase
                s_np++;
            } else {
                // the next phrase is a copy-phrase
                s_np += rlz->CPL(x_cp);
            }
        }
    } else {
        // copy-prhase
        s += rlz->R(x_r);
        s -= rlz->n;
        i++;
        x_r++;

        // there is a new phrase starting at i
        if (s_np <= i && i < rlz->n) [[unlikely]] {
            x_p++;
            x_cp++;
            x_r = rlz->SR(x_cp);

            if (rlz->PT(x_p)) {
                // the next phrase is a literal phrase
                s_np++;
            } else {
                // the next phrase is a copy-phrase
                s_np += rlz->CPL(x_cp);
            }
        }
    }

    return s;
}

template <typename pos_t>
template <typename report_fnc_t>
inline void rlzsa_opt<pos_t>::decode_context_t::report_left(pos_t b, report_fnc_t report)
{
    static constexpr bool report_pos = function_traits<report_fnc_t>::arity > 1;

    while (true) {
        // decode all copy-phrases after the previous literal phrase
        while (!rlz->PT(x_p)) {
            // decode the x_cp-th copy-phrase
            while (i >= s_np) {
                s += rlz->n;
                s -= rlz->R(x_r);
                i--;
                if constexpr (report_pos) report(i, s); else report(s);
                if (i == b) [[unlikely]] return;
                x_r--;
            }

            if (x_cp > 0) [[likely]] {
                x_cp--;
                x_r = rlz->SR(x_cp) + rlz->CPL(x_cp) - 1;
            }

            x_p--;
            s_np -= rlz->PT(x_p) ? 1 : rlz->CPL(x_cp);
        }

        // decode all literal phrases after the previous copy-phrase
        while (rlz->PT(x_p)) {
            // decode the x_lp-th literal phrase
            s += rlz->n;
            s -= rlz->LP(x_lp);
            i--;
            if constexpr (report_pos) report(i, s); else report(s);
            if (i == b) [[unlikely]] return;
            x_p--;
            x_lp--;
            s_np--;
        }

        // set s_np to the starting position of the last (the x_lp-th)
        // literal phrase before the current (the x_cp-th) copy-phrase
        s_np -= rlz->CPL(x_cp) - 1;
    }
}

template <typename pos_t>
template <typename report_fnc_t>
inline void rlzsa_opt<pos_t>::decode_context_t::report_right(pos_t e, report_fnc_t report)
{
    static constexpr bool report_pos = function_traits<report_fnc_t>::arity > 1;

    while (true) {
        // decode all copy-phrases before the next literal phrase
        while (!rlz->PT(x_p)) {
            // decode the x_cp-th copy-phrase
            while (i < s_np) {
                s += rlz->R(x_r);
                s -= rlz->n;
                if constexpr (report_pos) report(i, s); else report(s);
                if (i == e) [[unlikely]] return;
                i++;
                x_r++;
            }

            x_p++;
            x_cp++;
            x_r = rlz->SR(x_cp);
            s_np += rlz->PT(x_p) ? 1 : rlz->CPL(x_cp);
        }

        // decode all literal phrases before the next copy-phrase
        while (rlz->PT(x_p)) {
            // decode the x_lp-th literal phrase
            s += rlz->LP(x_lp);
            s -= rlz->n;
            if constexpr (report_pos) report(i, s); else report(s);
            if (i == e) [[unlikely]] return;
            i++;
            x_p++;
            x_lp++;
            s_np++;
        }

        // set s_np to the starting position of the next (the x_lp-th)
        // literal phrase after the current (the x_cp-th) copy-phrase
        s_np += rlz->CPL(x_cp) - 1;
    }
}
