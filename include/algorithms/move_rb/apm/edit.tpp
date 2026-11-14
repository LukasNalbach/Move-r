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

/**
 * @file edit.tpp
 * @brief the methods of apm_edit: the search, the errors of the found contexts and the occurrences at the text ends
 */

#ifndef MOVE_R_ALGORITHMS_MOVE_RB_APM_EDIT_TPP
#define MOVE_R_ALGORITHMS_MOVE_RB_APM_EDIT_TPP

#include <algorithms/move_rb/apm/edit.hpp>

// ############################# CLUSTERS #############################

template <typename idx_t>
template <typename sep_fnc_t>
auto apm_edit<idx_t>::cluster_t::local_minima(sep_fnc_t& at_separator) -> std::vector<node_t>
{
    std::vector<node_t> minima;
    minima.reserve(_last_cell + 1);

    for (uint64_t i = 0; i <= _last_cell; i++) {
        if (!_nodes[i].taken() && _dists[i] <= _k_max && is_local_min(i, at_separator)) {
            minima.emplace_back(take_cell(i));
        }
    }

    return minima;
}

template <typename idx_t>
template <typename sep_fnc_t>
auto apm_edit<idx_t>::cluster_t::choose_center(uint16_t k_min, std::vector<node_t>& descendants,
    std::vector<uint16_t>& init_dists, sep_fnc_t& at_separator) -> std::tuple<node_t, bool>
{
    descendants.reserve(_dists.size());
    init_dists.reserve(_dists.size());
    node_t node;
    uint64_t i = _last_cell + 1;

    for (uint64_t j = 0; j <= _last_cell; j++) {
        if (_dists[j] < k_min || _dists[j] > _k_max) continue;
        if (is_local_min(j, at_separator)) { i = j; break; }
    }

    if (i <= _last_cell) {
        // the cell is taken even if it was taken before (take_cell returns a new copy)
        node = take_cell(i);
        init_dists.emplace_back(_dists[i]);

        for (uint64_t j = i + 1; j <= _last_cell; j++) {
            descendants.emplace_back(_nodes[j]);
            init_dists.emplace_back(_dists[j]);
        }

        /* raise every local minimum x below k_min: [lo, hi] is the range of cells around x in which the distances
           decrease by 1 per cell towards x. It starts after cell 0 (the center, whose distance is the error count of
           node) and after the last cell before x with at least k_min errors whose string has an occurrence followed by
           the separator (see ends_at_separator; this occurrence is not an occurrence of the string of x, so its
           alignment cannot be continued through x). Its distances are replaced by distances that grow by 1 per cell
           from the cells next to the range (at most _k_max + 1); if the range reaches the last cell, they grow from
           the cell left of the range */
        const uint64_t last = init_dists.size() - 1;

        for (uint64_t x = 1; x <= last; x++) {
            if (init_dists[x] < k_min && init_dists[x] <= init_dists[x - 1] && (x == last || init_dists[x] <= init_dists[x + 1])) {
                uint64_t lo = 1;
                uint64_t hi = last;

                for (uint64_t l = x; l-- > 1;) {
                    if (init_dists[l] != init_dists[l + 1] + 1 || (init_dists[l] >= k_min && ends_at_separator(i + l, at_separator))) {
                        lo = l + 1;
                        break;
                    }
                }

                for (uint64_t l = x + 1; l <= last; l++) {
                    if (init_dists[l] != init_dists[l - 1] + 1) {
                        hi = l - 1;
                        break;
                    }
                }

                if (hi != last) {
                    const uint64_t left_dist = init_dists[lo - 1];
                    const uint64_t right_dist = init_dists[hi + 1];

                    for (uint64_t l = lo; l <= hi; l++) {
                        const uint64_t from_left = left_dist + (l - (lo - 1));
                        const uint64_t from_right = right_dist + ((hi + 1) - l);
                        init_dists[l] = std::min<uint64_t>(_k_max + 1, std::min(from_left, from_right));
                    }
                } else {
                    for (uint64_t l = lo; l <= last; l++) init_dists[l] = init_dists[l - 1] + 1;
                }
            }
        }

        return {node, true};
    }

    return {node, false};
}

// ############################# TEXT ENDS #############################

template <typename idx_t>
template <typename sub_t>
auto apm_edit<idx_t>::text_end(const sub_t& sub, uint64_t length) -> inp_t
{
    inp_t end(length);

    // the BWT position pos in the input interval x of M_LF; the BWT symbol at position 0 is the last symbol of the text
    uint64_t x = 0, pos = 0;

    for (uint64_t i = length; i-- > 0;) {
        end[i] = sub.unmap_symbol(sub.L_(x));
        if (i != 0) sub.M_LF().move(pos, x);
    }

    return end;
}

template <typename idx_t>
auto apm_edit<idx_t>::text_end_windows(uint64_t m, uint64_t k) const -> std::tuple<uint64_t, inp_t, inp_t>
{
    const uint64_t n = _index.n - 1;
    const uint64_t w = std::min<uint64_t>(n, m + k); // the length of a window
    const auto& seqs = _index.seq_data();
    const uint64_t left_beg = seqs.has_sequences() ? seqs.sequence_start(0) : 0;
    const uint64_t left_end = left_beg + (seqs.has_sequences() ? std::min<uint64_t>(seqs.sequence_length(0), w) : w);
    const uint64_t right_beg = seqs.has_sequences()
        ? std::max<uint64_t>(seqs.sequence_start(seqs.num_sequences() - 1), n - w) : n - w;

    // the first window is decoded from the end of T^R, without the symbols before the first sequence
    inp_t left = text_end(_index.backward_index(), left_end);
    std::reverse(left.begin(), left.end());
    left.erase(left.begin(), left.begin() + left_beg);
    inp_t right = text_end(_index.forward_index(), n - right_beg);

    // an occurrence contains no separator (of a masked run): the first window ends before its first separator, the
    // last window starts after its last one
    auto is_sep = [&](uint8_t c) { return this->is_separator(c); };
    left.erase(std::find_if(left.begin(), left.end(), is_sep), left.end());
    right.erase(right.begin(), std::find_if(right.rbegin(), right.rend(), is_sep).base());
    return {left_beg, std::move(left), std::move(right)};
}

template <typename idx_t>
template <cigar_mode_t mode, typename report_fnc_t>
void apm_edit<idx_t>::locate_at_text_ends(const inp_t& P, uint64_t m, uint64_t k, report_fnc_t& report,
    std::vector<cigar_t<uint8_t>>& cigars) const
{
    with_edit_distance_word(k, [&]<typename word_t>() { locate_at_text_ends_impl<mode, word_t>(P, m, k, report, cigars); });
}

template <typename idx_t>
template <cigar_mode_t mode, typename word_t, typename report_fnc_t>
void apm_edit<idx_t>::locate_at_text_ends_impl(const inp_t& P, uint64_t m, uint64_t k, report_fnc_t& report,
    std::vector<cigar_t<uint8_t>>& cigars) const
{
    const uint64_t n = _index.n - 1;
    if (m == 0 || n == 0) return;
    auto [left_beg, left, right] = text_end_windows(m, k);
    static thread_local edit_distance_matrix<word_t> mat, mat_any;
    const auto& ranks = _index.forward_index().map_int();

    // whether an alignment of Q against a prefix of W that starts with a deletion has fewer than err errors (W without
    // its first symbol then has fewer errors)
    auto dominated = [&](const inp_t& Q, const inp_t& W, bool first_matches, uint64_t err) {
        if (first_matches) return false;
        mat_any.set_input(directional_substring<inp_t>(Q, 0, m - 1, RIGHT), _index.sigma, ranks);
        return edit_dist_prefix(mat_any, W, m, k).first < err;
    };

    /* reports the alignment of Q against a prefix of W with the fewest errors that does not start with a deletion (see
       recompute_errors); pos_of maps its length to its position, and reverse_cigar is true for the reversed last window */
    auto report_window = [&](const inp_t& Q, const inp_t& W, auto pos_of, bool reverse_cigar) {
        const bool first_matches = !W.empty() && W[0] == Q[0];

        if constexpr (mode == CIGAR) {
            auto [err, len, cig] = edit_cigar_prefix(Q, m, W, W.size(), k);
            if (err > k || dominated(Q, W, first_matches, err)) return;
            if (reverse_cigar) std::reverse(cig.begin(), cig.end());
            cigars.push_back(std::move(cig));
            report({pos_of(len), len, err, uint32_t(cigars.size() - 1)});
        } else {
            mat.set_input(directional_substring<inp_t>(Q, 1, m - 1, RIGHT), _index.sigma, ranks);
            auto [err, len] = edit_dist_anchored_prefix(mat, W, m, k, first_matches);
            if (err > k || dominated(Q, W, first_matches, err)) return;
            report({pos_of(len), len, err});
        }
    };

    report_window(P, left, [left_beg](uint64_t) { return left_beg; }, false);
    inp_t P_rev(P.rbegin(), P.rend()), right_rev(right.rbegin(), right.rend());
    report_window(P_rev, right_rev, [n](uint64_t len) { return n - len; }, true);
}

// ############################# MATCHED STRINGS #############################

template <typename idx_t>
void apm_edit<idx_t>::dec_ref(uint64_t idx) const
{
    while (idx != no_piece) {
        uint16_t r = _matched_pieces[idx].refs;
        if (r == UINT16_MAX) break;
        if (r > 1) { _matched_pieces[idx].refs = r - 1; break; }
        _matched_pieces[idx].refs = 0;
        uint64_t prev = _matched_pieces[idx].prev();
        _free_pieces.emplace_back(idx);
        idx = prev;
    }
}

template <typename idx_t>
uint64_t apm_edit<idx_t>::append_literal(uint64_t prev, uint8_t sym, direction_t dir) const
{
    uint64_t i = allocate_piece();
    _matched_pieces[i] = matched_piece_t{prev, 0, 0, 0, sym, dir};
    if (prev != no_piece) inc_ref(prev);
    return i;
}

template <typename idx_t>
uint64_t apm_edit<idx_t>::append_copy(uint64_t prev, uint64_t pos, uint64_t len, direction_t dir) const
{
    uint64_t i = allocate_piece();
    _matched_pieces[i] = matched_piece_t{prev, pos, uint32_t(len), 0, uint8_t{}, dir};
    if (prev != no_piece) inc_ref(prev);
    return i;
}

template <typename idx_t>
void apm_edit<idx_t>::append_pattern_symbol(node_t& node, uint64_t pos, direction_t dir) const
{
    bool continues_copy = node.copy_len > 0 && node.copy_dir == dir &&
        (dir == RIGHT ? pos == node.copy_pos + node.copy_len : pos + 1 == node.copy_pos);
    if (continues_copy) { if (dir == LEFT) node.copy_pos = pos; node.copy_len++; }
    else { commit_copy(node); node.copy_pos = pos; node.copy_len = 1; node.copy_dir = dir; }
}

template <typename idx_t>
void apm_edit<idx_t>::append_last_symbol(node_t& node, const directional_substring<inp_t>& part, direction_t dir) const
{
    uint8_t sym = node.ctx.last_added_symbol();
    const inp_t& P = part.underlying();
    const uint64_t m = P.size();

    // append sym to the pending copy if sym equals the symbol of P that continues the copy
    if (node.copy_len > 0 && node.copy_dir == dir) {
        uint64_t next_pos = dir == RIGHT ? node.copy_pos + node.copy_len : node.copy_pos - 1;
        if (next_pos < m && P[next_pos] == sym) { if (dir == LEFT) node.copy_pos = next_pos; node.copy_len++; return; }
    }

    /* diag_pos: the position of P on the diagonal of the row of the node; j: the position closest to diag_pos with
       P[j] = sym and |j - diag_pos| <= _max_diagonal_shift (diag_pos + d before diag_pos - d), or no_piece if there is
       none */
    const int64_t diag_pos = int64_t(part.pos(node.match_len() - 1));
    auto match_at = [&](int64_t x) { return x >= 0 && x < int64_t(m) && P[x] == sym; };
    uint64_t j = no_piece;

    if (match_at(diag_pos)) {
        j = uint64_t(diag_pos);
    } else {
        for (int64_t d = 1; d <= int64_t(_max_diagonal_shift) && j == no_piece; d++) {
            if (match_at(diag_pos + d)) j = uint64_t(diag_pos + d);
            else if (match_at(diag_pos - d)) j = uint64_t(diag_pos - d);
        }
    }

    commit_copy(node);
    if (j != no_piece) { node.copy_pos = j; node.copy_len = 1; node.copy_dir = dir; }
    else node.set_last_piece(append_literal(node.last_piece_idx(), sym, dir));
}

template <typename idx_t>
void apm_edit<idx_t>::decode_matched_string(const inp_t& P, uint64_t last, inp_t& out) const
{
    // the indices of the pieces of the chain, from the last one to the first one
    static thread_local std::vector<uint64_t> pieces;
    pieces.clear();

    // total: the length of the string; n_left: the number of its symbols in prepended pieces
    uint64_t total = 0, n_left = 0;

    for (uint64_t i = last; i != no_piece; i = _matched_pieces[i].prev()) {
        const matched_piece_t& piece = _matched_pieces[i];
        pieces.emplace_back(i);
        uint64_t piece_len = piece.len > 0 ? piece.len : 1;
        total += piece_len;
        if (piece.dir == LEFT) n_left += piece_len;
    }

    // the prepended pieces fill out[0, n_left) from right to left, the appended ones out[n_left, total)
    out.resize(total);
    uint64_t lo = n_left, hi = n_left;

    for (auto it = pieces.rbegin(); it != pieces.rend(); ++it) {
        const matched_piece_t& piece = _matched_pieces[*it];
        if (piece.dir == LEFT) {
            if (piece.len > 0) { lo -= piece.len; for (uint64_t j = 0; j < piece.len; j++) out[lo + j] = P[piece.pos() + j]; }
            else out[--lo] = piece.sym;
        } else {
            if (piece.len > 0) { for (uint64_t j = 0; j < piece.len; j++) out[hi + j] = P[piece.pos() + j]; hi += piece.len; }
            else out[hi++] = piece.sym;
        }
    }
}

template <typename idx_t>
template <cigar_mode_t mode, typename word_t>
void apm_edit<idx_t>::recompute_errors(node_set_t& ctxts, const inp_t& P, uint64_t m, uint64_t k,
    std::vector<context_ref_t<mode>>& ctxts_sorted, std::vector<cigar_t<uint8_t>>& cigars) const
{
    static thread_local edit_distance_matrix<word_t> mat, mat_any;

    // mat: the alignments that do not start with a deletion (see edit_dist_anchored_prefix); mat_any: all alignments
    const auto& ranks = _index.forward_index().map_int();
    if constexpr (mode == NO_CIGAR) mat.set_input(directional_substring<inp_t>(P, 1, m - 1, RIGHT), _index.sigma, ranks);
    mat_any.set_input(directional_substring<inp_t>(P, 0, m - 1, RIGHT), _index.sigma, ranks);
    inp_t match;

    for (node_t& node : ctxts) {
        decode_matched_string(P, node.last_piece_idx(), match);
        const bool first_matches = !match.empty() && match[0] == P[0];

        /* a context is not reported if an alignment against a prefix of its string that starts with a deletion has
           fewer errors (its string without the first symbol has fewer errors); if P[0] matches the first symbol, an
           alignment that starts with a deletion has at least as many errors as one that does not */
        auto dominated = [&](uint64_t e) { return !first_matches && edit_dist_prefix(mat_any, match, m, k).first < e; };

        uint64_t err, len;
        if constexpr (mode == CIGAR) {
            auto [e, l, cig] = edit_cigar_prefix(P, m, match, match.size(), k);
            if (e > k || dominated(e)) continue;
            err = e; len = l;
            cigars.emplace_back(std::move(cig));
        } else {
            auto [e, l] = edit_dist_anchored_prefix(mat, match, m, k, first_matches);
            if (e > k || dominated(e)) continue;
            err = e; len = l;
        }

        node.ctx.set_errors(err);
        node.set_match_len(len);
        if constexpr (mode == CIGAR) ctxts_sorted.emplace_back(&node, uint32_t(cigars.size() - 1));
        else                         ctxts_sorted.emplace_back(&node);
    }
}

// ############################# REPORTING #############################

template <typename idx_t>
template <cigar_mode_t mode, typename report_fnc_t>
void apm_edit<idx_t>::report_contexts(std::vector<context_ref_t<mode>>& ctxts_sorted, report_fnc_t& report) const
{
    // sort the contexts by their suffix array intervals (see apm_context_order_t)
    std::vector<apm_context_order_t> order;
    order.reserve(ctxts_sorted.size());

    for (uint32_t i = 0; i < ctxts_sorted.size(); i++) {
        const node_t& node = node_of<mode>(ctxts_sorted[i]);
        auto [beg, end] = node.ctx.forward_sa_interval();
        order.emplace_back(apm_context_order_t { beg, end, node.len, i });
    }

    std::sort(order.begin(), order.end(), std::less<>());
    std::vector<context_ref_t<mode>> reordered;
    reordered.reserve(order.size());
    for (const apm_context_order_t& o : order) reordered.emplace_back(ctxts_sorted[o.idx]);
    ctxts_sorted.swap(reordered);

    auto ivs = effective_intervals<mode>(ctxts_sorted);
    uint64_t iv_idx = 0;

    // locate the occurrences of each context that is not nested in the previous one, and report each occurrence with
    // the interval that contains its suffix array position; the positions are reported from the center outwards
    for (auto it = ctxts_sorted.begin(); it != ctxts_sorted.end();) {
        const node_t& node = node_of<mode>(*it);
        auto [beg, end] = node.ctx.forward_sa_interval();
        auto loc_ctx = locate_phase_with_positions(node.ctx);
        uint64_t center = loc_ctx.compute_center();
        uint64_t iv_idx_max = ivs.size() - 1;
        uint64_t iv_idx_center = exp_search_max_leq<RIGHT>(center, iv_idx, iv_idx_max, [&](uint64_t x) {return ivs[x].beg();});
        iv_idx = iv_idx_center;

        loc_ctx.locate([&](uint64_t pos, uint64_t occ){
            if (pos <= center) {
                if (iv_idx > 0 && pos < ivs[iv_idx].beg()) [[unlikely]] iv_idx--;
            } else {
                if (iv_idx < iv_idx_center) [[unlikely]] iv_idx = iv_idx_center;
                if (ivs[iv_idx].end() < pos && iv_idx < iv_idx_max) [[unlikely]] iv_idx++;
            }

            if constexpr (mode == CIGAR) report({occ, ivs[iv_idx].len(), ivs[iv_idx].err(), ivs[iv_idx].cigar_idx});
            else                         report({occ, ivs[iv_idx].len(), ivs[iv_idx].err()});
        });

        iv_idx = std::max<uint64_t>(iv_idx, iv_idx_center);

        // skip the contexts nested in this one (their positions are reported with this one)
        do {it++;} while (it != ctxts_sorted.end() && std::get<0>(node_of<mode>(*it).ctx.forward_sa_interval()) <= end);
    }
}

template <typename idx_t>
template <cigar_mode_t mode>
auto apm_edit<idx_t>::effective_intervals(const std::vector<context_ref_t<mode>>& ctxts_sorted)
    -> std::vector<occurrence_interval_t>
{
    std::vector<occurrence_interval_t> ivs;
    ivs.reserve(ctxts_sorted.size());

    /* the suffix array intervals are nested or disjoint: in a sweep over their starting and end positions, the open
       intervals are kept on a stack, the innermost one on top; at each such position, the segment since the previous
       position becomes an interval with the best open context (the fewest errors, the longest on ties) */
    struct open_t {uint64_t end_excl; uint64_t len; uint64_t err; uint32_t cigar_idx;};
    std::vector<open_t> open;
    uint64_t seg_beg = 0;

    for (uint64_t i = 0; i < ctxts_sorted.size() || !open.empty();) {
        uint64_t next_start = i < ctxts_sorted.size()
            ? std::get<0>(node_of<mode>(ctxts_sorted[i]).ctx.forward_sa_interval()) : uint64_t(-1);
        uint64_t next_end = open.empty() ? uint64_t(-1) : open.back().end_excl;
        uint64_t pos = std::min(next_start, next_end);

        if (!open.empty() && seg_beg < pos) {
            const open_t* best = &open.front();

            for (const open_t& o : open)
                if (o.err < best->err || (o.err == best->err && o.len > best->len)) best = &o;

            ivs.emplace_back(occurrence_interval_t{seg_beg, pos - 1, uint32_t(best->len), uint32_t(best->err), best->cigar_idx});
        }

        while (!open.empty() && open.back().end_excl == pos) open.pop_back();
        while (i < ctxts_sorted.size() && std::get<0>(node_of<mode>(ctxts_sorted[i]).ctx.forward_sa_interval()) == pos) {
            const node_t& node = node_of<mode>(ctxts_sorted[i]);
            uint64_t end = std::get<1>(node.ctx.forward_sa_interval());
            open.emplace_back(open_t{end + 1, node.match_len(), node.ctx.errors(), cigar_idx_of<mode>(ctxts_sorted[i])});
            i++;
        }

        seg_beg = pos;
    }

    return ivs;
}

// ############################# SEARCH #############################

template <typename idx_t>
template <typename push_fnc_t>
bool apm_edit<idx_t>::extend_node(node_t& node, direction_t dir, push_fnc_t& push) const
{
    apm_dfs_children<search_context_t<LOCATE>> children;
    auto push_child = [&](search_context_t<LOCATE>&& ctx) { push(child_of(node, ctx)); };
    auto ext_ctx = node.ctx.extend_phase();
    ext_ctx.template prepare_extend_all<NO_TERM>(dir);

    while (ext_ctx.can_extend()) {
        if (next_is_separator(ext_ctx)) { ext_ctx.extend_next(); continue; }
        children.add(ext_ctx.extend_next(), push_child);
    }

    if (!children.has_one()) return false;
    node.ctx = std::move(children.only());
    descend(node);
    return true;
}

template <typename idx_t>
template <typename matrix_word_t>
void apm_edit<idx_t>::set_matrix_input(edit_distance_matrix<matrix_word_t>& matrix,
    const directional_substring<inp_t>& part) const
{
    if (_params.costs == nullptr) {
        matrix.set_input(part, _index.sigma, _index.forward_index().map_int());
        return;
    }

    matrix.set_input(part, _index.sigma, _index.forward_index().map_int(), [&](uint64_t pos, auto report_rank) {
        for (uint8_t c : _params.costs->free_symbols(pos)) {
            uint8_t r = _index.forward_index().map_int()[c];
            if (r != 0) report_rank(r);
        }
    });
}

template <typename idx_t>
template <typename matrix_word_t>
auto apm_edit<idx_t>::search(const inp_t& P, const search_scheme_t& scheme) const -> node_set_t
{
    uint64_t m = P.size();
    uint64_t p = scheme.p;
    node_set_t result;
    const node_t empty_node { search_context_t<LOCATE>(_index) };
    if (m == 0) { result.insert(empty_node.key(), empty_node); return result; }
    if (!fits_search_scheme(m, scheme)) throw std::invalid_argument(too_short_for_scheme("the pattern", m, scheme));

    _piece_owner = this; _matched_pieces.clear(); _free_pieces.clear();
    _max_diagonal_shift = scheme.k;

    std::vector<apm_search_steps> searches;
    std::vector<directional_substring<inp_t>> parts;
    searches.reserve(scheme.S.size());
    parts.reserve(p);

    // partition P into the p parts of the scheme or of its mirror (see apm_partition::mirrored) and compute the search
    // context of each part
    apm_partition<idx_t, LOCATE> partition(_index, P, scheme, EDIT_DISTANCE, _params);
    const search_scheme_t mirror = partition.mirrored() ? mirrored_scheme(scheme) : search_scheme_t();
    const search_scheme_t& oriented = partition.mirrored() ? mirror : scheme;

    for (uint64_t part = 0; part < p; part++) {
        parts.emplace_back(directional_substring<inp_t>(P, partition.beg(part), partition.end(part) - 1));
    }

    for (uint64_t s_idx = 0; s_idx < oriented.S.size(); s_idx++) {
        searches.emplace_back(apm_search_steps(oriented, s_idx));
    }

    std::vector<std::vector<node_t>> stacks;
    stacks.resize(p);

    std::vector<edit_distance_matrix<matrix_word_t>> matrices_fwd, matrices_bwd;
    matrices_fwd.resize(p);
    matrices_bwd.resize(p);

    for (const apm_search_steps& steps : searches) {
        search_data_t<matrix_word_t> data{steps, parts, stacks, matrices_fwd, matrices_bwd, result};

        for (uint64_t p_idx = 0; p_idx < p; p_idx++) {
            parts[steps.part(p_idx)].set_direction(steps.dir(p_idx));
        }

        node_t node = empty_node;
        uint64_t p_idx = 0;

        // the first steps without errors match their parts exactly
        if (steps.k_max(0) == 0) {
            bool match = true;

            for (; match && p_idx < parts.size(); p_idx++) {
                if (steps.k_max(p_idx) > 0) break;
                if (!partition.has_context(steps.part(p_idx))) break; // a part with free symbols is matched by its matrix
                const auto& part = parts[steps.part(p_idx)];

                // the first part starts at its search context (see apm_params_t::start_at_part_contexts)
                if (p_idx == 0 && _params.start_at_part_contexts) {
                    if (!partition.occurs(steps.part(0))) { match = false; break; }
                    node.ctx = partition.context(steps.part(0));
                    node.set_match_len(node.match_len() + part.size());
                    for (uint64_t i = 0; i < part.size(); i++) append_pattern_symbol(node, part.pos(i), part.direction());
                    continue;
                }

                for (uint64_t i = 0; match && i < part.size(); i++) {
                    auto [next_ctx, occurs] = node.ctx.extend(part[i], part.direction());
                    if (occurs) {
                        node.ctx = next_ctx;
                        node.set_match_len(node.match_len() + 1);
                        append_pattern_symbol(node, part.pos(i), part.direction());
                    }
                    match = occurs;
                }
            }

            if (!match) continue;
        }

        // the steps without errors matched every part (all upper bounds of the search are 0)
        if (p_idx == p) {
            add_result(result, with_committed_copy(node));
            continue;
        }

        search_part<matrix_word_t>(node, p_idx, data, {}, {}, {}, {});
    }

    return result;
}

template <typename idx_t>
template <typename matrix_word_t>
void apm_edit<idx_t>::search_part(node_t& node, uint64_t p_idx, search_data_t<matrix_word_t>& data,
    const std::vector<node_t>& descendants_prev_dir, const std::vector<node_t>& descendants_other_dir,
    const std::vector<uint16_t>& dists_prev_dir, const std::vector<uint16_t>& dists_other_dir) const
{
    auto& [steps, parts, stacks, matrices_fwd, matrices_bwd, result] = data;

    const auto& part = parts[steps.part(p_idx)];
    uint64_t k_max = steps.k_max(p_idx);
    direction_t dir = steps.dir(p_idx);
    bool dir_switch = steps.switches_dir(p_idx);
    auto& stack = stacks[p_idx];
    auto& matrix = dir == RIGHT ? matrices_fwd[steps.part(p_idx)] : matrices_bwd[steps.part(p_idx)];

    // the descendants of the direction of this step, and those of the other direction
    const auto& dists = dir_switch ? dists_other_dir : dists_prev_dir;
    const auto& descendants = dir_switch ? descendants_other_dir : descendants_prev_dir;
    const auto& dists_rev = dir_switch ? dists_prev_dir : dists_other_dir;
    const auto& descendants_rev = dir_switch ? descendants_prev_dir : descendants_other_dir;

    // the first column of the matrix: the error count of node, or the distances of the descendants increased by the
    // errors added since they were computed
    std::vector<uint16_t> init_dists;

    if (dists.empty()) {
        init_dists.emplace_back(node.ctx.errors());
    } else {
        uint16_t prev_dist = dir_switch ? *min_element(dists.begin(), dists.end()) : dists[0];
        uint16_t added_errors = node.ctx.errors() - prev_dist;
        init_dists.reserve(dists.size());

        for (uint64_t i = 0; i < dists.size(); i++) {
            init_dists.emplace_back(dists[i] + added_errors);
        }
    }

    if (!matrix.is_initialized()) set_matrix_input(matrix, part);
    matrix.init(k_max, init_dists);
    cluster_t cluster(matrix.last_col_size(), k_max, node.match_len());

    if (matrix.is_in_final_column(0)) {
        auto node_copy = node;
        node_copy.set_match_len(0);
        cluster.set_cell(0, with_committed_copy(node_copy), matrix(0, part.size()));
    }

    // the row of node at the start of the depth-first search (the row of the last descendant, if there is one)
    uint64_t node_row = 0;

    // the descendants are the first rows of the matrix
    if (!descendants.empty()) {
        for (uint64_t i = 0; i < descendants.size() && descendants[i].match_len() < matrix.num_rows(); i++) {
            if (branch_and_bound<matrix_word_t>(matrix, cluster, descendants[i], p_idx, data,
                dists_rev, descendants_rev, {descendants.begin() + i + 1, descendants.end()})
            ) return;
        }

        if (descendants.back().match_len() == uint64_t(matrix.num_rows() - 1)) return;
        if (!dir_switch) node = descendants.back();
        node_row = descendants.back().match_len();
    }

    // the node is in the row node_row; without an extension, the part ends there
    node.set_match_len(node_row);
    auto push = [&](node_t&& child) { stack.emplace_back(std::move(child)); };
    bool has_node = extend_node(node, dir, push);

    if (!has_node && stack.empty() && matrix.is_in_final_column(0)) {
        flush_cluster<matrix_word_t>(cluster, p_idx, data, dists_rev, descendants_rev, {});
        return;
    }

    /* depth-first search over the extensions: compute the row of each node and extend the nodes that are not bounded;
       a node with one extension is continued with it in place (see extend_node) */
    while (has_node || !stack.empty()) {
        if (!has_node) {
            node = std::move(stack.back());
            stack.pop_back();
        }

        append_last_symbol(node, part, dir);
        has_node = false;
        if (branch_and_bound<matrix_word_t>(matrix, cluster, node, p_idx, data, dists_rev, descendants_rev, {})) continue;
        has_node = extend_node(node, dir, push);
    }
}

template <typename idx_t>
void apm_edit<idx_t>::add_result(node_set_t& result, const node_t& node)
{
    const string_key_t key = node.key();
    node_t* node_old = result.find(key);

    if (node_old == nullptr) {
        result.insert(key, node);
        return;
    }

    uint16_t cur_err = node.ctx.errors();
    uint16_t old_err = node_old->ctx.errors();

    if (cur_err  < old_err ||
       (cur_err == old_err && node_old->ctx.s.d() == RIGHT && node.ctx.s.d() == LEFT)
    ) {
        *node_old = node;
    }
}

template <typename idx_t>
template <typename matrix_word_t>
void apm_edit<idx_t>::flush_cluster(cluster_t& cluster, uint64_t p_idx, search_data_t<matrix_word_t>& data,
    const std::vector<uint16_t>& dists_rev, const std::vector<node_t>& descendants_rev,
    const std::vector<node_t>& remaining_descendants) const
{
    auto& [steps, parts, stacks, matrices_fwd, matrices_bwd, result] = data;
    const uint16_t k_min = steps.k_min(p_idx);
    const direction_t dir = steps.dir(p_idx);
    auto at_separator = [&](const node_t& node) { return this->extends_with_separator(node.ctx, dir); };

    // adds the local minima with at least k_min errors to the result (see add_result)
    auto report_at_end = [&]() {
        auto minima = cluster.local_minima(at_separator);
        for (auto& node_cur : minima) if (node_cur.ctx.errors() >= k_min) add_result(result, node_cur);
    };

    // after the last step, at an outer part or before a direction switch, continue from each local minimum
    bool next_switches = (p_idx + 1 < parts.size()) && steps.switches_dir(p_idx + 1);

    if (steps.is_outer_part(p_idx) || next_switches) {
        if (p_idx + 1 == parts.size()) {
            report_at_end();
        } else {
            auto minima = cluster.local_minima(at_separator);

            for (auto& node : minima)
                if (node.ctx.errors() >= k_min)
                    search_part<matrix_word_t>(node, p_idx + 1, data, {}, descendants_rev, {}, dists_rev);
        }

        return;
    }

    // else continue in the same direction from one cell, with the cells after it as descendants
    std::vector<node_t> descendants;
    std::vector<uint16_t> dists;
    auto [center, found] = cluster.choose_center(k_min, descendants, dists, at_separator);

    // no local minimum with at least k_min errors: the alignments on this path are found by another search
    if (!found) return;

    descendants.insert(descendants.end(), remaining_descendants.begin(), remaining_descendants.end());

    for (uint64_t i = 0; i < descendants.size(); i++) {
        descendants[i].set_match_len(i + 1);
    }

    // remove the distances above the upper bound of the next step from the end of the first column of the next matrix
    uint64_t k_max = steps.k_max(p_idx + 1);
    while (dists.back() > k_max) dists.pop_back();

    search_part<matrix_word_t>(center, p_idx + 1, data, descendants, descendants_rev, dists, dists_rev);
}

template <typename idx_t>
template <typename matrix_word_t>
bool apm_edit<idx_t>::branch_and_bound(edit_distance_matrix<matrix_word_t>& matrix, cluster_t& cluster,
    const node_t& node, uint64_t p_idx, search_data_t<matrix_word_t>& data, const std::vector<uint16_t>& dists_rev,
    const std::vector<node_t>& descendants_rev, const std::vector<node_t>& remaining_descendants) const
{
    auto& [steps, parts, stacks, matrices_fwd, matrices_bwd, result] = data;
    uint64_t row = node.match_len();

    if (row >= matrix.num_rows()) return true;
    bool alive = matrix.compute_row(row, node.ctx.last_added_symbol());

    if (matrix.is_in_final_column(row)) {
        uint64_t cluster_idx = row - (matrix.num_rows() - cluster.size());
        cluster.set_cell(cluster_idx, with_committed_copy(node), matrix(row, matrix.num_cols() - 1));

        // the cell is terminal if no cell of its row has at most k_max errors, if only deletions remain, or if the
        // context cannot be extended (its only extension can be the separator)
        bool terminal = !alive || matrix.only_vertical_gaps_left(row);

        if (!terminal) {
            node_t node_ext = node;
            auto ext_ctx = node_ext.ctx.extend_phase();
            ext_ctx.template prepare_extend_all<NO_TERM>(steps.dir(p_idx));
            if (ext_ctx.can_extend() && next_is_separator(ext_ctx)) ext_ctx.extend_next();
            terminal = !ext_ctx.can_extend();
        }

        if (terminal) {
            flush_cluster<matrix_word_t>(cluster, p_idx, data, dists_rev, descendants_rev, remaining_descendants);
            return true;
        }
    }

    return !alive;
}

// ############################# QUERY METHODS #############################

template <typename idx_t>
template <cigar_mode_t mode, typename report_fnc_t>
std::vector<cigar_t<uint8_t>> apm_edit<idx_t>::locate(const inp_t& P, const search_scheme_t& scheme,
    report_fnc_t report) const
{
    uint64_t m = P.size();
    require_storable_occurrence_length(m + scheme.k);
    _piece_owner = this;
    if (_params.costs != nullptr) return locate_model<mode>(P, scheme, report);
    if (scheme.k == 0) return this->template locate_exact<mode>(P, report);

    std::vector<cigar_t<uint8_t>> cigars;
    auto ctxts = search(P, scheme);
    std::vector<context_ref_t<mode>> ctxts_sorted;
    ctxts_sorted.reserve(ctxts.size());
    if constexpr (mode == CIGAR) cigars.reserve(ctxts.size());
    with_edit_distance_word(scheme.k, [&]<typename word_t>() {
        recompute_errors<mode, word_t>(ctxts, P, m, scheme.k, ctxts_sorted, cigars);
    });

    // the occurrences at the text ends (the search misses some of them), at most one per position (the one with fewer
    // errors)
    std::vector<aprx_occ_t<mode>> end_occ;

    if constexpr (can_decode_text_ends) {
        auto collect = [&](const aprx_occ_t<mode>& o) {
            if (!end_occ.empty() && end_occ.back().pos() == o.pos()) {
                if (o.err() < end_occ.back().err()) end_occ.back() = o;
            } else {
                end_occ.emplace_back(o);
            }
        };

        locate_at_text_ends<mode>(P, m, scheme.k, collect, cigars);
    }

    // an occurrence in end_occ is replaced by an occurrence of the search at the same position with fewer errors
    auto report_unless_at_text_end = [&](const aprx_occ_t<mode>& o) {
        for (aprx_occ_t<mode>& e : end_occ) {
            if (e.pos() == o.pos()) {
                if (o.err() < e.err()) e = o;
                return;
            }
        }

        report(o);
    };

    report_contexts<mode>(ctxts_sorted, report_unless_at_text_end);
    for (const aprx_occ_t<mode>& e : end_occ) report(e);
    return cigars;
}

#endif // MOVE_R_ALGORITHMS_MOVE_RB_APM_EDIT_TPP
