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
 * @file edit_distance.tpp
 * @brief edit-distance approximate pattern matching on a bidirectional move-structure index, via search schemes.
 *
 * The pattern is split into the parts of a search scheme; each scheme search runs a banded bit-parallel
 * edit-distance-matrix DFS over the index extensions of each part, carrying a part's final-column band forward to
 * seed the next part's matrix. No in-text verification is performed.
 *
 * Techniques: search schemes over bidirectional FM-indexes (Lam et al.; Kucherov, Salikhov, Tsur; Kianfar et al.);
 * branch-and-bound over a banded bit-parallel matrix (Columba; Renders, Depuydt, Fostier).
 */

#pragma once

#include <algorithm>
#include <variant>
#include <vector>

#include <ips4o.hpp>
#include <gtl/phmap.hpp>
#include <misc/apm.hpp>

/**
 * @brief precomputes, for a single search of a search scheme, the per-step information needed to execute it
 *        (search directions, direction switches and the range of previously processed parts)
 */
class edit_dist_search
{
  protected:
    const search_scheme_t& scheme;
    const std::vector<search_step_t>& search_arr;
    uint8_t search_idx;
    std::vector<direction_t> dirs;
    std::vector<bool> is_dir_switch;
    std::vector<uint8_t> leftmost_prev_part;
    std::vector<uint8_t> rightmost_prev_part;

  public:
    /**
     * @brief precomputes the per-step information for the search with index search_idx in scheme
     * @param scheme the search scheme
     * @param search_idx the index of the search within the scheme
     */
    edit_dist_search(const search_scheme_t& scheme, uint8_t search_idx)
        : scheme(scheme), search_arr(scheme.S[search_idx]), search_idx(search_idx)
    {
        dirs.reserve(scheme.p);
        dirs.emplace_back((search_arr[1].part > search_arr[0].part) ? RIGHT : LEFT);

        for (uint8_t i = 1; i < scheme.p; i++) {
            dirs.emplace_back((search_arr[i].part > search_arr[i - 1].part) ? RIGHT : LEFT);
        }

        is_dir_switch.reserve(scheme.p);
        is_dir_switch.emplace_back(false);

        for (uint8_t i = 1; i < dirs.size(); i++) {
            is_dir_switch.emplace_back(dirs[i] != dirs[i - 1]);
        }

        leftmost_prev_part.reserve(scheme.p);
        rightmost_prev_part.reserve(scheme.p);
        leftmost_prev_part.emplace_back(search_arr[0].part);
        rightmost_prev_part.emplace_back(search_arr[0].part);

        for (uint8_t i = 1; i < scheme.p; i++) {
            uint8_t cur_part = search_arr[i].part;

            if (cur_part < leftmost_prev_part[i - 1]) {
                leftmost_prev_part.emplace_back(cur_part);
                rightmost_prev_part.emplace_back(rightmost_prev_part[i - 1]);
            } else {
                leftmost_prev_part.emplace_back(leftmost_prev_part[i - 1]);
                rightmost_prev_part.emplace_back(cur_part);
            }
        }
    }

    /**
     * @brief returns the minimum number of errors allowed after the i-th search step
     * @param i index of the search step
     * @return the lower error bound of the i-th search step
     */
    uint8_t lower_bound(uint8_t i) const
    {
        return search_arr[i].k_min;
    }

    /**
     * @brief returns the maximum number of errors allowed after the i-th search step
     * @param i index of the search step
     * @return the upper error bound of the i-th search step
     */
    uint8_t upper_bound(uint8_t i) const
    {
        return search_arr[i].k_max;
    }

    /**
     * @brief returns the pattern part processed in the i-th search step
     * @param i index of the search step
     * @return the index of the pattern part processed in the i-th search step
     */
    uint8_t part(uint8_t i) const
    {
        return search_arr[i].part;
    }

    /**
     * @brief returns the leftmost pattern part processed before the idx-th search step
     * @param idx index of the search step
     * @return the index of the leftmost previously processed part
     */
    uint8_t leftmost_previous_part(uint8_t idx) const
    {
        return leftmost_prev_part[idx - 1];
    }

    /**
     * @brief returns the rightmost pattern part processed before the idx-th search step
     * @param idx index of the search step
     * @return the index of the rightmost previously processed part
     */
    uint8_t rightmost_previous_part(uint8_t idx) const
    {
        return rightmost_prev_part[idx - 1];
    }

    /**
     * @brief returns the direction in which the i-th search step extends the pattern
     * @param i index of the search step
     * @return the search direction of the i-th search step
     */
    direction_t part_dir(uint8_t i) const
    {
        return dirs[i];
    }

    /**
     * @brief returns whether the i-th search step switches the search direction
     * @param i index of the search step
     * @return whether the i-th search step switches the search direction
     */
    bool does_part_switch_dir(uint8_t i) const
    {
        return is_dir_switch[i];
    }

    /**
     * @brief returns the number of pattern parts
     * @return the number of pattern parts
     */
    uint8_t num_parts() const
    {
        return scheme.p;
    }

    /**
     * @brief returns the index of this search within its search scheme
     * @return the index of this search within its search scheme
     */
    uint8_t search_index() const
    {
        return search_idx;
    }

    /**
     * @brief returns whether the part processed in the i-th search step is an outer (edge) part of the pattern
     * @param i index of the search step
     * @return whether the i-th part is the first or the last pattern part
     */
    bool is_edge(uint8_t i) const
    {
        return search_arr[i].part == 0 || search_arr[i].part == scheme.p - 1;
    }
};

/**
 * @brief approximate pattern matching (APM) helper for a bidirectional index; holds a reference to the index
 *        and implements the edit-distance locate algorithm outside the index class
 * @tparam idx_t the bidirectional index type (move_rb or a compatible index, e.g. b_move_adapter); it must
 *         expose pos_type, inp_t, supports_locate, search_context_t<>, empty_context<>(), forward_index(),
 *         sigma and n
 * @tparam mode whether the search maintains each context's matched string (mode == CIGAR) so that locate can
 *         additionally report a CIGAR alignment (and the exact ed(P, matched string)) per occurrence
 */
template <typename idx_t, cigar_mode_t mode = NO_CIGAR>
class apm_edit {
protected:
    using pos_t = typename idx_t::pos_type;
    const idx_t& index; // the index to perform approximate pattern matching on

public:
    /**
     * @brief constructs an approximate-pattern-matching helper for the index
     * @param index the index to perform approximate pattern matching on
     */
    apm_edit(const idx_t& index) : index(index) {}

protected:
    using inp_t = typename idx_t::inp_t;
    using sym_t = typename idx_t::sym_type;

    // whether the index can decode both text ends (bidirectional revert_range) for the boundary-clipped occurrence pass;
    // move_rb can, diagnostic adapters (columba, b-move) cannot and simply skip it. The designated initializer selects
    // the retrieve_params overload unambiguously (a bare {} would also match the std::string / file-name overload)
    static constexpr bool can_decode_boundary = requires(const idx_t& i) {
        i.backward_index().revert_range({.l = pos_t(0), .r = pos_t(0)});
        i.forward_index().revert_range({.l = pos_t(0), .r = pos_t(0)});
    };

    template <query_support_t query_support>
    using search_context_t = typename idx_t::template search_context_t<query_support>;

    struct match_cell_t {
        pos_t prev;      // preceding cell index in match_pool, or no_match at the root
        pos_t pos;       // run: start index of the pattern substring P[pos, pos+len); unused for a literal
        pos_t len;       // run: length (>= 1); 0 marks a literal (whose symbol is sym)
        sym_t sym;       // literal: the matched text symbol; unused for a run
        direction_t dir; // extension direction: LEFT (prepended) or RIGHT (appended)
        uint16_t refs;   // reference count: cells pointing here via prev, plus contexts whose str head is here
    };

    // sentinel index marking the empty (root) matched-string list
    static constexpr pos_t no_match = pos_t(-1);

    // the apm_edit whose match_pool owns the cells referenced by the str_ref handles alive on this thread. Set for the
    // duration of a locate (a single query is single-threaded); the handles use it to reach inc_ref/dec_ref.
    static inline thread_local const apm_edit* active_pool = nullptr;

    // a reference-counted handle to a match_pool cell (or no_match). Copy increments, move steals, destroy decrements
    // the cell's refcount; a cell reaching zero references is freed (its slot recycled) and its prev released. This
    // ensures match_pool does not store discarded cells.
    struct str_ref {
        pos_t idx = no_match;

        str_ref() = default;
        explicit str_ref(pos_t i) : idx(i) { if (idx != no_match) active_pool->inc_ref(idx); }
        str_ref(const str_ref& o) : idx(o.idx) { if (idx != no_match) active_pool->inc_ref(idx); }
        str_ref(str_ref&& o) noexcept : idx(o.idx) { o.idx = no_match; }

        str_ref& operator=(const str_ref& o)
        {
            if (this != &o) {
                if (o.idx != no_match) active_pool->inc_ref(o.idx);
                if (idx != no_match) active_pool->dec_ref(idx);
                idx = o.idx;
            }
            return *this;
        }

        str_ref& operator=(str_ref&& o) noexcept
        {
            if (this != &o) {
                if (idx != no_match) active_pool->dec_ref(idx);
                idx = o.idx;
                o.idx = no_match;
            }
            return *this;
        }

        ~str_ref() { if (idx != no_match) active_pool->dec_ref(idx); }
    };

    struct node_t {
        search_context_t<LOCATE> ctx;
        pos_t len = 0;
        str_ref str;
        pos_t run_pos = 0;
        pos_t run_len = 0;
        direction_t run_dir = RIGHT;
        bool rep = false;

        node_t() = default;
        node_t(const search_context_t<LOCATE>& c) : ctx(c) {}

        pos_t match_len() const { return len; }
        void set_match_len(pos_t d) { len = d; }
        bool reported() const { return rep; }
        void set_reported(bool r) { rep = r; }
        pos_t match_str() const { return str.idx; }
        void set_match_str(pos_t m) { str = str_ref(m); }

        bool operator==(const node_t& o) const {
            auto [b, e] = ctx.forward_sa_interval();
            auto [ob, oe] = o.ctx.forward_sa_interval();
            return b == ob && e == oe && len == o.len;
        }

        bool operator<(const node_t& o) const {
            auto [b, e] = ctx.forward_sa_interval();
            auto [ob, oe] = o.ctx.forward_sa_interval();
            if (b != ob) return b < ob;
            if (e != oe) return e > oe;
            return len < o.len;
        }

        struct hash {
            pos_t operator()(const node_t& n) const {
                auto [b, e] = n.ctx.forward_sa_interval();
                pos_t h = pos_hash<pos_t>(b);
                hash_combine<pos_t>(h, pos_hash<pos_t>(e));
                hash_combine<pos_t>(h, pos_hash<pos_t>(n.len));
                return h;
            }
        };
    };

    using node_set_t = tsl::sparse_set<node_t, typename node_t::hash>;
    mutable std::vector<match_cell_t> match_pool;
    mutable std::vector<pos_t> match_free;
    mutable pos_t band_k = 0;

    // allocates a cell slot (reusing a freed one when available)
    pos_t alloc_cell() const
    {
        if (!match_free.empty()) { pos_t i = match_free.back(); match_free.pop_back(); return i; }
        match_pool.emplace_back();
        return match_pool.size() - 1;
    }

    // increments a cell's reference count (saturating, so an extremely shared hub cell is pinned, never wrapping)
    void inc_ref(pos_t idx) const
    {
        if (match_pool[idx].refs != UINT16_MAX) match_pool[idx].refs++;
    }

    // decrements a cell's reference count; on reaching zero the cell is freed and its prev released in turn
    // (iteratively -- chains can be long). A saturated (pinned) cell is never freed.
    void dec_ref(pos_t idx) const
    {
        while (idx != no_match) {
            uint16_t& r = match_pool[idx].refs;
            if (r == UINT16_MAX) break; // pinned
            if (r > 1) { r--; break; }  // still referenced
            r = 0;                      // last reference: recycle the slot and release its prev
            pos_t prev = match_pool[idx].prev;
            match_free.emplace_back(idx);
            idx = prev;
        }
    }

    // appends a literal matched symbol to the chain prev (dir == LEFT <=> prepended), returning the new list head
    pos_t push_literal(pos_t prev, sym_t sym, direction_t dir) const
    {
        pos_t i = alloc_cell();
        match_pool[i] = match_cell_t{prev, 0, 0, sym, dir, 0};
        if (prev != no_match) inc_ref(prev); // the new cell references its predecessor
        return i;
    }

    // appends a run (a reference to P[pos, pos+len)) to the chain prev (dir == LEFT <=> prepended), returning the head
    pos_t push_run(pos_t prev, pos_t pos, pos_t len, direction_t dir) const
    {
        pos_t i = alloc_cell();
        match_pool[i] = match_cell_t{prev, pos, len, sym_t{}, dir, 0};
        if (prev != no_match) inc_ref(prev);
        return i;
    }

    // a child inherits its parent's committed matched-string chain and pending run (its own symbol is recorded when
    // the child is popped from the DFS stack, so fanned-out children that get pruned never create a pool cell)
    void inherit_match(node_t& child, const node_t& parent) const
    {
        child.set_match_str(parent.match_str());
        child.run_pos = parent.run_pos;
        child.run_len = parent.run_len;
        child.run_dir = parent.run_dir;
    }

    // commits a node's pending run (if any) into its matched-string chain, resetting it
    void flush_run(node_t& node) const
    {
        if (node.run_len > 0) {
            node.set_match_str(push_run(node.match_str(), node.run_pos, node.run_len, node.run_dir));
            node.run_len = 0;
        }
    }

    // returns a copy of node whose pending run has been committed (so its matched string is complete for a cell that
    // is about to be stored in a cluster and later reported/carried)
    node_t flushed(node_t node) const
    {
        flush_run(node);
        return node;
    }

    // records a matched symbol known to equal P[pp] (dir == LEFT <=> prepended): extends the pending run if pp
    // continues it contiguously, otherwise flushes the old run and starts a new one at pp
    void run_match(node_t& node, pos_t pp, direction_t dir) const
    {
        bool contig = node.run_len > 0 && node.run_dir == dir &&
            (dir == RIGHT ? pp == node.run_pos + node.run_len : pp == node.run_pos - 1);
        if (contig) { if (dir == LEFT) node.run_pos = pp; node.run_len++; }
        else { flush_run(node); node.run_pos = pp; node.run_len = 1; node.run_dir = dir; }
    }

    // records the node's last-added matched symbol (a new extension of part `part` in direction `dir`): if it equals
    // the pattern at its diagonal position, extend the pending run contiguously (or start a new one); otherwise flush
    // the run and append the symbol as a literal. Correctness is independent of the diagonal position's accuracy: a
    // run is only formed when P[pos] equals the symbol, so reconstruct is exact; a stale position just yields a literal.
    void run_update(node_t& node, const directional_substring<pos_t, inp_t>& part, direction_t dir) const
    {
        sym_t s = node.ctx.last_added_symbol();
        const inp_t& P = part.underlying();
        const pos_t m = P.size();

        if (node.run_len > 0 && node.run_dir == dir) {
            pos_t cont = dir == RIGHT ? node.run_pos + node.run_len : node.run_pos - 1;
            if (cont < m && P[cont] == s) { if (dir == LEFT) node.run_pos = cont; node.run_len++; return; }
        }

        const int64_t pp = int64_t(part.pos(node.match_len() - 1));
        auto match_at = [&](int64_t x) { return x >= 0 && x < int64_t(m) && P[x] == s; };
        pos_t j = no_match;
        if (match_at(pp)) j = pos_t(pp);
        else for (int64_t d = 1; d <= int64_t(band_k) && j == no_match; d++) {
            if (match_at(pp + d)) j = pos_t(pp + d);
            else if (match_at(pp - d)) j = pos_t(pp - d);
        }

        flush_run(node);
        if (j != no_match) { node.run_pos = j; node.run_len = 1; node.run_dir = dir; }
        else node.set_match_str(push_literal(node.match_str(), s, dir));
    }

    // reconstructs the matched text string (in text order) from a matched-string list head index into out (which is
    // reused across contexts to avoid a per-context allocation), expanding runs into their referenced pattern
    // substrings P[pos, pos+len)
    void reconstruct_match(const inp_t& P, pos_t head, inp_t& out) const
    {
        static thread_local std::vector<pos_t> cells;
        cells.clear();
        pos_t total = 0, n_left = 0;

        for (pos_t i = head; i != no_match; i = match_pool[i].prev) {
            const match_cell_t& c = match_pool[i];
            cells.emplace_back(i);
            pos_t l = c.len > 0 ? c.len : 1;
            total += l;
            if (c.dir == LEFT) n_left += l;
        }

        out.resize(total);
        pos_t lo = n_left, hi = n_left;

        for (auto it = cells.rbegin(); it != cells.rend(); ++it) {
            const match_cell_t& c = match_pool[*it];
            if (c.dir == LEFT) {
                if (c.len > 0) { lo -= c.len; for (pos_t j = 0; j < c.len; j++) out[lo + j] = P[c.pos + j]; }
                else out[--lo] = c.sym;
            } else {
                if (c.len > 0) { for (pos_t j = 0; j < c.len; j++) out[hi + j] = P[c.pos + j]; hi += c.len; }
                else out[hi++] = c.sym;
            }
        }
    }

    using sorted_ctx_t = std::conditional_t<mode == CIGAR, std::pair<const node_t*, pos_t>, const node_t*>;
    static const node_t& ctx_of(const sorted_ctx_t& e) { if constexpr (mode == CIGAR) return *e.first; else return *e; }
    static pos_t cig_of(const sorted_ctx_t& e) { if constexpr (mode == CIGAR) return e.second; else return no_match; }

    // recomputes each context's exact error and best-prefix length by aligning its matched string against P
    template <typename word_t>
    void recompute_errors(node_set_t& ctxts, const inp_t& P, pos_t m, pos_t k,
        std::vector<sorted_ctx_t>& ctxts_sorted, std::vector<cigar_t>& cigars) const
    {
        static thread_local edit_distance_matrix<word_t> mat;
        if constexpr (mode == NO_CIGAR) {
            directional_substring<pos_t, inp_t> P_col(P, 0, m - 1, RIGHT);
            mat.set_input(P_col, index.sigma, index.forward_index().map_int());
        }

        inp_t match;
        for (const node_t& node : ctxts) {
            reconstruct_match(P, node.match_str(), match);

            pos_t err, len;
            if constexpr (mode == CIGAR) {
                auto [e, l, cig] = edit_cigar_prefix<pos_t>(P, m, match, match.size(), k);
                if (e > k) continue;
                err = e; len = l;
                cigars.emplace_back(std::move(cig));
            } else {
                auto [e, l] = edit_dist_prefix(mat, match, m, k);
                if (e > k) continue;
                err = e; len = l;
            }

            node_t& mut = const_cast<node_t&>(node);
            mut.ctx.set_errors(err);
            mut.set_match_len(len);
            if constexpr (mode == CIGAR) ctxts_sorted.emplace_back(&node, cigars.size() - 1);
            else                         ctxts_sorted.emplace_back(&node);
        }
    }

    // references to the structures shared across one search's rec_search recursion
    template <typename matrix_word_t>
    using data_struct_ref_t = std::tuple<
        const edit_dist_search&,                                  // search
        const std::vector<directional_substring<pos_t, inp_t>>&,  // parts
        std::vector<std::vector<node_t>>&,                        // stacks
        std::vector<edit_distance_matrix<matrix_word_t>>&,        // matrices_fwd
        std::vector<edit_distance_matrix<matrix_word_t>>&,        // matrices_bwd
        node_set_t&>;                                             // result

    // an effective interval: every SA-position in [beg, end] is reported as the same occurrence (len, err). In
    // CIGAR mode, str is the winning context's matched-string head (shared by all positions in the interval).
    // idx is the CIGAR index (in the cigars array) of the winning context, in CIGAR mode (unused otherwise)
    struct interval_t {pos_t beg; pos_t end; pos_t err; pos_t len; pos_t idx = no_match;};

    // collects the cells of one part's final matrix column. Each cell is a final-column search context.
    class cluster_t {
      protected:
        pos_t start_depth;
        uint16_t k_max;
        pos_t last_cell;
        std::vector<node_t> nodes;
        std::vector<uint16_t> dists;

      public:
        // last_cell == pos_t(-1) marks an empty cluster; the first set_cell raises it to the written index
        cluster_t(pos_t size, pos_t k_max, pos_t start_depth)
            : start_depth(start_depth), k_max(k_max), last_cell(pos_t(-1)),
              nodes(size), dists(size, k_max + 1) {}

        // records the cell at idx: its context (marked not-yet-reported) and its edit distance
        void set_cell(pos_t idx, node_t node, uint16_t dist)
        {
            nodes[idx] = std::move(node);
            nodes[idx].set_reported(false);
            dists[idx] = dist;
            last_cell = idx;
        }

        pos_t size() const { return dists.size(); }

        // whether cell i is a local minimum of the column (<= both of its in-bounds neighbours)
        bool is_local_min(pos_t i) const
        {
            return (i == 0 || dists[i] <= dists[i - 1]) && (i == last_cell || dists[i] <= dists[i + 1]);
        }

        // marks cell i reported and returns its reportable node: a copy carrying its distance as the error and its
        // match_len rebased by start_depth (the matched length accumulated before this part)
        node_t emit(pos_t i)
        {
            nodes[i].set_reported(true);
            node_t node = nodes[i];
            node.ctx.set_errors(dists[i]);
            node.set_match_len(node.match_len() + start_depth);
            return node;
        }

        // returns the cells to report as occurrences: every in-bounds cell that is a local minimum of the column's
        // distances (no greater than either neighbour, within the error budget, not already reported)
        std::vector<node_t> centers()
        {
            std::vector<node_t> centers;
            centers.reserve(last_cell + 1);

            for (pos_t i = 0; i <= last_cell; i++) {
                if (!nodes[i].reported() && dists[i] <= k_max && is_local_min(i)) {
                    centers.emplace_back(emit(i));
                }
            }

            return centers;
        }

        // selects the cluster center and the next band's seed: the center is the first local-minimum cell within
        // [k_min, k_max]; the deeper cells become its descendants and their distances become the next matrix' first
        // column (init_dists). Returns {center, true} on success, {.., false} otherwise.
        std::tuple<node_t, bool> cluster_center(uint16_t k_min,
            std::vector<node_t>& desc, std::vector<uint16_t>& init_dists
        ) {
            desc.reserve(dists.size());
            init_dists.reserve(dists.size());
            node_t node;
            pos_t i = last_cell + 1;

            for (pos_t j = 0; j <= last_cell; j++) {
                if (dists[j] < k_min || dists[j] > k_max) continue;
                if (is_local_min(j)) { i = j; break; }
            }

            if (i <= last_cell) {
                // carry whenever a valid center exists: emit() rebases a fresh copy (idempotent on the cell's
                // match_len), so the already-reported state must not suppress the carry here -- it only guards
                // against double-emitting an occurrence elsewhere
                node = emit(i);
                bool res = true;

                init_dists.emplace_back(dists[i]);

                for (pos_t j = i + 1; j <= last_cell; j++) {
                    desc.emplace_back(nodes[j]);
                    init_dists.emplace_back(dists[j]);
                }

                for (pos_t k = 1; k < init_dists.size(); k++) {
                    if (init_dists[k] < k_min && init_dists[k] <= init_dists[k - 1] &&
                       (k == init_dists.size() - 1 || init_dists[k] <= init_dists[k + 1])
                    ) {
                        pos_t high = 0;
                        pos_t low = init_dists.size() - 1;

                        for (pos_t l = k; l-- > 0;) {
                            if (init_dists[l] != init_dists[l + 1] + 1) {
                                high = l + 1;
                                break;
                            }
                        }

                        for (pos_t l = k + 1; l < init_dists.size(); l++) {
                            if (init_dists[l] != init_dists[l - 1] + 1) {
                                low = l - 1;
                                break;
                            }
                        }

                        if (high != 0 && low != init_dists.size() - 1) {
                            pos_t left_base = init_dists[high - 1];
                            pos_t right_base = init_dists[low + 1];

                            for (pos_t l = high; l <= low; l++) {
                                pos_t wall_l = left_base + (l - (high - 1));
                                pos_t wall_r = right_base + ((low + 1) - l);
                                init_dists[l] = std::min<pos_t>(k_max + 1, std::min(wall_l, wall_r));
                            }
                        } else if (high == 0 && low != init_dists.size() - 1) {
                            for (pos_t l = low; l-- > 0;) {
                                init_dists[l] = init_dists[l + 1] + 1;
                            }
                        } else if (high != 0 && low == init_dists.size() - 1) {
                            for (pos_t l = high; l < init_dists.size(); l++) {
                                init_dists[l] = init_dists[l - 1] + 1;
                            }
                        }
                    }
                }

                return {node, res};
            }

            return {node, false};
        }
    };

    // The bidirectional search only finds occurrences fitting ENTIRELY within the text -- a part reaching a text end
    // stops there and never covers a pattern edge running off it. Those boundary-clipped occurrences can only touch the
    // first / last (m + k) characters. Both windows are decoded once and P matched directly
    template <typename report_fnc_t>
    void locate_boundary(const inp_t& P, pos_t m, pos_t k, report_fnc_t report) const
    {
        if (k <= edit_distance_matrix<uint64_t>::k_limit) locate_boundary_impl<uint64_t>(P, m, k, report);
        else                                              locate_boundary_impl<__uint128_t>(P, m, k, report);
    }

    template <typename word_t, typename report_fnc_t>
    void locate_boundary_impl(const inp_t& P, pos_t m, pos_t k, report_fnc_t& report) const
    {
        pos_t n = index.n - 1;
        if (m == 0 || n == 0) return;
        pos_t w = std::min<pos_t>(n, m + k);

        pos_t left_hi  = index.has_sequences() ? std::min<pos_t>(index.sequence_length(0), w) : w;
        pos_t right_lo = index.has_sequences()
            ? std::max<pos_t>(index.sequence_start(index.num_sequences() - 1), n - w) : n - w;

        inp_t left = index.backward_index().revert_range({.l = n - left_hi, .r = n - 1});
        inp_t right = index.forward_index().revert_range({.l = right_lo, .r = n - 1});
        std::reverse(left.begin(), left.end());

        static thread_local edit_distance_matrix<word_t> mat;
        
        auto emit = [&](const inp_t& Q, const inp_t& W, auto pos_of, bool rev_cig) {
            if constexpr (mode == CIGAR) {
                auto [err, len, cig] = edit_cigar_prefix<pos_t>(Q, m, W, W.size(), k);
                if (err > k) return;
                if (rev_cig) std::reverse(cig.begin(), cig.end());
                report({pos_of(len), len, err, std::move(cig)});
            } else {
                directional_substring<pos_t, inp_t> Q_col(Q, 0, m - 1, RIGHT);
                mat.set_input(Q_col, index.sigma, index.forward_index().map_int());
                auto [err, len] = edit_dist_prefix(mat, W, m, k);
                if (err > k) return;
                report({pos_of(len), len, err});
            }
        };

        emit(P, left, [](pos_t) { return pos_t(0); }, false);
        inp_t P_rev(P.rbegin(), P.rend()), right_rev(right.rbegin(), right.rend());
        emit(P_rev, right_rev, [n](pos_t len) { return n - len; }, true);
    }

public:
    /**
     * @brief locates a pattern with at most k errors (edit distance). Guarantees the relaxed contract: every
     *        reported (pos, len, err) has some p' in [pos-2k-1, pos+2k+1] and some length L with ed(P, T[p',p'+L)) <= k.
     * @tparam report_fnc_t type of the function report (called with {pos, len, err})
     * @param P the pattern to search
     * @param scheme the search scheme to use (provides k)
     * @param report the occurrence callback
     */
    template <typename report_fnc_t>
    void locate(const inp_t& P, const search_scheme_t& scheme, report_fnc_t report) const requires(idx_t::supports_locate)
    {
        pos_t m = P.size();
        active_pool = this; // owner of the str_ref handles created during this query

        if (scheme.k == 0) {
            if constexpr (mode == CIGAR) {
                cigar_t cig = m > 0 ? cigar_t{MATCH(m)} : cigar_t{};
                index.forward_index().locate(P, [&](pos_t occ){report({occ, m, 0, cig});});
            } else {
                index.forward_index().locate(P, [&](pos_t occ){report({occ, m, 0});});
            }
            return;
        }

        auto ctxts = search(P, scheme);
        std::vector<sorted_ctx_t> ctxts_sorted;
        ctxts_sorted.reserve(ctxts.size());
        std::vector<cigar_t> cigars; // one CIGAR per surviving context; empty and untouched in NO_CIGAR mode
        if constexpr (mode == CIGAR) cigars.reserve(ctxts.size());
        if (scheme.k <= edit_distance_matrix<uint64_t>::k_limit)
             recompute_errors<uint64_t>(ctxts, P, m, scheme.k, ctxts_sorted, cigars);
        else recompute_errors<__uint128_t>(ctxts, P, m, scheme.k, ctxts_sorted, cigars);

        // sort the contexts by SA-interval (required by effective_intervals). In CIGAR mode each element carries its
        // CIGAR's index, so this single in-place sort leaves the cigars array untouched -- the index rides along
        ips4o::sort(ctxts_sorted.begin(), ctxts_sorted.end(),
            [](const sorted_ctx_t& a, const sorted_ctx_t& b){ return ctx_of(a) < ctx_of(b); });

        auto ivs = effective_intervals(ctxts_sorted);
        pos_t iv_idx = 0;

        for (auto it = ctxts_sorted.begin(); it != ctxts_sorted.end();) {
            const node_t& node = ctx_of(*it);
            auto [beg, end] = node.ctx.forward_sa_interval();
            auto loc_ctx = node.ctx.locate_phase();
            pos_t center = loc_ctx.compute_center();
            pos_t iv_idx_max = ivs.size() - 1;
            pos_t iv_idx_center = exp_search_max_leq<pos_t, RIGHT>(center, iv_idx, iv_idx_max, [&](pos_t x) {return ivs[x].beg;});
            iv_idx = iv_idx_center;

            loc_ctx.locate([&](pos_t pos, pos_t occ){
                if (pos <= center) {
                    if (iv_idx > 0 && pos < ivs[iv_idx].beg) [[unlikely]] iv_idx--;
                } else {
                    if (iv_idx < iv_idx_center) [[unlikely]] iv_idx = iv_idx_center;
                    if (ivs[iv_idx].end < pos && iv_idx < iv_idx_max) [[unlikely]] iv_idx++;
                }

                if constexpr (mode == CIGAR) report({occ, ivs[iv_idx].len, ivs[iv_idx].err, cigars[ivs[iv_idx].idx]});
                else                         report({occ, ivs[iv_idx].len, ivs[iv_idx].err});
            });

            iv_idx = std::max<pos_t>(iv_idx, iv_idx_center);
            do {it++;} while (it != ctxts_sorted.end() && std::get<0>(ctx_of(*it).ctx.forward_sa_interval()) <= end);
        }

        // the search above omits occurrences clipped at a text end; recover them from the boundary windows (only for
        // indexes that can decode both text ends -- move_rb; diagnostic adapters without bidirectional revert skip this)
        if constexpr (can_decode_boundary) locate_boundary(P, m, scheme.k, report);
    }

    /**
     * @brief the search phase in isolation: the deduplicated search contexts
     *        (carrying match_len bookkeeping; and, when mode == CIGAR, each context's matched string). Picks the matrix
     *        word type fitting k.
     * @param P the pattern to search
     * @param scheme the search scheme to use
     * @return all relaxed search contexts of P in T with at most k errors
     */
    auto search(const inp_t& P, const search_scheme_t& scheme) const
        -> node_set_t requires(idx_t::supports_locate)
    {
        band_k = scheme.k;
        if (scheme.k <= edit_distance_matrix<uint64_t>::k_limit)
             return search<uint64_t>(P, scheme);
        else return search<__uint128_t>(P, scheme);
    }

    /**
     * @brief computes the maximal sub-intervals of the union of the given (sorted) SA-intervals such that all
     *        suffixes in a sub-interval are reported as the same minimum-error (longest on ties) occurrence
     * @param ctxts_sorted the LOCATE contexts sorted by their forward SA-interval starting position
     * @return the effective intervals, each storing its [beg, end] range and its (len, err)
     */
    static auto effective_intervals(const std::vector<sorted_ctx_t>& ctxts_sorted)
        -> std::vector<interval_t>
    {
        std::vector<interval_t> ivs;
        ivs.reserve(ctxts_sorted.size());

        // SA-intervals are laminar (nested or disjoint), so a left-to-right sweep keeps the currently open
        // intervals on a stack (LIFO by their end). At each event position the segment since the last event is
        // closed off into one effective interval -- the minimum-error (longest on ties) of the open contexts.
        // end_excl = one past the last SA-position; idx = the winner's CIGAR index (CIGAR mode)
        struct open_t {pos_t end_excl; pos_t len; pos_t err; pos_t idx;}; 
        std::vector<open_t> open;
        open.reserve(16);
        pos_t seg_beg = 0;

        for (uint64_t i = 0; i < ctxts_sorted.size() || !open.empty();) {
            pos_t next_start = i < ctxts_sorted.size()
                ? std::get<0>(ctx_of(ctxts_sorted[i]).ctx.forward_sa_interval()) : pos_t(-1);
            pos_t next_end = open.empty() ? pos_t(-1) : open.back().end_excl;
            pos_t pos = std::min(next_start, next_end);

            if (!open.empty() && seg_beg < pos) {
                const open_t* best = &open.front();
                for (const open_t& o : open)
                    if (o.err < best->err || (o.err == best->err && o.len > best->len)) best = &o;
                ivs.emplace_back(interval_t{seg_beg, pos - 1, best->err, best->len, best->idx});
            }

            while (!open.empty() && open.back().end_excl == pos) open.pop_back();
            while (i < ctxts_sorted.size() && std::get<0>(ctx_of(ctxts_sorted[i]).ctx.forward_sa_interval()) == pos) {
                const node_t& node = ctx_of(ctxts_sorted[i]);
                pos_t end = std::get<1>(node.ctx.forward_sa_interval());
                open.emplace_back(open_t{end + 1, node.match_len(), node.ctx.errors(), cig_of(ctxts_sorted[i])});
                i++;
            }

            seg_beg = pos;
        }

        return ivs;
    }

protected:
    // whether the extend context's next symbol is the index's separator. The separator is the largest symbol of the
    // alphabet, so the O(sigma) extend loop reaches it last; a search must never extend into it (cross a sequence
    // boundary), so the loop stops there. Always false for an index that defines no separator (e.g. a non-FASTA index
    // or a benchmark adapter), where the separator query is not even instantiated.
    template <typename ext_ctx_t>
    bool next_is_separator(const ext_ctx_t& ext_ctx) const
    {
        if constexpr (requires (const idx_t& ix) { ix.separator_sym(); ix.has_sequences(); }) {
            return index.has_sequences() && ext_ctx.can_extend() && ext_ctx.next_symbol() == index.separator_sym();
        } else {
            return false;
        }
    }

    /**
     * @brief executes a search scheme (band-sharing), choosing the matrix word type
     * @tparam matrix_word_t word type of the edit-distance matrix
     * @param P the pattern to search
     * @param scheme the search scheme to use
     * @return all relaxed search contexts of P in T with at most k errors
     */
    template <typename matrix_word_t>
    auto search(const inp_t& P, const search_scheme_t& scheme) const
        -> node_set_t requires(idx_t::supports_locate)
    {
        pos_t m = P.size();
        pos_t p = scheme.p;
        if (m == 0) return {index.template empty_context<LOCATE>()};
        if (p >= m) throw std::runtime_error("p >= m");

        active_pool = this; match_pool.clear(); match_free.clear();

        node_set_t result;
        std::vector<edit_dist_search> searches;
        std::vector<directional_substring<pos_t, inp_t>> parts;
        searches.reserve(scheme.S.size());
        parts.reserve(p);

        // cut P into the p contiguous parts addressed by the scheme
        for (pos_t p_idx = 0; p_idx < p; p_idx++) {
            pos_t part = scheme.S[0][p_idx].part;
            pos_t beg = (part * m) / p;
            pos_t end = (part == p - 1 ? m : (((part + 1) * m) / p)) - 1;
            parts.emplace_back(directional_substring<pos_t, inp_t>(P, beg, end));
        }

        for (pos_t s_idx = 0; s_idx < scheme.S.size(); s_idx++) {
            searches.emplace_back(edit_dist_search(scheme, s_idx));
        }

        std::vector<std::vector<node_t>> stacks;
        stacks.resize(p);

        std::vector<edit_distance_matrix<matrix_word_t>> matrices_fwd, matrices_bwd;
        matrices_fwd.resize(p);
        matrices_bwd.resize(p);

        for (const edit_dist_search& search : searches) {
            data_struct_ref_t<matrix_word_t> ds_ref(search, parts, stacks, matrices_fwd, matrices_bwd, result);

            for (pos_t p_idx = 0; p_idx < p; p_idx++) {
                parts[search.part(p_idx)].set_direction(search.part_dir(p_idx));
            }

            node_t node(index.template empty_context<LOCATE>());
            pos_t p_idx = 0;

            if (search.upper_bound(0) == 0) {
                bool match = true;

                for (; match && p_idx < parts.size(); p_idx++) {
                    if (search.upper_bound(p_idx) > 0) break;
                    const auto& part = parts[search.part(p_idx)];

                    for (pos_t i = 0; match && i < part.size(); i++) {
                        auto [ctx_nxt, match_nxt] = node.ctx.extend(part[i], part.direction());
                        if (match_nxt) {
                            node.ctx = ctx_nxt;
                            node.set_match_len(node.match_len() + 1);
                            // an exact-part symbol is a known match at P[part.pos(i)]
                            run_match(node, part.pos(i), part.direction());
                        }
                        match = match_nxt;
                    }
                }

                if (!match) continue;
            }

            rec_search<matrix_word_t>(node, p_idx, ds_ref, {}, {}, {}, {});
        }

        return result;
    }

    /**
     * @brief recursively processes the p_idx-th part of the search by extending the current context in
     *        the part's direction and branching/bounding on the edit-distance matrix, carrying the previous
     *        part's band forward as descendants/init_dists when the direction does not switch
     * @tparam matrix_word_t word type of the edit-distance matrix
     * @param node the current search context
     * @param p_idx index of the part to process
     * @param ds_ref references to the data structures shared during the search
     * @param desc_prev_dir descendant contexts from the previous part in the same direction
     * @param desc_diff_dir descendant contexts from the previous part in the opposite direction
     * @param dists_prev_dir edit distances corresponding to desc_prev_dir
     * @param dists_diff_dir edit distances corresponding to desc_diff_dir
     */
    template <typename matrix_word_t>
    void rec_search(
        node_t& node, pos_t p_idx,
        data_struct_ref_t<matrix_word_t>& ds_ref,
        const std::vector<node_t>& desc_prev_dir,
        const std::vector<node_t>& desc_diff_dir,
        const std::vector<uint16_t>& dists_prev_dir,
        const std::vector<uint16_t>& dists_diff_dir
    ) const requires(idx_t::supports_locate)
    {
        auto& [search, parts, stacks, matrices_fwd, matrices_bwd, result] = ds_ref;

        const auto& part = parts[search.part(p_idx)];
        pos_t k_max = search.upper_bound(p_idx);
        direction_t dir = search.part_dir(p_idx);
        bool dir_switch = search.does_part_switch_dir(p_idx);
        auto& stack = stacks[p_idx];
        auto& matrix = dir == RIGHT ? matrices_fwd[search.part(p_idx)] : matrices_bwd[search.part(p_idx)];

        // on a direction switch the roles of the two carried bands swap: the band accumulated in the now-current
        // direction (the "diff" band from before the switch) becomes the same-direction band that seeds this matrix,
        // and the previous part's band becomes the reverse one
        const auto& dists_dir = dir_switch ? dists_diff_dir : dists_prev_dir;
        const auto& desc_dir = dir_switch ? desc_diff_dir : desc_prev_dir;
        const auto& dists_rev_dir = dir_switch ? dists_prev_dir : dists_diff_dir;
        const auto& desc_rev_dir = dir_switch ? desc_prev_dir : desc_diff_dir;

        // seed the matrix' first column: a single cell on a restart, otherwise the previous frontier's distances
        // shifted by the error increment of the current context
        std::vector<uint16_t> dists_new;

        if (dists_dir.empty()) {
            dists_new.emplace_back(node.ctx.errors());
        } else {
            uint16_t prev_dist = dir_switch ? *min_element(dists_dir.begin(), dists_dir.end()) : dists_dir[0];
            uint16_t add = node.ctx.errors() - prev_dist;
            dists_new.reserve(dists_dir.size());

            for (pos_t i = 0; i < dists_dir.size(); i++) {
                dists_new.emplace_back(dists_dir[i] + add);
            }
        }

        if (!matrix.is_initialized()) matrix.set_input(part, index.sigma, index.forward_index().map_int());
        matrix.init(k_max, dists_new);
        cluster_t cluster(matrix.last_col_size(), k_max, node.match_len());

        if (matrix.is_in_final_column(0)) {
            auto node_copy = node;
            node_copy.set_match_len(0);
            cluster.set_cell(0, flushed(node_copy), matrix(0, part.size()));
        }

        pos_t old_depth = 0;

        // replay the carried-over descendant rows of the band before fanning out into new text symbols
        if (!desc_dir.empty()) {
            for (pos_t i = 0; i < desc_dir.size() && desc_dir[i].match_len() < matrix.num_rows(); i++) {
                if (branch_and_bound<matrix_word_t>(matrix, cluster, desc_dir[i], p_idx, ds_ref,
                    dists_rev_dir, desc_rev_dir, {desc_dir.begin() + i + 1, desc_dir.end()})
                ) return;
            }

            if (desc_dir.back().match_len() == matrix.num_rows() - 1) return;
            if (!dir_switch) node = desc_dir.back();
            old_depth = desc_dir.back().match_len();
        }

        auto ext_ctx = node.ctx.extend_phase();
        ext_ctx.prepare_extend_all(dir);
        bool extended = false;

        while (ext_ctx.can_extend()) {
            if (next_is_separator(ext_ctx)) break; // the separator is the largest symbol, hence last: nothing follows
            stack.emplace_back(ext_ctx.extend_next());
            stack.back().set_match_len(old_depth + 1);

            // inherit the parent's chain + pending run; the symbol is recorded when this child is popped (see the DFS loop)
            inherit_match(stack.back(), node);
            extended = true;
        }

        if (!extended && matrix.is_in_final_column(0)) {
            flush_cluster<matrix_word_t>(cluster, p_idx, ds_ref, dists_rev_dir, desc_rev_dir, {});
            return;
        }

        // DFS over the text continuations: compute each row, prune via branch_and_bound, extend the remaining contexts
        while (!stack.empty()) {
            auto node_cur = stack.back(); stack.pop_back();
            // record this new extension's matched symbol (run or literal) before it is reported/carried/extended
            run_update(node_cur, part, dir);
            if (branch_and_bound<matrix_word_t>(matrix, cluster, node_cur, p_idx, ds_ref, dists_rev_dir, desc_rev_dir, {})) continue;
            auto ext_ctx_cur = node_cur.ctx.extend_phase();
            ext_ctx_cur.prepare_extend_all(dir);

            while (ext_ctx_cur.can_extend()) {
                if (next_is_separator(ext_ctx_cur)) break; // do not extend across a sequence boundary
                stack.emplace_back(ext_ctx_cur.extend_next());
                stack.back().set_match_len(node_cur.match_len() + 1);

                // children inherit the parent's matched-string chain and pending run; their own symbol is recorded
                // when they are popped (above), so no pool cell is created for the fanned-out (mostly pruned) children
                inherit_match(stack.back(), node_cur);
            }
        }
    }

    /**
     * @brief flushes a completed cluster (a part's final matrix column): for the last part it reports the
     *        occurrences, for an inner edge/non-edge part it recurses into the next part with the appropriate
     *        cluster center(s)
     * @tparam matrix_word_t word type of the edit-distance matrix
     * @param cluster the cluster collecting the cells of the matrix' final column
     * @param p_idx index of the current part
     * @param ds_ref references to the data structures shared during the search
     * @param dists_rev_dir edit distances of the descendants in the reverse direction
     * @param desc_rev_dir descendant contexts in the reverse direction
     * @param desc_dir_rem remaining descendant contexts in the current direction
     */
    template <typename matrix_word_t>
    void flush_cluster(
        cluster_t& cluster, pos_t p_idx,
        data_struct_ref_t<matrix_word_t>& ds_ref,
        const std::vector<uint16_t>& dists_rev_dir,
        const std::vector<node_t>& desc_rev_dir,
        const std::vector<node_t>& desc_dir_rem
    ) const requires(idx_t::supports_locate)
    {
        auto& [search, parts, stacks, matrices_fwd, matrices_bwd, result] = ds_ref;
        const auto& k_min = search.lower_bound(p_idx);

        // last part: report the winning cells, keeping the minimum-error (LEFT-preferred on ties) one
        auto report_at_end = [&]() {
            auto centers = cluster.centers();

            for (auto& node_cur : centers) {
                if (node_cur.ctx.errors() >= k_min) {
                    auto it = result.find(node_cur);

                    if (it != result.end()) {
                        const auto& node_old = *it;
                        uint16_t cur_err = node_cur.ctx.errors();
                        uint16_t old_err = node_old.ctx.errors();

                        if (cur_err  < old_err ||
                           (cur_err == old_err && node_old.ctx.s.d == RIGHT && node_cur.ctx.s.d == LEFT)
                        ) {
                            it = result.erase(it);
                            result.emplace_hint(it, node_cur);
                        }
                    } else {
                        result.emplace(node_cur);
                    }
                }
            }
        };

        // at a direction switch or edge, continue from the local-minimum anchor cells; otherwise (same
        // direction) carry the cluster center + its band forward
        bool next_switches = (p_idx + 1 < parts.size()) && search.does_part_switch_dir(p_idx + 1);

        if (search.is_edge(p_idx) || next_switches) {
            if (p_idx + 1 == parts.size()) {
                report_at_end();
            } else {
                auto centers = cluster.centers();

                for (auto& node : centers)
                    if (node.ctx.errors() >= k_min)
                        rec_search<matrix_word_t>(node, p_idx + 1, ds_ref, {}, desc_rev_dir, {}, dists_rev_dir);
            }

            return;
        }

        std::vector<node_t> desc_dir;
        std::vector<uint16_t> dists_dir;
        auto [node_cntr, res] = cluster.cluster_center(k_min, desc_dir, dists_dir);

        // no strict local minimum >= k_min: this subtree is covered by another search, so it is not continued
        if (!res) return;

        desc_dir.insert(desc_dir.end(), desc_dir_rem.begin(), desc_dir_rem.end());

        for (pos_t i = 0; i < desc_dir.size(); i++) {
            desc_dir[i].set_match_len(i + 1);
        }

        pos_t k_max = search.upper_bound(p_idx + 1);
        while (dists_dir.back() > k_max) dists_dir.pop_back();

        // same direction continues: carry the center + descendants forward into the next part
        rec_search<matrix_word_t>(node_cntr, p_idx + 1, ds_ref, desc_dir, desc_rev_dir, dists_dir, dists_rev_dir);
    }

    /**
     * @brief computes the matrix row for node and, when the final column is reached, records the cell in the
     *        cluster and (if the cell is terminal) flushes the cluster; bounds the search if no valid distance remains
     * @tparam matrix_word_t word type of the edit-distance matrix
     * @param matrix the edit-distance matrix of the current part
     * @param cluster the cluster collecting the cells of the matrix' final column
     * @param node the current search context
     * @param p_idx index of the current part
     * @param ds_ref references to the data structures shared during the search
     * @param dists_rev_dir edit distances of the descendants in the reverse direction
     * @param desc_rev_dir descendant contexts in the reverse direction
     * @param desc_dir_rem remaining descendant contexts in the current direction
     * @return true <=> the subtree rooted at node has been fully processed (it need not be extended further)
     */
    template <typename matrix_word_t>
    bool branch_and_bound(
        edit_distance_matrix<matrix_word_t>& matrix, cluster_t& cluster,
        const node_t& node, pos_t p_idx,
        data_struct_ref_t<matrix_word_t>& ds_ref,
        const std::vector<uint16_t>& dists_rev_dir,
        const std::vector<node_t>& desc_rev_dir,
        const std::vector<node_t>& desc_dir_rem
    ) const requires(idx_t::supports_locate)
    {
        auto& [search, parts, stacks, matrices_fwd, matrices_bwd, result] = ds_ref;
        pos_t row = node.match_len();

        if (row >= matrix.num_rows()) return true;
        bool valid_dist = matrix.compute_row(row, node.ctx.last_added_symbol());

        if (matrix.is_in_final_column(row)) {
            pos_t cluster_idx = row - (matrix.num_rows() - cluster.size());
            cluster.set_cell(cluster_idx, flushed(node), matrix(row, matrix.num_cols() - 1));

            // a cell is terminal when the band holds no further distance, only vertical gaps remain, or the index
            // cannot be extended any further
            bool terminal = !valid_dist || matrix.only_vertical_gaps_left(row);

            if (!terminal) {
                // a cell whose only remaining extension is the separator (a sequence boundary the search must not
                // cross) cannot be extended further, so it is terminal and is reported here
                node_t node_ext = node;
                auto ext_ctx = node_ext.ctx.extend_phase();
                ext_ctx.prepare_extend_all(search.part_dir(p_idx));
                terminal = !ext_ctx.can_extend() || next_is_separator(ext_ctx);
            }

            if (terminal) {
                flush_cluster<matrix_word_t>(cluster, p_idx, ds_ref, dists_rev_dir, desc_rev_dir, desc_dir_rem);
                return true;
            }
        }

        return !valid_dist;
    }
};

// deduction guide so that apm_edit{index} deduces the index type
template <typename idx_t>
apm_edit(const idx_t&) -> apm_edit<idx_t>;
