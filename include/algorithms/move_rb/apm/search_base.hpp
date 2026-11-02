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

#ifndef MOVE_R_ALGORITHMS_MOVE_RB_APM_SEARCH_BASE_HPP
#define MOVE_R_ALGORITHMS_MOVE_RB_APM_SEARCH_BASE_HPP

#include <cstdint>
#include <optional>
#include <vector>

#include <gtl/phmap.hpp>
#include <misc/apm.hpp>
#include <misc/cigar.hpp>

/**
 * @brief returns the direction in which each step of a search matches its part: RIGHT if its part lies right of the
 *        part of the previous step, else LEFT; the first step has the direction of the second one (RIGHT if p = 1)
 * @param S the search
 * @return [0..p-1] the direction of step i
 */
inline std::vector<direction_t> apm_step_directions(const search_t& S)
{
    std::vector<direction_t> dirs(S.size());
    dirs[0] = S.size() == 1 || S[1].part > S[0].part ? RIGHT : LEFT;
    for (uint64_t i = 1; i < S.size(); i++) dirs[i] = S[i].part > S[i - 1].part ? RIGHT : LEFT;
    return dirs;
}

/**
 * @brief a set of search contexts with distinct keys, in the order of their insertion: the contexts are stored in a
 *        vector and a hash map maps the key of each one to its index, so that an insertion moves no context
 * @tparam ctx_t type of a context
 * @tparam key_t type of the key of a context (the values that identify it)
 * @tparam key_hash_t hash function of the keys
 */
template <typename ctx_t, typename key_t, typename key_hash_t = gtl::Hash<key_t>>
class apm_context_set {
protected:
    std::vector<ctx_t> _contexts; // [0..size()-1] the contexts in the order of their insertion
    gtl::flat_hash_map<key_t, uint32_t, key_hash_t> _index_of; // maps the key of each context to its index in _contexts

public:
    // returns the context with a key, or nullptr if there is none
    ctx_t* find(const key_t& key)
    {
        auto it = _index_of.find(key);
        return it == _index_of.end() ? nullptr : &_contexts[it->second];
    }

    // inserts a context with a key that no context of the set has
    void insert(const key_t& key, const ctx_t& ctx)
    {
        _index_of.emplace(key, uint32_t(_contexts.size()));
        _contexts.emplace_back(ctx);
    }

    // returns the number of contexts
    uint64_t size() const { return _contexts.size(); }

    // removes all contexts
    void clear()
    {
        _contexts.clear();
        _index_of.clear();
    }

    // returns the contexts in the order of their insertion
    std::vector<ctx_t>& contexts() { return _contexts; }

    auto begin() { return _contexts.begin(); }
    auto end() { return _contexts.end(); }
    auto begin() const { return _contexts.begin(); }
    auto end() const { return _contexts.end(); }
};

/**
 * @brief the children of a node of a depth-first search in the order of their creation: the first one is kept, and a
 *        second one pushes both onto the stack in this order, as every further one is pushed; a node with one child
 *        continues with it without the stack, and the nodes are visited in the same order
 * @tparam child_t type of a child
 */
template <typename child_t>
class apm_dfs_children {
protected:
    std::optional<child_t> _first; // the first child, until a second one is added
    uint64_t _num = 0; // the number of children

public:
    /**
     * @brief adds a child
     * @tparam push_fnc_t type of the function push
     * @param child the child
     * @param push called as push(child) for each child that is pushed onto the stack
     */
    template <typename push_fnc_t>
    void add(child_t&& child, push_fnc_t&& push)
    {
        if (_num++ == 0) { _first.emplace(std::move(child)); return; }
        if (_num == 2) push(std::move(*_first));
        push(std::move(child));
    }

    // returns whether there is exactly one child
    bool has_one() const { return _num == 1; }

    // returns the only child (if has_one())
    child_t& only() { return *_first; }
};

/**
 * @brief the members of apm_hamming and apm_edit: the index, the query parameters, the locate query without errors
 *        and the test for the separator of the sequences
 * @tparam idx_t the bidirectional index type (move_rb, or a benchmark adapter with the same interface)
 */
template <typename idx_t>
class apm_search_base {
protected:
    using inp_t = typename idx_t::inp_t;

    template <query_support_t query_support>
    using search_context_t = typename idx_t::template search_context_t<query_support>;

    const idx_t& _index; // the index
    apm_params_t _params; // the query parameters

    apm_search_base(const idx_t& index, apm_params_t params) : _index(index), _params(params) {}

    /**
     * @brief returns whether a symbol is the separator of the sequences of the index (always false for an index
     *        without sequences); a search does not extend with it
     * @param sym a symbol
     */
    bool is_separator(uint64_t sym) const
    {
        if constexpr (requires (const idx_t& ix) { ix.seq_data(); }) {
            return _index.seq_data().has_sequences() && sym == _index.seq_data().separator_sym();
        } else {
            return false;
        }
    }

    /**
     * @brief returns whether the next symbol of an extension context is the separator of the sequences of the index
     *        (see is_separator)
     * @tparam ext_ctx_t type of the extension context
     * @param ext_ctx an extension context
     */
    template <typename ext_ctx_t>
    bool next_is_separator(const ext_ctx_t& ext_ctx) const
    {
        if constexpr (requires (const idx_t& ix) { ix.seq_data(); }) {
            return ext_ctx.can_extend() && is_separator(ext_ctx.next_symbol());
        } else {
            return false;
        }
    }

    /**
     * @brief locates the occurrences of P without errors
     * @tparam mode whether each occurrence is reported with a CIGAR index
     * @tparam report_fnc_t type of the report function
     * @param P the pattern
     * @param report called as report({pos, m, 0}) or, in CIGAR mode, report({pos, m, 0, 0}) for each occurrence
     * @return in CIGAR mode, the CIGAR of the occurrences (m matches); empty in NO_CIGAR mode
     */
    template <cigar_mode_t mode, typename report_fnc_t>
    std::vector<cigar_t<uint8_t>> locate_exact(const inp_t& P, report_fnc_t& report) const
    {
        const uint64_t m = P.size();
        std::vector<cigar_t<uint8_t>> cigars;

        if constexpr (mode == CIGAR) {
            cigars.push_back(m > 0 ? cigar_t<uint8_t>{MATCH<uint8_t>(m)} : cigar_t<uint8_t>{});
            _index.forward_index().locate(P, [&](uint64_t, uint64_t occ){report({occ, m, 0, uint32_t(0)});});
        } else {
            _index.forward_index().locate(P, [&](uint64_t, uint64_t occ){report({occ, m, 0});});
        }

        return cigars;
    }
};

#endif // MOVE_R_ALGORITHMS_MOVE_RB_APM_SEARCH_BASE_HPP
