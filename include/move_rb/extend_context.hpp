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

#ifndef MOVE_R_MOVE_RB_EXTEND_CONTEXT_HPP
#define MOVE_R_MOVE_RB_EXTEND_CONTEXT_HPP

#include <move_rb/move_rb.hpp>

/**
 * @brief the state of extending a search context with every symbol in one direction: prepare_extend_all
 *        computes the symbols c that occur in the range of input intervals of the pattern, together with the
 *        first and the last input interval with c (O(sigma) time, or O(e_ - b_) for a narrow range); then each
 *        extend_next extends with the next symbol with O(1) move queries
 * @tparam query_support whether the search context supports COUNT or LOCATE queries
 */
template <move_r_support support>
template <query_support_t query_support>
struct move_rb<support>::extend_context_t {
    friend class move_rb;
    friend struct search_context_t<query_support>;

    // ############################# TYPES #############################

protected:
    // how the prepared extensions are enumerated:
    // - _prev_next: from the arrays prev and next (see build_prev_next) in prev_next_buf; cur is the symbol
    // - _symbol_list: from syms; cur is an index in syms
    // - _table_group: from a range of entries of a level of the context tables (see context_tables::extension_range); cur
    //   is the entry
    enum form_t : uint8_t { _prev_next, _symbol_list, _table_group };

    // ############################# DATA MEMBERS #############################

    const move_rb<support>* idx = nullptr; // the index this context queries
    search_context_t<query_support>* ctx = nullptr; // the search context that is extended

    // keys of the level of the extensions in the context tables (see context_tables); the key of an extension is
    // (entry of the pattern) * sigma + symbol
    const rank_select_set* tbl_keys = nullptr;
    uint64_t tbl_base = 0; // key of the extension with the symbol 0 (entry of the pattern * sigma)
    uint64_t b_R_nxt; // left limit of the interval of the other direction of the next extension
    uint64_t cur; // the next symbol (_prev_next), its index in syms (_symbol_list) or its entry (_table_group)

    uint32_t tbl_sym = 0; // symbol of the next extension (_table_group)
    uint32_t tbl_len = 0; // length of the extended patterns (their level in the context tables)
    uint32_t tbl_end = 0; // last entry of the extensions + 1 (_table_group)
    form_t form = _prev_next; // how the prepared extensions are enumerated
    direction_t dir; // extend direction (set by prepare_extend_all)

    /* the symbols of the extensions (_symbol_list): [0..k-1], one entry per distinct symbol c of L'[b_..e_], sorted
       by c (see symbol_runs_t); prepare_extend_all fills it temporarily with one entry per input interval */
    std::vector<symbol_runs_t> syms;
    // prev[0..sigma-1] followed by next[0..sigma-1] (_prev_next, see build_prev_next)
    std::vector<int64_t> prev_next_buf;

    // ############################# INTERNAL METHODS #############################

    // returns the pool of unused buffers of the current thread; syms and prev_next_buf are taken from the
    // pool and returned to it by the destructor
    template <typename val_t>
    static std::vector<std::vector<val_t>>& buffer_pool()
    {
        static thread_local std::vector<std::vector<val_t>> pool;
        return pool;
    }

    // takes a buffer from the pool of the current thread, unless buffer already has capacity
    template <typename val_t>
    static void take_buffer(std::vector<val_t>& buffer)
    {
        std::vector<std::vector<val_t>>& pool = buffer_pool<val_t>();
        if (buffer.capacity() != 0 || pool.empty()) return;

        buffer = std::move(pool.back());
        pool.pop_back();
    }

    // returns a buffer to the pool of the current thread
    template <typename val_t>
    static void return_buffer(std::vector<val_t>& buffer)
    {
        if (buffer.capacity() == 0) return;
        buffer_pool<val_t>().emplace_back(std::move(buffer));
    }

    // return the arrays prev and next (see prev_next_buf)
    inline int64_t* prev_arr() { return prev_next_buf.data(); }
    inline int64_t* next_arr() { return prev_next_buf.data() + idx->sigma; }
    inline const int64_t* prev_arr() const { return prev_next_buf.data(); }
    inline const int64_t* next_arr() const { return prev_next_buf.data() + idx->sigma; }

    // resizes prev_next_buf to store prev and next
    inline void resize_prev_next_buf() { take_buffer(prev_next_buf); prev_next_buf.resize(2 * idx->sigma); }

    // returns the symbol of the next extension in the internal alphabet
    inline uint64_t next_internal_symbol() const
    { return form == _table_group ? tbl_sym : (form == _symbol_list ? syms[cur].sym : cur); }

    // sets tbl_sym to the symbol of the extension at entry cur
    inline void read_table_symbol() { if (cur < tbl_end) tbl_sym = tbl_keys->select_1(cur + 1) - tbl_base; }

    // returns the first input interval with the symbol of the next extension in the range of the pattern
    inline uint64_t next_first_run() const
    { return form == _symbol_list ? syms[cur].first : uint64_t(next_arr()[cur]); }

    // returns the last input interval with the symbol of the next extension in the range of the pattern
    inline uint64_t next_last_run() const
    { return form == _symbol_list ? syms[cur].last : uint64_t(prev_arr()[cur]); }

    /**
     * @brief advances cur to the next symbol the pattern can be extended with (or past the last one)
     * @tparam dir extend direction
     */
    template <direction_t dir>
    void advance_symbol();

    // ############################# CONSTRUCTORS #############################

public:
    extend_context_t() {}

    // constructs the extend context of the search context ctx of the index idx
    extend_context_t(const move_rb<support>& idx, search_context_t<query_support>& ctx)
        : idx(&idx), ctx(&ctx)
    { }

    // returns the buffers to the pool (see buffer_pool)
    ~extend_context_t() { return_buffer(syms); return_buffer(prev_next_buf); }

    extend_context_t(extend_context_t&&) = default;
    extend_context_t& operator=(extend_context_t&&) = default;
    extend_context_t(const extend_context_t&) = default;
    extend_context_t& operator=(const extend_context_t&) = default;

    // ############################# QUERY METHODS #############################

    // returns the symbol the next extend_next() extends with (only valid while can_extend())
    uint64_t next_symbol() const { return idx->unmap_symbol(next_internal_symbol()); }

    // returns whether there is a symbol left to extend with
    bool can_extend() const
    {
        if (form == _table_group) return cur < tbl_end;
        return form == _symbol_list ? cur < syms.size() : cur < idx->sigma;
    }

    /**
     * @brief prepares the extensions of the search context in direction dir (see extend_next)
     * @tparam scope which symbols are enumerated (see extension_scope_t)
     * @param dir extend direction
     * @return this extend context
     */
    template <extension_scope_t scope = TERM>
    extend_context_t& prepare_extend_all(direction_t dir);

    /**
     * @brief prepares the extensions of the search context in direction dir (see extend_next)
     * @tparam dir extend direction
     * @tparam scope which symbols are enumerated (see extension_scope_t)
     */
    template <direction_t dir, extension_scope_t scope = TERM>
    void prepare_extend_all();

    // returns the search context extended with the next symbol (in the direction set by prepare_extend_all)
    search_context_t<query_support> extend_next();

    /**
     * @brief returns the search context extended with the next symbol
     * @tparam dir extend direction
     */
    template <direction_t dir>
    search_context_t<query_support> extend_next();
};

#include <move_rb/extend_context.tpp>

#endif // MOVE_R_MOVE_RB_EXTEND_CONTEXT_HPP
