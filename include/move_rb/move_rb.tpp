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

#ifndef MOVE_R_MOVE_RB_MOVE_RB_TPP
#define MOVE_R_MOVE_RB_MOVE_RB_TPP

#include <move_rb/move_rb.hpp>

// ############################# MISC METHODS #############################

template <move_r_support support>
uint64_t move_rb<support>::size_in_bytes() const
{
    uint64_t size =
        idx_fwd.size_in_bytes() +
        idx_bwd.size_in_bytes() +
        _S_MLF_p_fwd.size_in_bytes() +
        _S_MLF_p_bwd.size_in_bytes() +
        _S_MPhi_p_fwd.size_in_bytes() +
        _S_MPhi_inv_p_fwd.size_in_bytes() +
        opt_size(_SA_pos) +
        _seq_data.size_in_bytes() +
        _context_tables.size_in_bytes();

    return size;
}

template <move_r_support support>
void move_rb<support>::log_construction_summary(uint64_t time_construction, uint64_t peak_mem_usage,
    std::ostream* mf) const
{
    // the overall peak, without subtracting the memory allocated before the construction; no throughput is logged,
    // if the pfp output was read from a cache (see move_r::rlbwt_from_cache)
    (void) peak_mem_usage;
    ::log_construction_summary(time_construction, n,
        is_malloc_count_linked ? overall_peak() : peak_not_measured, _pfp_cached_texts == 0);

    if (mf != nullptr)
        *mf << " throughput_in_mb_per_s="
            << throughput_in_mb_per_s(n, time_construction)
            << " bigbwt_cached_texts=" << int(_pfp_cached_texts);
}

template <move_r_support support>
template <direction_t dir>
void move_rb<support>::log_direction_sizes() const
{
    std::cout << (dir == LEFT ? "forward index:" : "backward index:") << std::endl;
    listing_indent++;
    index<dir>().log_data_structure_sizes();

    log_component("S_MLF_p", S_MLF_p<dir>().size_in_bytes());

    if constexpr (dir == LEFT && has_mphi) {
        log_component("S_MPhi_p", _S_MPhi_p_fwd.size_in_bytes());
        log_component("S_MPhi^{-1}_p", _S_MPhi_inv_p_fwd.size_in_bytes());
    }

    if constexpr (dir == RIGHT && has_SA_pos) log_component("SA_pos", opt_size(_SA_pos));
    listing_indent--;
}

template <move_r_support support>
void move_rb<support>::log_built_context_tables() const
{
    if (_context_tables.empty()) {
        std::cout << "no context tables built (the budget is too small for one level)" << std::endl;
    } else if (_context_tables.single_length() != 0) {
        std::cout << "context tables of lengths 1 to " << _context_tables.max_length() << " and "
                  << _context_tables.single_length() << " alone, using " << _context_tables.bits_per_pattern()
                  << " bits/context" << std::endl;
    } else {
        std::cout << "context tables (up to length " << _context_tables.max_length() << "), using "
                  << _context_tables.bits_per_pattern() << " bits/context" << std::endl;
    }
}

template <move_r_support support>
void move_rb<support>::log_data_structure_sizes() const
{
    const bool owner = begin_listing(size_in_bytes());
    log_direction_sizes<LEFT>();
    log_direction_sizes<RIGHT>();
    log_context_tables_size();
    if (_seq_data.has_sequences()) log_component("FASTA sequences", _seq_data.size_in_bytes());
    if (owner) end_listing();
}

template <move_r_support support>
void move_rb<support>::add_parameters(index_facts& facts) const
{
    idx_fwd.add_parameters(facts);
    idx_bwd.add_parameters(facts, "^R");
}

template <move_r_support support>
void move_rb<support>::log_data_structures() const
{
    if (n > max_logged_input_length) {
        std::cout << "cannot log the contents for n > " << max_logged_input_length << std::endl;
        return;
    }

    std::cout << "Forward Index:" << std::endl;
    idx_fwd.log_data_structures();

    std::cout << std::endl << "Backward Index:" << std::endl;
    idx_bwd.log_data_structures();

    if constexpr (has_SA_pos) {
        std::cout << std::endl;
        aligned_table log_pos;
        log_pos.add_row("SA_pos:", _SA_pos.size(), [&](uint64_t i) { return _SA_pos[i]; });
        log_pos.print();
    }
}

template <move_r_support support>
void move_rb<support>::log_data_structure_sizes(std::ostream& out) const
{
    out << " sigma=" << sigma;
    out << " r=" << idx_fwd.num_bwt_runs();
    out << " r_=" << idx_fwd.M_LF().num_intervals();

    out << " r_rev=" << idx_bwd.num_bwt_runs();
    out << " r_rev_=" << idx_bwd.M_LF().num_intervals();

    if constexpr (supports_multiple_locate) {
        if constexpr (has_locate_move) {
            out << " r__=" << idx_fwd.M_Phi_inv().num_intervals();
            out << " r___=" << idx_fwd.M_Phi().num_intervals();
        } else if constexpr (has_SA_rlz) {
            out << " z=" << idx_fwd.num_SA_rlz_phrases();
            out << " z_l=" << idx_fwd.num_SA_rlz_literal_phrases();
            out << " z_c=" << idx_fwd.num_SA_rlz_copy_phrases();
        }
    }

    out << " size_index=" << size_in_bytes();
    idx_fwd.log_data_structure_sizes(out, "_fwd");
    out << " size_s_mlf_p_fwd=" << _S_MLF_p_fwd.size_in_bytes();

    if constexpr (has_mphi) {
        out << " size_s_mphi_p=" << _S_MPhi_p_fwd.size_in_bytes();
        out << " size_s_mphi_inv_p=" << _S_MPhi_inv_p_fwd.size_in_bytes();
    }

    idx_bwd.log_data_structure_sizes(out, "_bwd");
    out << " size_s_mlf_p_bwd=" << _S_MLF_p_bwd.size_in_bytes();
    if constexpr (has_SA_pos) out << " size_sa_pos=" << opt_size(_SA_pos);

    if (has_context_tables()) {
        out << " context_tables_max_length=" << context_tables_max_length();
        out << " size_context_tables=" << _context_tables.size_in_bytes();
    }

}

// ############################# SERIALIZATION METHODS #############################

template <move_r_support support>
void move_rb<support>::serialize(std::ostream& out) const
{
    move_r_support _support = support;
    // the type of locate support (see read_move_rb_support)
    out.write((char*) &_support, sizeof(move_r_support));

    out.write((char*) &n, sizeof(uint64_t));
    out.write((char*) &sigma, sizeof(uint64_t));

    idx_fwd.serialize(out);
    idx_bwd.serialize(out);

    _S_MLF_p_fwd.serialize(out);
    _S_MLF_p_bwd.serialize(out);

    _S_MPhi_p_fwd.serialize(out);
    _S_MPhi_inv_p_fwd.serialize(out);

    opt_serialize(_SA_pos, out);

    _seq_data.serialize(out);
    _context_tables.serialize(out);
}

template <move_r_support support>
void move_rb<support>::load(std::istream& in)
{
    move_r_support _support;
    in.read((char*) &_support, sizeof(move_r_support));

    if (_support != support) {
        throw std::runtime_error("move_rb: the index was built with another locate support (" +
            std::to_string(int(_support)) + ") than the index type it is loaded into (" +
            std::to_string(int(support)) + ")");
    }

    in.read((char*) &n, sizeof(uint64_t));
    in.read((char*) &sigma, sizeof(uint64_t));

    idx_fwd.load(in);
    idx_bwd.load(in);

    _S_MLF_p_fwd.load(in);
    _S_MLF_p_bwd.load(in);

    _S_MPhi_p_fwd.load(in);
    _S_MPhi_inv_p_fwd.load(in);

    opt_load(_SA_pos, in);

    _seq_data.load(in);
    _context_tables.load(in);
}

#endif // MOVE_R_MOVE_RB_MOVE_RB_TPP
