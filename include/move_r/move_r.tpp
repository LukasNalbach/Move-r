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

#ifndef MOVE_R_MOVE_R_MOVE_R_TPP
#define MOVE_R_MOVE_R_MOVE_R_TPP

#include <move_r/move_r.hpp>
#include <sstream>

// ############################# INTERNAL METHODS #############################

template <move_r_support support, move_pos_encoding_t mlf_enc>
template <typename move_r_t>
void move_r<support, mlf_enc>::assign_from(move_r_t&& other)
{
    n = other.n; sigma = other.sigma; r = other.r; r_ = other.r_; r__ = other.r__; r___ = other.r___;
    a = other.a; omega_idx = other.omega_idx;
    symbols_remapped = other.symbols_remapped;
    byte_alphabet = other.byte_alphabet;
    _map_int = std::forward<move_r_t>(other)._map_int; _map_ext = std::forward<move_r_t>(other)._map_ext;
    _alphabet = std::forward<move_r_t>(other)._alphabet;
    _alphabet_has_max = other._alphabet_has_max;
    _alphabet_order = std::forward<move_r_t>(other)._alphabet_order;
    _alphabet_order_inv = std::forward<move_r_t>(other)._alphabet_order_inv;
    _M_LF = std::forward<move_r_t>(other)._M_LF;
    _PSB_L_ = std::forward<move_r_t>(other)._PSB_L_;
    _PSO_L_ = std::forward<move_r_t>(other)._PSO_L_;
    _WTPS_L_ = std::forward<move_r_t>(other)._WTPS_L_;
    _L_prev_min = std::forward<move_r_t>(other)._L_prev_min;
    // PSB_L' and PSO_L' read L' through a pointer to the index (see L_reader)
    _PSB_L_.set_read(L_reader { this });
    _PSO_L_.set_read(L_reader { this });
    _M_Phi_inv = std::forward<move_r_t>(other)._M_Phi_inv; _M_Phi = std::forward<move_r_t>(other)._M_Phi;
    _SA_Phi_inv = std::forward<move_r_t>(other)._SA_Phi_inv;
    _SA_run_pos = std::forward<move_r_t>(other)._SA_run_pos; _SA_run = std::forward<move_r_t>(other)._SA_run;
    _SA_rlz = std::forward<move_r_t>(other)._SA_rlz; _seq_data = std::forward<move_r_t>(other)._seq_data;
    _revert_starts = std::forward<move_r_t>(other)._revert_starts;
    _rlbwt_from_cache = other._rlbwt_from_cache;
    _plcp_row_width = other._plcp_row_width;
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
template <typename input_t>
void move_r<support, mlf_enc>::construct_from_input(input_t& input, bool delete_input, move_r_params params)
{
    if (!params.input_file.empty()) {
        move_r_input_buffer_t<input_t> file_buffer;
        move_r_input_view view(file_buffer);
        construction(*this, view, true, std::move(params));
    } else {
        move_r_with_view(input, [&](move_r_input_view& view) {
            construction(*this, view, delete_input, std::move(params));
        });
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::assign_alphabet(const std::vector<uint64_t>& symbols, const std::vector<uint64_t>& order)
{
    _alphabet_has_max = !symbols.empty() && symbols.back() == std::numeric_limits<uint64_t>::max();
    uint64_t num_stored = symbols.size() - uint64_t(_alphabet_has_max);
    uint64_t universe = num_stored == 0 ? 1 : symbols[num_stored - 1] + 1;
    _alphabet = elias_fano_bit_vector(
        [&](uint64_t i) { return symbols[i]; }, num_stored, universe);
    set_alphabet_order(order);
}

// ############################# MISC PUBLIC METHODS #############################

template <move_r_support support, move_pos_encoding_t mlf_enc>
uint64_t move_r<support, mlf_enc>::input_length(const move_r_input_view& T, const move_r_params& params)
{
    if (params.input_file.empty()) return T.size() + 1;
    // the length of the text of a FASTA file (without headers and line breaks, see fasta_input_t)
    if (params.fasta.text_len != 0) return params.fasta.text_len + 1;
    const uint64_t bytes = std::filesystem::file_size(params.input_file);
    if (bytes % T.symbol_width() != 0)
        throw std::runtime_error("move_r: input-file size is not a multiple of its symbol width");
    return bytes / T.symbol_width() + 1;
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::set_alphabet_maps(const std::vector<uint8_t>& map_int, const std::vector<uint8_t>& map_ext)
    requires(is_bidirectional)
{
    symbols_remapped = true;
    byte_alphabet = true;

    _map_int = map_int;
    _map_ext = map_ext;
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::set_alphabet(const std::vector<uint64_t>& symbols,
    const std::vector<uint64_t>& order) requires(is_bidirectional)
{
    symbols_remapped = true;
    byte_alphabet = false;
    assign_alphabet(symbols, order);
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::set_alphabet_order(const std::vector<uint64_t>& order)
{
    if (order.empty()) {
        _alphabet_order = interleaved_bit_aligned_vectors();
        _alphabet_order_inv = interleaved_bit_aligned_vectors();
        return;
    }

    const uint8_t width = std::max<uint8_t>(1, std::bit_width(order.size()));
    _alphabet_order = interleaved_bit_aligned_vectors({ width });
    _alphabet_order_inv = interleaved_bit_aligned_vectors({ width });
    _alphabet_order.resize_no_init(order.size());
    _alphabet_order_inv.resize_no_init(order.size());

    for (uint64_t i = 0; i < order.size(); i++) {
        _alphabet_order.template set<0, uint64_t>(i, order[i]);
        _alphabet_order_inv.template set<0, uint64_t>(order[i] - 1, i + 1);
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
uint64_t move_r<support, mlf_enc>::size_in_bytes() const
{
    return sizeof(*this) +
        _M_LF.size_in_bytes() +
        map_ext_bytes() +
        l_data_structures_bytes() +
        _M_Phi.size_in_bytes() +
        _M_Phi_inv.size_in_bytes() +
        _SA_Phi_inv.size_in_bytes() +
        _SA_run_pos.size_in_bytes() +
        _SA_run.size_in_bytes() +
        _SA_rlz.size_in_bytes() +
        _seq_data.size_in_bytes() +
        _revert_starts.size_in_bytes();
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::log_data_structures() const
{
    if (n > max_logged_input_length) {
        std::cout << "cannot log the contents for n > " << max_logged_input_length << std::endl;
        return;
    }

    std::cout << "n = " << format_number(n) << std::endl;
    std::cout << "sigma = " << sigma << std::endl;
    if (M_LF().num_intervals() != 0) std::cout << "a = " << a << std::endl;

    // formats an external symbol, printing the terminator (0) as '$' and a symbol of a byte alphabet as a character
    auto external_sym = [&](uint64_t sym) {
        std::ostringstream s;
        if (sym == 0) s << '$';
        else if (byte_alphabet) s << char(sym);
        else s << sym;
        return s.str();
    };

    std::cout << std::endl;
    {
        auto T = revert<no_init_vector<uint64_t>>();

        // capture T's symbols (followed by the terminator) before T is remapped in place below
        std::vector<std::string> T_syms;
        T_syms.reserve(n);
        for (auto sym : T) T_syms.push_back(external_sym(sym));
        T_syms.push_back("$");

        for (uint64_t i = 0; i < T.size(); i++) {
            T[i] = map_symbol(T[i]);
        }

        auto SA = move_r<_locate_move>(T).SA_range();
        auto L = BWT_range<no_init_vector<uint64_t>>();

        aligned_table log;
        log.add_row("T:", T_syms.size(), [&](uint64_t i) { return T_syms[i]; });
        log.add_row("SA:", SA.size(), [&](uint64_t i) { return SA[i]; });
        log.add_row("L:", L.size(), [&](uint64_t i) { return external_sym(L[i]); });
        log.print();
    }
    std::cout << std::endl;

    if (M_LF().num_intervals() != 0) {
        std::cout << "M_LF: " << "r' = " << r_ << std::endl;

        // L' is logged as the first row of the M_LF block, aligned above its p, q and idx arrays
        aligned_table log;
        log.add_row("L':", r_, [&](uint64_t i) { return logged_symbol(L_(i)); });
        log.add_row("p:", r_, [&](uint64_t i) { return M_LF().p(i); });
        log.add_row("q:", r_, [&](uint64_t i) { return M_LF().q(i); });
        log.add_row("idx:", r_, [&](uint64_t i) { return M_LF().idx(i); });
        log.print();
    }

    // M_Phi^{-1} and M_Phi are empty, if they have been released (see move_rb::release_phi_fwd)
    if constexpr (has_locate_move) {
        if (!_M_Phi_inv.empty()) {
            std::cout << std::endl;
            std::cout << "M_Phi^{-1}: " << "r'' = " << r__ << std::endl;
            M_Phi_inv().log_data_structures();
        }
    }

    if constexpr (has_mphi) {
        if (!_M_Phi.empty()) {
            std::cout << std::endl;
            std::cout << "M_Phi: " << "r''' = " << r___ << std::endl;
            M_Phi().log_data_structures();
        }
    }

    // for a bidirectional (move_rb) sub-index, the run samples are logged by move_rb instead
    if constexpr (!is_bidirectional && supports_locate) {
        // n is logged for the input intervals of M_LF that do not start a BWT run
        std::cout << std::endl;
        aligned_table log;
        log.add_row("SA_run_start:", r_, [&](uint64_t x) { return run_starts_at(x) ? SA_run_start(x) : n; });
        log.print();
    }

    if (!_SA_Phi_inv.empty()) {
        std::cout << std::endl;
        log_indexed("SA_Phi^{-1}:", _SA_Phi_inv);
    }

    if (!_SA_rlz.empty()) {
        std::cout << std::endl;
        _SA_rlz.log_data_structures();
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::log_data_structure_sizes() const
{
    const bool owner = begin_listing(size_in_bytes());
    uint64_t size_l_ = (_M_LF.template row_width<0>() * (r_ + 1)) / 8;
    log_component("M_LF", _M_LF.size_in_bytes() - size_l_);
    log_component("L'", size_l_);
    if (!_PSB_L_.empty()) log_component("PSB_L'", _PSB_L_.size_in_bytes());
    if (!_PSO_L_.empty()) log_component("PSO_L'", _PSO_L_.size_in_bytes());
    if (!_WTPS_L_.empty()) log_component("WTPS_L'", _WTPS_L_.size_in_bytes());
    if (!_L_prev_min.empty()) log_component("L_prev_min", _L_prev_min.size_in_bytes());

    if (!byte_alphabet && symbols_remapped) {
        log_component("alphabet", _alphabet.size_in_bytes() + _alphabet_order.size_in_bytes()
            + _alphabet_order_inv.size_in_bytes());
    }

    if (!_M_Phi.empty()) log_component("M_Phi", _M_Phi.size_in_bytes());
    if (!_M_Phi_inv.empty()) log_component("M_Phi^{-1}", _M_Phi_inv.size_in_bytes());
    if (!_SA_Phi_inv.empty()) log_component("SA_Phi^{-1}", _SA_Phi_inv.size_in_bytes());
    if (!_SA_run.empty()) log_component("SA_run", _SA_run.size_in_bytes());
    if (!_SA_run_pos.empty()) log_component("SA_run_pos", _SA_run_pos.size_in_bytes());
    if (!_SA_rlz.empty()) _SA_rlz.log_data_structure_sizes();
    if (_seq_data.has_sequences()) log_component("FASTA sequences", _seq_data.size_in_bytes());
    if (!_revert_starts.empty()) log_component("revert starts", _revert_starts.size_in_bytes());
    if (owner) end_listing();
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::add_parameters(index_facts& facts, const std::string& suffix) const
{
    if (suffix.empty()) facts.add("n", n);
    facts.add("r" + suffix, r);
    facts.add("n/r" + suffix + " = " + format_ratio(r == 0 ? 0.0 : double(n) / double(r)));
    if (suffix.empty()) facts.add("sigma", sigma);
    facts.add("k' (M_LF" + suffix + ")", _M_LF.num_intervals());
    if (!_M_Phi.empty()) facts.add("k' (M_Phi" + suffix + ")", _M_Phi.num_intervals());
    if (!_M_Phi_inv.empty()) facts.add("k' (M_Phi^{-1}" + suffix + ")", _M_Phi_inv.num_intervals());
    if (!_SA_rlz.empty()) facts.add("z(SA_rlz" + suffix + ")", _SA_rlz.num_phrases());

    if (_seq_data.has_sequences()) facts.add("FASTA sequences", _seq_data.num_sequences());
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::log_data_structure_sizes(std::ostream& out, std::string suffix) const
{
    out << " size_index" << suffix << "=" << size_in_bytes();
    out << " a" << suffix << "=" << a;

    uint64_t size_l_ = (_M_LF.template row_width<0>() * (r_ + 1)) / 8;
    out << " size_m_lf" << suffix << "=" << _M_LF.size_in_bytes() - size_l_;
    out << " size_l_" << suffix << "=" << size_l_;

    out << " size_l_data_structures_" << suffix << "=" << l_data_structures_bytes();
    out << " size_revert_starts" << suffix << "=" << _revert_starts.size_in_bytes();
    out << " num_revert_starts" << suffix << "=" << _revert_starts.size();
    out << " size_L_prev_min" << suffix << "=" << _L_prev_min.size_in_bytes();

    if (!byte_alphabet && symbols_remapped) {
        out << " size_alphabet" << suffix << "=" << (_alphabet.size_in_bytes()
            + _alphabet_order.size_in_bytes() + _alphabet_order_inv.size_in_bytes());
        out << " size_map_ext" << suffix << "=" << map_ext_bytes();
    }

    if constexpr (support == _locate_one || support == _locate_rlzsa || support == _locate_rlzsa_bi ||
        stores_only_SA_samples) {
        out << " size_sa_run" << suffix << "=" << _SA_run.size_in_bytes() + _SA_run_pos.size_in_bytes();
    }

    if constexpr (has_locate_move) {
        if constexpr (has_mphi) {
            out << " size_m_phi" << suffix << "=" << _M_Phi.size_in_bytes();
        }

        out << " size_m_phi_inv" << suffix << "=" << _M_Phi_inv.size_in_bytes();
        out << " size_sa_phi_inv" << suffix << "=" << _SA_Phi_inv.size_in_bytes();
    }
    if constexpr (has_SA_rlz) {
        _SA_rlz.log_data_structure_sizes(out, suffix);
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
std::pair<std::vector<uint64_t>, std::vector<uint64_t>> move_r<support, mlf_enc>::split_into_run_parts(
    uint64_t num_parts, uint64_t p) const
{
    std::vector<uint64_t> part_begin(num_parts + 1, r_);
    part_begin[0] = 0;

    // part t starts at the first run start at or after the input interval r' * t / num_parts
    for (uint64_t t = 1; t < num_parts; t++) {
        uint64_t x = std::max(part_begin[t - 1], r_ * t / num_parts);
        while (x < r_ && !run_starts_at(x)) x++;
        part_begin[t] = x;
    }

    // count the runs of each part, then compute the prefix sums
    std::vector<uint64_t> part_run(num_parts + 1, 0);

    #pragma omp parallel for num_threads(p) schedule(dynamic, 1)
    for (uint64_t t = 0; t < num_parts; t++) {
        uint64_t runs = 0;
        for (uint64_t x = part_begin[t]; x < part_begin[t + 1]; x++) if (run_starts_at(x)) runs++;
        part_run[t + 1] = runs;
    }

    for (uint64_t t = 0; t < num_parts; t++) part_run[t + 1] += part_run[t];
    return { std::move(part_begin), std::move(part_run) };
}

// ############################# ACCESS METHODS #############################

template <move_r_support support, move_pos_encoding_t mlf_enc>
inline uint64_t move_r<support, mlf_enc>::SA_run_start(uint64_t x) const
    requires(builds_SA_samples)
{
    if constexpr (support == _locate_move) {
        return M_Phi_inv().q(SA_Phi_inv(x));
    } else {
        // the samples are stored in SA_Phi^{-1}, if SA_run_pos is empty
        if constexpr (has_locate_move) {
            if (_SA_run_pos.empty()) return M_Phi_inv().q(_SA_Phi_inv.template get<0, uint64_t>(x));
        }

        return _SA_run.template get<0, uint64_t>(_SA_run_pos.start_sample_index(x));
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
inline uint64_t move_r<support, mlf_enc>::SA_run_end(uint64_t x) const
    requires(support != _count && support != _count_bi && support != _locate_rlzsa && support != _locate_one)
{
    if constexpr (support == _locate_move) {
        if (x == r_ - 1) [[unlikely]] {
            return M_Phi_inv().p(SA_Phi_inv(0));
        } else {
            return M_Phi_inv().p(SA_Phi_inv(x + 1));
        }
    } else {
        // the input interval SA_Phi^{-1}[(x + 1) mod r'] of M_Phi^{-1} starts with SA[M_LF.p[x + 1] - 1]
        if constexpr (has_locate_move) {
            if (_SA_run_pos.empty())
                return M_Phi_inv().p(_SA_Phi_inv.template get<0, uint64_t>(x + 1 == r_ ? 0 : x + 1));
        }

        return _SA_run.template get<0, uint64_t>(_SA_run_pos.end_sample_index(x));
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
template <typename sym_t>
inline uint64_t move_r<support, mlf_enc>::map_symbol(sym_t sym) const
{
    if (symbol_idx(sym) == TERMINATOR_SYMBOL) [[unlikely]] return 0;
    if (byte_alphabet) return _map_int[symbol_idx(sym)];
    if (!symbols_remapped) return symbol_idx(sym);

    uint64_t value = symbol_idx(sym);
    if (value == std::numeric_limits<uint64_t>::max()) {
        // UINT64_MAX is not marked in _alphabet and has the largest rank (see _alphabet_has_max)
        if (!_alphabet_has_max) return 0;
        const uint64_t rank = _alphabet.num_ones() + 1;
        return alphabet_reordered()
            ? _alphabet_order.template get<0, uint64_t>(rank - 1) : rank;
    }

    // rank of the symbol among the symbols marked in _alphabet (0 is returned if it is not marked)
    uint64_t rank = _alphabet.rank_1(value + 1);
    if (rank == 0 || _alphabet.select_1(rank) != value) return 0;
    return alphabet_reordered() ? _alphabet_order.template get<0, uint64_t>(rank - 1) : rank;
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
bool move_r<support, mlf_enc>::alphabet_in_input_order() const
{
    if (alphabet_reordered()) return false;
    if (!byte_alphabet) return true;
    // a byte alphabet is in the order of the input iff the internal symbols of the occurring bytes increase with them
    uint64_t last = 0;

    for (uint64_t b = 1; b < 256; b++) {
        const uint64_t sym = map_symbol(uint8_t(b));
        if (sym == 0) continue;
        if (sym < last) return false;
        last = sym;
    }

    return true;
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
inline uint64_t move_r<support, mlf_enc>::unmap_symbol(uint64_t sym) const
{
    if (sym == 0) [[unlikely]] return TERMINATOR_SYMBOL;
    if (byte_alphabet) return _map_ext[sym];
    if (!symbols_remapped) return sym;
    uint64_t rank = sym;

    if (alphabet_reordered()) {
        // UINT64_MAX is returned for a value that is not an internal symbol
        if (sym > _alphabet_order_inv.size()) return std::numeric_limits<uint64_t>::max();
        rank = _alphabet_order_inv.template get<0, uint64_t>(sym - 1);
    }

    return rank <= _alphabet.num_ones()
        ? _alphabet.select_1(rank) : std::numeric_limits<uint64_t>::max();
}

// ############################# SERIALIZATION METHODS #############################

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::serialize_SA_samples(std::ostream& out) const
{
    const sa_samples_layout_t layout = support != _locate_move && !_SA_Phi_inv.empty()
        ? _as_phi_inv_intervals : _at_run_boundaries;
    out.write((char*) &layout, 1);

    if (layout == _at_run_boundaries) {
        _SA_run_pos.serialize(out);
        _SA_run.serialize(out);
    } else {
        out.write((char*) &omega_idx, 1);
        _SA_Phi_inv.serialize(out);
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::read_SA_samples(std::istream& in)
{
    uint8_t layout = 0;
    in.read((char*) &layout, 1);

    if (layout == _at_run_boundaries) {
        _SA_run_pos.load(in);
        _SA_run.load(in);
    } else if (layout == _as_phi_inv_intervals) {
        in.read((char*) &omega_idx, 1);
        _SA_Phi_inv.load(in);
    } else {
        throw std::runtime_error("move_r: unsupported format of the suffix array samples, "
            "the index has to be rebuilt");
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
move_r_support move_r<support, mlf_enc>::peek_support(std::istream& in)
{
    std::streampos pos = in.tellg();
    move_r_support _support;
    in.read((char*) &_support, sizeof(move_r_support));
    in.seekg(pos);
    return _support;
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::serialize(std::ostream& out) const
{
    move_r_support _support = support;
    out.write((char*) &_support, sizeof(move_r_support));

    out.write((char*) &n, sizeof(uint64_t));
    out.write((char*) &sigma, sizeof(uint32_t));
    out.write((char*) &r, sizeof(uint64_t));
    out.write((char*) &r_, sizeof(uint64_t));
    uint16_t a_ = a;
    out.write((char*) &a_, sizeof(uint16_t));
    out.write((char*) &byte_alphabet, 1);
    out.write((char*) &symbols_remapped, 1);
    if (symbols_remapped) {
        if (byte_alphabet) {
            out.write((char*) &_map_int[0], 256);
            out.write((char*) &_map_ext[0], 256);
        } else {
            _alphabet.serialize(out);
            out.write((char*) &_alphabet_has_max, 1);
            _alphabet_order.serialize(out);
            _alphabet_order_inv.serialize(out);
        }
    }

    _M_LF.serialize(out);
    // the data structures over L' are also serialized if they are empty
    _PSB_L_.serialize(out);
    _PSO_L_.serialize(out);
    _WTPS_L_.serialize(out);
    _L_prev_min.serialize(out);

    if constexpr (support == _locate_one) {
        serialize_SA_samples(out);
    } else if constexpr (has_locate_move) {
        if constexpr (has_mphi) {
            out.write((char*) &r___, sizeof(uint64_t));
            _M_Phi.serialize(out);
        }

        out.write((char*) &r__, sizeof(uint64_t));
        _M_Phi_inv.serialize(out);

        if constexpr (support == _locate_move) {
            out.write((char*) &omega_idx, 1);
            _SA_Phi_inv.serialize(out);
        } else {
            serialize_SA_samples(out);
        }
    } else if constexpr (has_SA_rlz) {
        serialize_SA_samples(out);
        _SA_rlz.serialize(out);
    } else if constexpr (stores_only_SA_samples) {
        serialize_SA_samples(out);
    }

    _seq_data.serialize(out);
    _revert_starts.serialize(out);
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::load(std::istream& in)
{
    move_r_support _support;
    in.read((char*) &_support, sizeof(move_r_support));

    if (_support != support) {
        throw std::runtime_error("move_r: the index was built with another locate support (" +
            std::to_string(int(_support)) + ") than the index type it is loaded into (" +
            std::to_string(int(support)) + ")");
    }

    in.read((char*) &n, sizeof(uint64_t));
    in.read((char*) &sigma, sizeof(uint32_t));
    in.read((char*) &r, sizeof(uint64_t));
    in.read((char*) &r_, sizeof(uint64_t));
    uint16_t a_;
    in.read((char*) &a_, sizeof(uint16_t));
    a = a_;
    in.read((char*) &byte_alphabet, 1);
    in.read((char*) &symbols_remapped, 1);
    if (symbols_remapped) {
        if (byte_alphabet) {
            _map_int.resize(256);
            in.read((char*) &_map_int[0], 256);

            _map_ext.resize(256);
            in.read((char*) &_map_ext[0], 256);
        } else {
            _alphabet.load(in);
            in.read((char*) &_alphabet_has_max, 1);
            _alphabet_order.load(in);
            _alphabet_order_inv.load(in);
        }
    }

    _M_LF.load(in, mlf_num_rows);
    _PSB_L_.load(in);
    _PSB_L_.set_read(L_reader { this });
    _PSO_L_.load(in);
    _PSO_L_.set_read(L_reader { this });
    _WTPS_L_.load(in);
    _L_prev_min.load(in);

    if constexpr (support == _locate_one) {
        read_SA_samples(in);
    } else if constexpr (has_locate_move)
    {
        if constexpr (has_mphi) {
            in.read((char*) &r___, sizeof(uint64_t));
            _M_Phi.load(in, 0);
        }

        in.read((char*) &r__, sizeof(uint64_t));
        _M_Phi_inv.load(in, mphi_inv_num_rows);

        if constexpr (support == _locate_move) {
            in.read((char*) &omega_idx, 1);
            _SA_Phi_inv.load(in);
        } else {
            read_SA_samples(in);
        }
    } else if constexpr (has_SA_rlz)
    {
        read_SA_samples(in);
        _SA_rlz.load(in);
    } else if constexpr (stores_only_SA_samples) {
        read_SA_samples(in);
    }

    _seq_data.load(in);
    _revert_starts.load(in);
}

#endif // MOVE_R_MOVE_R_MOVE_R_TPP
