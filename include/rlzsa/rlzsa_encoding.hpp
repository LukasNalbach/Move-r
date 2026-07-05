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

#include <misc/utils.hpp>
#include <misc/files.hpp>
#include <misc/search.hpp>
#include <misc/log.hpp>

#include <ips4o.hpp>
#include <libsais.h>
#include <libsais64.h>
#include <ankerl/unordered_dense.h>
#include <sdsl/int_vector.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/csa_bitcompressed.hpp>
#include <sdsl/suffix_array_algorithm.hpp>

/**
 * @brief a standalone relative-Lempel-Ziv encoding of the differential suffix array (unoptimized; the move-r
 *        optimized variant is rlzsa_opt in include/rlzsa_opt/)
 * @tparam int_t signed integer type of the suffix array entries
 */
template <typename int_t = int32_t>
class rlzsa_encoding {

protected:
    static constexpr uint64_t PS_sample_rate = 64;

    uint64_t n; // input size
    uint64_t z; // number of rlzsa phrases

    sdsl::int_vector<> R; // rlzsa reference
    sdsl::int_vector<> PS; // rlzsa phrase starting positions (every PS_sample_rate-th phrase is sampled)
    sdsl::int_vector<16> PL; // rlzsa phrase lengths
    sdsl::int_vector<> S; // starting positions of the rlzsa phrases in R

public:
    rlzsa_encoding() = default;

    /**
     * @brief builds the rlzsa encoding from the suffix array and the differential suffix array
     * @param SA the suffix array
     * @param DSA the differential suffix array (consumed during the construction)
     * @param reference_size the target size of the reference R
     * @param add_literal_phrase_after_each_copy_phrase whether to add a literal phrase after each copy phrase
     * @param log whether to print log messages
     */
    rlzsa_encoding(
        const std::vector<int_t>& SA,
        std::vector<int_t>&& DSA,
        uint64_t reference_size,
        bool add_literal_phrase_after_each_copy_phrase = false,
        bool log = false
    ) {
        n = DSA.size();
        auto time = now();
        log_phase_start(log, time, "sorting DSA");

        for (uint64_t i = 1; i < n; i++) {
            DSA[i] += n;
        }

        std::vector<int_t> alphabet; // alphabet of DSA
        alphabet.resize(n);

        for (uint64_t i = 0; i < n; i++) {
            alphabet[i] = DSA[i];
        }

        ips4o::sort(alphabet.begin(), alphabet.end()); // sort the values in DSA

        log_phase_end(log, time);
        log_phase_start(log, time, "building DSA alphabet map");

        auto it = std::unique(alphabet.begin(), alphabet.end()); // remove duplicates in DSA
        alphabet.resize(std::distance(alphabet.begin(), it));
        ankerl::unordered_dense::map<int_t, int_t> alphabet_map; // map to the effective alphabet of DSA
        sdsl::int_vector<> alphabet_map_inv; // map from the effective alphabet of DSA
        alphabet_map_inv.width(8 * byte_width(2 * n));
        alphabet_map_inv.resize(alphabet.size() + 1);
        alphabet_map_inv[0] = 0;

        {
            uint64_t j = 1;

            // build alphabet_map and alphabet_map_inv
            for (uint64_t i = 0; i < alphabet.size(); i++) {
                alphabet_map.insert({ alphabet[i], j });
                alphabet_map_inv[j] = alphabet[i];
                j++;
            }
        }

        log_phase_end(log, time);
        log_phase_start(log, time, "mapping DSA to its effective alphabet");

        // map DSA to its effective alphabet
        for (uint64_t i = 0; i < n; i++) {
            DSA[i] = alphabet_map[DSA[i]];
        }

        log_phase_end(log, time);
        log_phase_start(log, time, "building suffix array of DSA");

        std::vector<int_t> sa_dsa;
        sa_dsa.resize(n);

        if constexpr (std::is_same_v<int_t, int32_t>) {
            libsais_int(DSA.data(), sa_dsa.data(), n, alphabet.size() + 1, 0);
        } else {
            libsais64_long(DSA.data(), sa_dsa.data(), n, alphabet.size() + 1, 0);
        }

        sa_dsa.resize(n);
        log_phase_end(log, time);
        log_phase_start(log, time, "building inverse suffix array of DSA");

        sdsl::int_vector<> ISA_sad;
        ISA_sad.width(8 * byte_width(n));
        ISA_sad.resize(n);

        for (uint64_t i = 0; i < n; i++) {
            ISA_sad[sa_dsa[i]] = i;
        }

        log_phase_end(log, time);
        log_phase_start(log, time, "counting frequencies of k-mers in DSA");

        uint64_t segment_size = std::max<uint64_t>(1, std::min<uint64_t>(4096, reference_size)); // >= 1 to avoid a division by zero for tiny inputs
        uint64_t num_segments = n / segment_size;
        std::vector<uint64_t> selected_segments;
        uint64_t k = std::min<uint64_t>(8, reference_size);
        uint64_t num_groups;

        sdsl::int_vector<> group;
        sdsl::int_vector<> freq;

        group.width(8 * byte_width(n));
        freq.width(8 * byte_width(n));

        group.resize(n);
        freq.resize(n);

        {
            uint64_t i = 1;
            uint64_t g = 0;

            while (i < n) {
                uint64_t j = 1;
                bool equal = true;

                while (equal && i + j < n) {
                    for (uint64_t l = 0; l < k; l++) {
                        if (sa_dsa[i + j] + l == n ||
                            sa_dsa[i + j - 1] + l == n ||
                            DSA[sa_dsa[i + j] + l] != DSA[sa_dsa[i + j - 1] + l]
                        ) {
                            equal = false;
                            break;
                        }
                    }

                    if (equal) {
                        j++;
                    }
                }

                for (uint64_t p = 0; p < j; p++) {
                    group[i + p] = g;
                    freq[i + p] = j;
                }

                g++;
                i += j;
            }

            group[0] = g;
            freq[0] = 1;
            num_groups = g + 1;
        }

        log_phase_end(log, time);
        log_phase_start(log, time, "counting frequencies of k-mers per segment of DSA");

        {
            uint64_t num_segments_to_select = reference_size / segment_size;
            std::vector<double> segment_freqs;
            segment_freqs.resize(num_segments + 1, 0);

            for (uint64_t s = 0; s < num_segments; s++) {
                sdsl::bit_vector is_group_processed(num_groups + 1, 0);

                for (uint64_t i = s * segment_size; i < (s + 1) * segment_size - k; i++) {
                    uint64_t p = ISA_sad[i];

                    if (is_group_processed[group[p]] == 0) {
                        segment_freqs[s] += std::sqrt((double) freq[p]);
                        is_group_processed[group[p]] = 1;
                    }
                }
            }

            log_phase_end(log, time);
            log_phase_start(log, time, "building R");

            selected_segments.resize(num_segments_to_select);
            sdsl::bit_vector is_segment_selected(num_segments + 1, 0);
            sdsl::bit_vector is_group_processed(num_groups + 1, 0);

            for (uint64_t s = 0; s < num_segments_to_select; s++) {
                uint64_t best_segment = 0;

                for (uint64_t i = 1; i < num_segments; i++) {
                    if (is_segment_selected[i] == 0 && segment_freqs[i] > segment_freqs[best_segment]) {
                        best_segment = i;
                    }
                }

                selected_segments[s] = best_segment;
                is_segment_selected[best_segment] = 1;
                segment_freqs[best_segment] = 0;

                for (uint64_t i = best_segment * segment_size; i < (best_segment + 1) * segment_size - k; i++) {
                    uint64_t p = ISA_sad[i];

                    if (p < n && is_group_processed[group[p]] == 0) {
                        uint64_t g = group[p];
                        double freq_sqrt = std::sqrt((double) freq[p]);

                        while (p > 0 && group[p - 1] == g) {
                            p--;
                        }

                        sdsl::bit_vector is_segment_processed(num_segments + 1, 0);
                        is_segment_processed[best_segment] = 1;

                        do {
                            uint64_t sa_p = sa_dsa[p];
                            uint64_t seg = sa_p / segment_size;

                            if (sa_p < (seg + 1) * segment_size - k &&
                                is_segment_processed[seg] == 0 &&
                                is_segment_selected[seg] == 0
                            ) {
                                segment_freqs[seg] -= freq_sqrt;
                                is_segment_processed[seg] = 1;
                            }

                            p++;
                        } while (p < n && group[p] == g);

                        is_group_processed[g] = 1;
                    }
                }
            }

            ips4o::sort(selected_segments.begin(), selected_segments.end());
        }

        {
            R.width(8 * byte_width(2 * n));
            R.resize(segment_size * selected_segments.size());
            uint64_t j = 0;

            for (uint64_t i = 0; i < selected_segments.size(); i++) {
                uint64_t segment_start = selected_segments[i] * segment_size;

                for (uint64_t l = 0; l < segment_size; l++) {
                    R[j] = alphabet_map_inv[DSA[segment_start + l]];
                    j++;
                }
            }
        }

        selected_segments.clear();
        selected_segments.shrink_to_fit();

        selected_segments.clear();
        selected_segments.shrink_to_fit();

        group.resize(0);
        freq.resize(0);

        ISA_sad.resize(0);

        log_phase_end(log, time);
        log_phase_start(log, time, "building FM-Index for R");

        sdsl::int_vector<> revR;
        revR.width(8 * byte_width(2 * n));
        revR.resize(R.size());

        for (uint64_t i = 0; i < R.size(); i++) {
            revR[i] = R[R.size() - i - 1];
        }

        sdsl::csa_bitcompressed<sdsl::int_alphabet<>> FM_rrev;
        sdsl::construct_im(FM_rrev, revR, 0);
        revR.resize(0);

        log_phase_end(log, time);
        log_phase_start(log, time, "building rlzsa parsing");

        S.width(8 * byte_width(n));
        PS.width(8 * byte_width(n));

        S.resize(1);
        PS.resize(1);
        PL.resize(1);
        S[0] = SA[0];
        PS[0] = 0;
        PL[0] = 0;
        uint64_t z_l = 1;

        {
            z = 1;
            uint64_t i = 1;
            sdsl::int_vector<> DSA_val;
            DSA_val.resize(1);

            // rlz-encode DSA
            while (i < n) {
                uint64_t b = 0;
                uint64_t e = FM_rrev.size() - 1;
                uint64_t b_last;
                uint64_t l = 0;

                // find the longest prefix DSA[i,i+d) of DSA[i,n) that occurs in R
                while (true) {
                    b_last = b;
                    DSA_val[0] = alphabet_map_inv[DSA[i + l]];

                    if (sdsl::backward_search(FM_rrev, b, e, DSA_val.begin(), DSA_val.end(), b, e) == 0) {
                        break;
                    }

                    l++;

                    if (l == UINT16_MAX || i + l == n) {
                        b_last = b;
                        break;
                    }
                }

                z++;
                PL.resize(z);
                PL[z - 1] = l;
                S.resize(z);

                if (z % PS_sample_rate == 1) {
                    PS.resize(PS.size() + 1);
                    PS[PS.size() - 1] = i;
                }

                if (l == 0) {
                    // literal phrase
                    S[z - 1] = SA[i];
                    i++;
                    z_l++;
                } else {
                    // copy-phrase
                    S[z - 1] = (R.size() - FM_rrev[b_last]) - l;
                    i += l;

                    if (add_literal_phrase_after_each_copy_phrase && i < n) {
                        // add a literal phrase after each copy-phrase
                        z++;
                        PL.resize(z);
                        PL[z - 1] = 0;
                        S.resize(z);
                        S[z - 1] = SA[i];

                        if (z % PS_sample_rate == 1) {
                            PS.resize(PS.size() + 1);
                            PS[PS.size() - 1] = i;
                        }

                        i++;
                        z_l++;
                    }
                }
            }
            
            PL.resize(z+1);
			PL[z] = 1;
        }
        
        log_phase_end(log, time);
        if (log) {
            std::cout << "z: " << z << ", ";
            std::cout << "n/z: " << n / (double) z << ", ";
            std::cout << "z_l/z: " << z_l / (double) z << std::endl;
            std::cout << "peak memory usage: " << format_size(malloc_count_peak()) << std::endl;
        }
    }
    
    /**
     * @brief returns the number of phrases in the rlzsa
     * @return the number of phrases
     */
    uint64_t num_phrases() const
    {
        return z;
    }

    /**
     * @brief returns the size of the input
     * @return the size of the input
     */
    uint64_t input_size() const
    {
        return n;
    }

    /**
     * @brief extracts the suffix array values in the range [l, r], reporting each value
     * @tparam out_t output value type
     * @tparam report_fnc_t type of the report function
     * @param l left range limit
     * @param r right range limit
     * @param report function that is called with every extracted value (or (position, value))
     * @param sa_l SA[l], if already known (-1 if it has to be decoded; default: -1)
     */
    template <typename out_t, typename report_fnc_t>
    void extract(uint64_t l, uint64_t r, report_fnc_t report, int64_t sa_l = -1) const
    {
        static constexpr bool report_pos = function_traits<report_fnc_t>::arity > 1;

        if (r < l) return;

        // index in PS of the last phrase starting before or at position l		
        uint64_t x_ps = bin_search_max_leq<uint64_t>(l, 0, PS.size() - 1, [&](uint64_t x) { return PS[x]; });
        uint64_t s_cp = PS[x_ps]; // starting position of the current phrase
        uint64_t x_p = PS_sample_rate * x_ps; // index of the current phrase
        
        // find the phrase containing l
        while (s_cp + std::max<uint64_t>(1,PL[x_p]) <= l) {
            s_cp += std::max<uint64_t>(1,PL[x_p]);
            x_p++;
        }

        uint64_t i; // current position in the suffix array
        uint64_t s; // current suffix

        if (sa_l == -1) {
            if (PL[x_p] == 0) {
                // the x_p-th phrase is a literal phrase
                i = l;
            } else {
                // the x_p-th phrase is a copy phrase, hence the x_p-1-th phrase is a
                // literal phrase, because there is a literal phrase before each copy-phrase;
                // so store its SA-value in s and decode the x_p-th phrase up to position l
    
                i = s_cp;
                s = S[x_p - 1];
                uint64_t p_r = S[x_p]; // current position in R
    
                while (i < l) {
                    s += R[p_r];
                    s -= n;
                    p_r++;
                    i++;
                }
            }
        } else {
            s = sa_l;
            i = l;
            if constexpr (report_pos) report(i, s); else report(s);
            i++;

            if (i == s_cp + std::max<uint64_t>(1, PL[x_p])) {
                x_p++;
                s_cp = i;
            }
        }

        uint64_t s_np = s_cp + std::max<uint64_t>(1, PL[x_p]); // starting position of the next phrase
    
        while (true) {
            // decode all literal phrases until the next copy phrase or until i > r
            while (PL[x_p] == 0 && i <= r) {
                s = S[x_p];
                if constexpr (report_pos) report(i, s); else report(s);
                x_p++;
                s_cp++;
                s_np += std::max<uint64_t>(1, PL[x_p]);
                i++;
            }

            if (i > r) {
                break;
            }

            uint64_t p_r = S[x_p] + (i - s_cp); // current position in R

            // decode the copy phrase and stop as soon as i > r
            do {
                s += R[p_r];
                s -= n;
                if constexpr (report_pos) report(i, s); else report(s);
                p_r++;
                i++;
            } while (i < s_np && i <= r);

            if (i > r) {
                break;
            }

            x_p++;
            s_cp = s_np;
            s_np += std::max<uint64_t>(1, PL[x_p]);
        }
    }

    /**
     * @brief extracts and returns the suffix array values in the range [l, r]
     * @tparam out_t output value type
     * @param l left range limit
     * @param r right range limit
     * @param sa_l SA[l], if already known (-1 if it has to be decoded; default: -1)
     * @return the suffix array values in the range [l, r]
     */
    template <typename out_t>
    std::vector<out_t> extract(uint64_t l, uint64_t r, int64_t sa_l = -1) const
    {
        std::vector<out_t> result;
        result.reserve(r - l + 1);
        extract<out_t>(l, r, [&](out_t v){result.emplace_back(v);}, sa_l);
        return result;
    }

    /**
     * @brief serializes the encoding to an output stream
     * @param out output stream
     */
    void serialize(std::ostream& out) const
    {
        out.write((char*) &n, sizeof(uint64_t));
        out.write((char*) &z, sizeof(uint64_t));

        R.serialize(out);
        PS.serialize(out);
        PL.serialize(out);
        S.serialize(out);
    }

    /**
     * @brief loads the encoding from an input stream
     * @param in input stream
     */
    void load(std::istream& in)
    {
        in.read((char*) &n, sizeof(uint64_t));
        in.read((char*) &z, sizeof(uint64_t));

        R.load(in);
        PS.load(in);
        PL.load(in);
        S.load(in);
    }

    /**
     * @brief returns the size of the data structure in bytes
     * @return the size of the data structure in bytes
     */
    uint64_t size_in_bytes() const
    {
        uint64_t size = 0;

        size += sdsl::size_in_bytes(R);
        size += sdsl::size_in_bytes(PS);
        size += sdsl::size_in_bytes(PL);
        size += sdsl::size_in_bytes(S);

        return size;
    }
};