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

#include <random>

#include <sdsl/int_vector_buffer.hpp>
#include <gtl/btree.hpp>
#include <gtl/phmap.hpp>
#include <hash_table5.hpp>

#include <move_r/move_r.hpp>
#include <rlzsa_opt/rlzsa_opt.hpp>

/** optional parameters for the rlzsa construction */
struct rlzsa_construction_params {
    uint16_t num_threads = omp_get_max_threads(); // maximum number of threads to use during the construction
    move_r_construction_mode mode = _suffix_array; // construction mode to use
    bool log = false; // controls, whether to print log messages
    /** path to the suffix-array file output by bigbwt; if non-empty, SA^d is read from this file (every
     * 40-bit entry is an SA value), otherwise SA^d is obtained from the caller-provided SAd function */
    std::string sa_file_path = "";
};

template <typename pos_t>
template <typename sad_func_t>
class rlzsa_opt<pos_t>::construction {
public:
    // ############################# CONSTRUCTION INPUTS #############################

    /** the rlzsa to build */
    rlzsa_opt<pos_t>& rlz;
    /** caller-provided function returning SA^d[i] (used when no suffix-array file is given); SAd_func(i_p, i)
     * is evaluated by thread i_p */
    sad_func_t SAd_func;
    /** The disjoint interval sequence for Phi^{-1}; consumed (cleared) while computing the value frequencies */
    std::vector<std::pair<pos_t, pos_t>>& I_Phi_m1;

    pos_t n; // the length of T
    pos_t r; // r, the number of runs in L
    uint16_t p; // the number of threads to use
    move_r_construction_mode mode; // the construction mode
    bool log; // controls, whether to print log messages
    /** path to the suffix-array file; if non-empty, SA^d is read from this file */
    std::string sa_file_path;
    /** true <=> SA^d values are read from sa_file_path (instead of from SAd_func) */
    bool from_file;

    /**
     * @brief constructs the rlzsa construction
     * @param rlz the rlzsa to build
     * @param SAd caller-provided function returning SA^d[i] (callable as SAd(i_p, i)); only used if
     *        params.sa_file_path is empty
     * @param I_Phi_m1 the disjoint interval sequence for Phi^{-1}
     * @param n the length of T
     * @param r the number of runs in L
     * @param params optional construction parameters
     */
    construction(
        rlzsa_opt<pos_t>& rlz,
        sad_func_t SAd,
        std::vector<std::pair<pos_t, pos_t>>& I_Phi_m1,
        pos_t n,
        pos_t r,
        rlzsa_construction_params params = {})
        : rlz(rlz)
        , SAd_func(std::move(SAd))
        , I_Phi_m1(I_Phi_m1)
        , n(n)
        , r(r)
        , p(params.num_threads)
        , mode(params.mode)
        , log(params.log)
        , sa_file_path(params.sa_file_path)
        , from_file(!params.sa_file_path.empty())
    {
        prefix_tmp_files = "rlzsa_" + random_alphanumeric_string(10);
        baseline_mem_usage = malloc_count_current();
    }

    // ############################# rlzsa CONSTRUCTION STATE OWNED BY THE CONSTRUCTION #############################

    /** [0..p] n_p[0] < n_p[1] < ... < n_p[p] = n; n_p[i] = start position of thread i's section in SA;
     * built by the construction from n and p */
    std::vector<pos_t> n_p;
    /** [0..p-1] file buffers (of each thread) for reading the suffix array file; only used when from_file */
    std::vector<sdsl::int_vector_buffer<40>> SA_file_bufs;
    std::string prefix_tmp_files; // prefix of the temporary files (chosen by the construction)
    std::chrono::steady_clock::time_point time; // time of the start of the last build phase
    uint64_t baseline_mem_usage = 0; // memory allocation at the start of the construction

    /**
     * @brief returns SA^d[i]; reads from the suffix-array file (if from_file) or via the SAd function
     * @param i_p [0..p-1] thread index (selects the thread's suffix-array file buffer)
     * @param i index in SA^d
     * @return SA^d[i]
     */
    inline uint64_t SAd(uint16_t i_p, pos_t i)
    {
        if (from_file) {
            if (i <= 1) [[unlikely]] {
                if (i == 0) [[unlikely]] return uint64_t { n } - 1;
                return SA_file_bufs[i_p][0] + 1;
            }
            return (SA_file_bufs[i_p][i - 1] + uint64_t { n }) - SA_file_bufs[i_p][i - 2];
        } else {
            return SAd_func(i_p, i);
        }
    }

    /**
     * @brief logs the peak memory usage until now
     */
    void log_peak_mem_usage()
    {
        std::cout << "peak memory allocation until now: "
                  << format_size(malloc_count_peak() - baseline_mem_usage)
                  << std::endl;
    }

    // ############################# rlzsa CONSTRUCTION MISC VARIABLES #############################

    pos_t size_R = 0; // size of the reference (R) in the rlzsa
    pos_t size_R_target = 0; // target size for R
    pos_t seg_size = 0; // (maximum) size of each segment from SA^d to be included in R
    pos_t num_cand_segs = 0; // number of candidate segments that are considered in each iteration during the construction of R

    // ############################# rlzsa CONSTRUCTION DATA STRUCTURE VARIABLES #############################

    /** type of hash map for storing the frequencies of values in SA^d */
    template <typename sad_t>
    using sad_freq_t = emhash5::HashMap<sad_t, pos_t, std::identity>;

    /** hashmap (for sad_t = uint32_t) mapping each value in SA^d (+n) to its frequency in SA^d */
    sad_freq_t<uint32_t> SAd_freq_32;
    /** [0..p-1] hashmaps ... (for sad_t = uint64_t; see SAd_freq_32) */
    sad_freq_t<uint64_t> SAd_freq_64;

    /** @brief a segment SA^d[beg..end) of SA^d that is a candidate for inclusion in the reference R */
    struct segment {
        pos_t beg;
        pos_t end;
    };

    /** comparator for segment in T_s */
    struct cmp_ts {
        bool operator()(const segment& s1, const segment& s2) const { return s1.beg < s2.beg; }
    };

    using ts_t = gtl::btree_set<segment, cmp_ts>; // type of T_s
    using ts_it_t = ts_t::iterator; // type of iterator in T_s

    /** B-tree storing the selected segments from SA^d to be included in the reference (R) for the rlzsa; stores pairs (pos,len)
     * to represent segments SA^d[pos..pos+len) in ascending order of their starting position (pos) in SA^d */
    ts_t T_s;

    // gap between two consecutive segments
    struct gap {
        pos_t beg_prev;
        float score;
    };

    /** comparator for gaps in T_g */
    struct cmp_tg {
        bool operator()(const gap& g1, const gap& g2) const
        {
            return g1.score > g2.score || (g1.score == g2.score && g1.beg_prev < g2.beg_prev);
        }
    };

    /** B-tree storing the gaps between the selected segments; stores pairs (beg_prev,score), where beg_prev is the starting
     * position of the segment preceding the gap, and score is the length of the connected segment resulting from closing the
     * gap divided by the length of the gap; the gaps ordered descendingly by their score */
    gtl::btree_set<gap, cmp_tg> T_g;

    /** build buffer for the reference (R) of the rlzsa; moved into rlz once the rlzsa has been built */
    interleaved_bit_aligned_vectors<uint64_t> R;

    /** [0..size_R-1] rev(R) (for sad_t = uint32_t) */
    std::vector<uint32_t> revR_32;
    /** [0..size_R-1] rev(R) (for sad_t = uint64_t) */
    std::vector<uint64_t> revR_64;

    /** type of move-r index for finding copy phrases in the rlzsa factorization */
    template <typename sad_t, typename irr_pos_t>
    using idx_revr_t = move_r<_locate_one, sad_t, irr_pos_t>;

    /** index for finding maximum length copy phrases in the rlzsa factorization (for sad_t = uint32_t and irr_pos_t = uint32_t) */
    idx_revr_t<uint32_t, uint32_t> idx_revR_32_32;
    /** index for finding maximum length copy phrases in the rlzsa factorization (for sad_t = uint64_t and irr_pos_t = uint32_t) */
    idx_revr_t<uint64_t, uint32_t> idx_revR_64_32;
    /** index for finding maximum length copy phrases in the rlzsa factorization (for sad_t = uint64_t and irr_pos_t = uint64_t) */
    idx_revr_t<uint64_t, uint64_t> idx_revR_64_64;

    /** [0..p-1] stores at position i_p in ascending order the lengths of rlzsa
     * copy phrases in SA^d[n_p[i_p]..n_p[i_p+1]-1] */
    std::vector<std::vector<uint16_t>> _CPL;
    /** [0..p-1] stores at position i_p in ascending order the starting positions (in R) of rlzsa
     * copy phrases in SA^d[n_p[i_p]..n_p[i_p+1]-1] */
    std::vector<interleaved_bit_aligned_vectors<pos_t>> _SR;
    /** [0..p-1] stores at position i_p in ascending order the literal phrases of the rlzsa
     * in SA^d[n_p[i_p]..n_p[i_p+1]-1] */
    std::vector<interleaved_bit_aligned_vectors<pos_t>> _LP;
    /** [0..p-1] stores at position i_p in ascending order the types of phrases of the rlzsa in SA^d[n_p[i_p]..n_p[i_p+1]-1];
     * i.e, PT[i_p][i] = 0 <=> the i-th phrase in SA^d[n_p[i_p]..n_p[i_p+1]-1] is literal */
    std::vector<sdsl::bit_vector> _PT;
    /** [0..p-1] stores at position i_p a file buffer with the same content as CPL[i_p]*/
    std::vector<sdsl::int_vector_buffer<>> CPL_file_bufs;
    /** [0..p-1] stores at position i_p a file buffer with the same content as SR[i_p]*/
    std::vector<sdsl::int_vector_buffer<>> SR_file_bufs;
    /** [0..p-1] stores at position i_p a file buffer with the same content as LP[i_p]*/
    std::vector<sdsl::int_vector_buffer<>> LP_file_bufs;
    /** [0..p-1] stores at position i_p a file buffer with the same content as PT[i_p]*/
    std::vector<sdsl::int_vector_buffer<>> PT_file_bufs;
    /** [0..p] z_p[0] < z_p[1] < ... < z_p[p] = z; z_p[i_p] = number of phrases in the rlzsa starting before n_p[i_p] in SA^d */
    std::vector<pos_t> z_p;
    /** [0..p] zl_p[0] < zl_p[1] < ... < zl_p[p] = z_l; zl_p[i_p] = number of literal phrases in the rlzsa starting before n_p[i_p] in SA^d */
    std::vector<pos_t> zl_p;
    /** [0..p] zc_p[0] < zc_p[1] < ... < zc_p[p] = z_c; zc_p[i_p] = number of copy phrases in the rlzsa starting before n_p[i_p] in SA^d */
    std::vector<pos_t> zc_p;

    // ############################# rlzsa CONSTRUCTION MISC METHODS #############################

    /**
     * @brief returns SAd_freq
     * @tparam sad_t type of the values in SA^d
     * @return SAd_freq
     */
    template <typename sad_t>
    constexpr sad_freq_t<sad_t>& get_SAd_freq()
    {
        if constexpr (std::is_same_v<sad_t, uint32_t>) {
            return SAd_freq_32;
        } else {
            return SAd_freq_64;
        }
    }

    /**
     * @brief returns idx_revR
     * @tparam sad_t type of the values in SA^d
     * @tparam irr_pos_t position type (pos_t) for the index of rev(R)
     * @return idx_revR
     */
    template <typename sad_t, typename irr_pos_t>
    constexpr idx_revr_t<sad_t, irr_pos_t>& get_idx_revR()
    {
        if constexpr (std::is_same_v<sad_t, uint32_t>) {
            return idx_revR_32_32;
        } else {
            if constexpr (std::is_same_v<irr_pos_t, uint32_t>) {
                return idx_revR_64_32;
            } else {
                return idx_revR_64_64;
            }
        }
    }

    /**
     * @brief returns revR
     * @tparam sad_t type of the values in SA^d
     * @return revR
     */
    template <typename sad_t>
    constexpr std::vector<sad_t>& get_revR()
    {
        if constexpr (std::is_same_v<sad_t, uint32_t>) {
            return revR_32;
        } else {
            return revR_64;
        }
    }

    /**
     * @brief returns CPL[i_p][i]
     * @tparam space controls, whether to read CPL from a file
     * @param i_p [0..p-1] thread index
     * @param i [zc_p[i_p]..zc_p[i_p+1]-1] index in CPL[i_p]
     * @return CPL[i_p][i]
     */
    template <bool space>
    inline uint16_t CPL(uint16_t i_p, pos_t i)
    {
        if constexpr (space) {
            return CPL_file_bufs[i_p][i];
        } else {
            return _CPL[i_p][i];
        }
    }

    /**
     * @brief returns SR[i_p][i]
     * @tparam space controls, whether to read SR from a file
     * @param i_p [0..p-1] thread index
     * @param i [zc_p[i_p]..zc_p[i_p+1]-1] index in SR[i_p]
     * @return SR[i_p][i]
     */
    template <bool space>
    inline pos_t SR(uint16_t i_p, pos_t i)
    {
        if constexpr (space) {
            return SR_file_bufs[i_p][i];
        } else {
            return _SR[i_p][i];
        }
    }

    /**
     * @brief returns LP[i_p][i]
     * @tparam space controls, whether to read LP from a file
     * @param i_p [0..p-1] thread index
     * @param i [zl_p[i_p]..zl_p[i_p+1]-1] index in LP[i_p]
     * @return LP[i_p][i]
     */
    template <bool space>
    inline pos_t LP(uint16_t i_p, pos_t i)
    {
        if constexpr (space) {
            return LP_file_bufs[i_p][i];
        } else {
            return _LP[i_p][i];
        }
    }

    /**
     * @brief returns PT[i_p][i]
     * @tparam space controls, whether to read PT from a file
     * @param i_p [0..p-1] thread index
     * @param i [z_p[i_p]..z_p[i_p+1]-1] index in PT[i_p]
     * @return PT[i_p][i]
     */
    template <bool space>
    inline bool PT(uint16_t i_p, pos_t i)
    {
        if constexpr (space) {
            return PT_file_bufs[i_p][i];
        } else {
            return _PT[i_p][i];
        }
    }

    // ############################# rlzsa CONSTRUCTION DISPATCHERS #############################

    /**
     * @brief builds the rlzsa
     */
    void build()
    {
        size_R_target = std::min<pos_t>(std::max<pos_t>(1, n / 3), 5.2 * r);
        seg_size = std::min<pos_t>(3072, size_R_target);

        if constexpr (std::is_same_v<pos_t, uint32_t>) {
            build<uint32_t, uint32_t>();
        } else if (2 * n <= UINT_MAX) {
            build<uint32_t, uint32_t>();
        } else if (size_R_target + seg_size <= UINT_MAX) {
            build<uint64_t, uint32_t>();
        } else {
            build<uint64_t, uint64_t>();
        }
    }

    /**
     * @brief builds the rlzsa
     * @tparam sad_t type of the values in SA^d
     * @tparam irr_pos_t position type (pos_t) for the index of rev(R)
     */
    template <typename sad_t, typename irr_pos_t>
    void build()
    {
        if (from_file) {
            for (uint16_t i = 0; i < p; i++) {
                SA_file_bufs.emplace_back(sdsl::int_vector_buffer<40>(
                    sa_file_path, std::ios::in, 128 * 1024, 40, true));
            }
        }

        bool _space = from_file || mode == _suffix_array_space;
        build_freq_sad<sad_t>();
        build_r<sad_t>();
        if (_space) store_r();
        build_idx_rev_r<sad_t, irr_pos_t>();

        if (_space) {
            load_r();
            build_rlzsa_factorization<true, sad_t, irr_pos_t>();
        } else {
            build_rlzsa_factorization<false, sad_t, irr_pos_t>();
        }
    }

    // ############################# rlzsa CONSTRUCTION METHOD DECLARATIONS #############################

    /**
     * @brief builds freq_SAd
     * @tparam sad_t type of the values in SA^d
     */
    template <typename sad_t>
    void build_freq_sad();

    /**
     * @brief builds R and rev(R)
     * @tparam sad_t type of the values in SA^d
     */
    template <typename sad_t>
    void build_r();

    /**
     * @brief builds the move-r index of rev(R)
     * @tparam sad_t type of the values in SA^d
     * @tparam irr_pos_t position type (pos_t) for the index of rev(R)
     */
    template <typename sad_t, typename irr_pos_t>
    void build_idx_rev_r();

    /**
     * @brief builds the rlzsa factorization
     * @tparam space true <=> store rlzsa data structures in files during the construction
     * @tparam sad_t type of the values in SA^d
     * @tparam irr_pos_t position type (pos_t) for the index of rev(R)
     */
    template <bool space, typename sad_t, typename irr_pos_t>
    void build_rlzsa_factorization();

    /**
     * @brief stores R to disk
     */
    void store_r();

    /**
     * @brief loads R from disk
     */
    void load_r();
};
