#pragma once

#include <cmath>
#include <ips4o.hpp>
#include <move_r/move_r.hpp>

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::preprocess_t(bool in_memory, bool bigbwt, std::ifstream* T_ifile)
{
    if (log)
        std::cout << "preprocessing T" << std::flush;

    if constexpr (byte_alphabet) {
        // contains_uchar_thr[i_p][c] = true <=> thread i_p found the character c in its section of T[0..n-2].
        std::vector<std::vector<uint8_t>> contains_uchar_thr(p, std::vector<uint8_t>(256, 0));

        if (in_memory) {
            // Iterate over T[0..n-2] and report the occurrence of each found character in T[0..n-2] in contains_uchar_thr.
            #pragma omp parallel for num_threads(p)
            for (uint64_t i = 0; i < n - 1; i++) {
                contains_uchar_thr[omp_get_thread_num()][T<i_sym_t>(i)] = 1;
            }
        } else {
            T_ifile->seekg(0, std::ios::end);
            n = T_ifile->tellg() + (std::streamsize) + 1;
            idx.n = n;
            T_ifile->seekg(0, std::ios::beg);

            pos_t max_t_buf_size = std::max((pos_t)1, n / 500);
            std::string T_buf;
            no_init_resize(T_buf, max_t_buf_size);
            pos_t cur_t_buf_size;
            pos_t n_ = n - 1;

            while (n_ > 0) {
                cur_t_buf_size = std::min(n_, max_t_buf_size);
                read_from_file(*T_ifile, T_buf.c_str(), cur_t_buf_size);

                // Iterate over T[0..n-2] and report the occurrence of each found character in T[0..n-2] in contains_uchar_thr.
                #pragma omp parallel for num_threads(p)
                for (uint64_t i = 0; i < cur_t_buf_size; i++) {
                    contains_uchar_thr[omp_get_thread_num()][char_to_uchar(T_buf[i])] = 1;
                }

                n_ -= cur_t_buf_size;
            }

            T_ifile->seekg(0, std::ios::beg);
        }

        // contains_uchar[c] = 1 <=> c in T[0..n-2].
        std::vector<uint8_t> contains_uchar(256, 0);

        /* Combine the results of each thread's sections in contains_uchar_thr[0..i_p-1][0..255]
        into contains_uchar[0..255]. */
        for (uint16_t cur_uchar = 0; cur_uchar < 256; cur_uchar++) {
            for (uint16_t i_p = 0; i_p < p; i_p++) {
                if (contains_uchar_thr[i_p][cur_uchar] == 1) {
                    contains_uchar[cur_uchar] = 1;
                    break;
                }
            }
        }

        contains_uchar_thr.clear();
        contains_uchar_thr.shrink_to_fit();

        // The number of distinct characters in T[0..n-1].
        idx.sigma = 1;

        // Count the number of ones in contains_uchar[0..255] in idx._map_ext.
        for (uint16_t cur_uchar = 0; cur_uchar < 256; cur_uchar++) {
            if (contains_uchar[cur_uchar] == 1) {
                idx.sigma++;
            }
        }

        bool contains_invalid_char = false;

        for (uint8_t cur_uchar = 0; cur_uchar < min_valid_char; cur_uchar++) {
            if (contains_uchar[cur_uchar] == 1) {
                contains_invalid_char = true;
                break;
            }
        }

        // If the input contains too many distinct characters, we have to remap the characters in T[0..n-2].
        if (bigbwt || contains_invalid_char) {
            if (idx.sigma > 256 - min_valid_char) {
                /* If T[0..n-2] contains more than 256 - min_valid_char distinct characters, we cannot remap them into the
                range [0..255] without using a character less than min_valid_char, hence we cannot build an index for T. */
                std::cout << "Error: the input contains more than " << std::to_string(256 - min_valid_char) << " distinct characters" << std::endl;
                return;
            }

            idx.symbols_remapped = true;

            // build the mapping function map_symbol that remaps the characters of T, s.t. it does
            // not contain 0 or 1; also build its inverse function unmap_symbol

            idx._map_int.resize(256, 0);
            idx._map_ext.resize(256, 0);

            /* To preserve the order among characters in T[0..n-2], we start by mapping smallest
            character in T[0..n-2] to 2, the second smallest to 3, ... . */

            // The character, to map the currently next largest character in T[0..n-2] to.
            uint16_t next_uchar_to_remap_to = min_valid_char;
            max_remapped_uchar = 0;

            for (uint16_t cur_uchar = 0; cur_uchar < (bigbwt ? 256 : next_uchar_to_remap_to); cur_uchar++) {
                if (contains_uchar[cur_uchar] == 1) {
                    idx._map_int[cur_uchar] = next_uchar_to_remap_to;
                    idx._map_ext[next_uchar_to_remap_to] = cur_uchar;
                    max_remapped_uchar = cur_uchar;
                    next_uchar_to_remap_to++;

                    if (bigbwt && next_uchar_to_remap_to == 10) {
                        // make sure the input does not contain '\n' = 10
                        next_uchar_to_remap_to++;
                    }
                }
            }

            max_remapped_to_uchar = next_uchar_to_remap_to - 1;

            if (!bigbwt) {
                for (uint16_t cur_uchar = max_remapped_to_uchar + 1; cur_uchar < 256; cur_uchar++) {
                    if (contains_uchar[cur_uchar] == 1) {
                        idx._map_int[cur_uchar] = cur_uchar;
                        idx._map_ext[cur_uchar] = cur_uchar;
                    }
                }
            }

            // Apply map_symbol to T.
            if (in_memory) {
                #pragma omp parallel for num_threads(p)
                for (uint64_t i = 0; i < n - 1; i++) {
                    if (T<i_sym_t>(i) <= max_remapped_uchar) {
                        T<i_sym_t>(i) = idx._map_int[T<i_sym_t>(i)];
                    }
                }
            }
        }

        if (!in_memory) {
            std::ofstream T_ofile(prefix_tmp_files);
            pos_t max_t_buf_size = std::max((pos_t)1, n / 500);
            std::vector<uint8_t> T_buf;
            no_init_resize(T_buf, max_t_buf_size);
            pos_t cur_t_buf_size;
            pos_t n_ = n - 1;

            while (n_ > 0) {
                cur_t_buf_size = std::min(n_, max_t_buf_size);
                read_from_file(*T_ifile, (char*)&T_buf[0], cur_t_buf_size);

                if (idx.symbols_remapped) {
                    #pragma omp parallel for num_threads(p)
                    for (uint64_t i = 0; i < cur_t_buf_size; i++) {
                        if (T_buf[i] <= max_remapped_uchar) {
                            T_buf[i] = idx._map_int[T_buf[i]];
                        }
                    }
                }

                write_to_file(T_ofile, (char*)&T_buf[0], cur_t_buf_size);
                n_ -= cur_t_buf_size;
            }

            T_ifile->close();
            T_ofile.close();
        }

        p_ = p;
    } else if (idx.sigma == 0) {
        idx.symbols_remapped = true;
        uint64_t alloc_before = malloc_count_current();

        for (pos_t i = 0; i < n - 1; i++) {
            idx._map_int.try_emplace(T<sym_t>(i), 0);
        }

        idx.sigma = idx._map_int.size() + 1;
        idx.size_map_int = malloc_count_current() - alloc_before;
        p_ = 1;
        no_init_resize(idx._map_ext, idx.sigma);
        idx._map_ext[0] = 0;
        pos_t sym_cur = 1;

        for (auto p : idx._map_int) {
            idx._map_ext[sym_cur] = p.first;
            sym_cur++;
        }

        if (p == 1) {
            ips4o::sort(idx._map_ext.begin() + 1, idx._map_ext.end());
        } else {
            ips4o::parallel::sort(idx._map_ext.begin() + 1, idx._map_ext.end());
        }

        if (idx.sigma > p) {
            #pragma omp parallel num_threads(p)
            {
                uint16_t i_p = omp_get_thread_num();

                i_sym_t b = 1 + i_p * ((idx.sigma - 1) / p);
                i_sym_t e = i_p == p - 1 ? idx.sigma - 1 : ((i_p + 1) * ((idx.sigma - 1) / p));

                for (i_sym_t i = b; i <= e; i++) {
                    idx._map_int[idx._map_ext[i]] = i;
                }
            }
        } else {
            for (i_sym_t i = 1; i < idx.sigma; i++) {
                idx._map_int[idx._map_ext[i]] = i;
            }
        }
    }

    if (log) {
        if (mf_idx != NULL)
            *mf_idx << " time_preprocess_t=" << time_diff_ns(time, now());
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
template <rlbwt_build_mode mode, typename sa_sint_t>
void move_r<support, sym_t, pos_t>::construction::build_rlbwt_c()
{
    if (log) {
        time = now();
        std::cout << "building RLBWT" << std::flush;
    }

    std::vector<sa_sint_t>& SA = get_sa<sa_sint_t>(); // [0..n-1] The suffix array

    r_p.resize(p_ + 1, 0);
    uint8_t width_bwt = byte_alphabet ? 1 : (uint8_t)std::ceil(std::log2(idx.sigma + 1) / (double)8);
    RLBWT.resize(p_, interleaved_vectors<uint32_t, uint32_t>({ width_bwt, 4 }));

    if constexpr (byte_alphabet) {
        C.resize(p_, std::vector<pos_t>(byte_alphabet ? 256 : idx.sigma, 0));
    }

    if constexpr (mode == _bwt_file) {
        for (uint16_t i = 0; i < p; i++) {
            BWT_file_bufs.emplace_back(sdsl::int_vector_buffer<8>(
                prefix_tmp_files + ".bwt", std::ios::in,
                128 * 1024, 8, true));
        }
    }

    for (uint16_t i = 0; i < p_; i++) {
        n_p.emplace_back(i * (n / p_));
    }

    n_p.emplace_back(n);

    #pragma omp parallel num_threads(p_)
    {
        // Index in [0..p_-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // Iteration range start position of thread i_p.
        pos_t b = n_p[i_p];
        // Iteration range end position of thread i_p.
        pos_t e = n_p[i_p + 1];

        // Store in C[i_p][c] the number of occurrences of c in L[b..e], for each c in [0..255].

        i_sym_t prev_sym; // previous symbol in L.
        i_sym_t cur_sym; // current symbol in L.
        pos_t i_ = b; // start position of the last seen run in L.

        if constexpr (mode == _sa) {
            prev_sym = SA[b] == 0 ? 0 : T<i_sym_t>(SA[b] - 1);
        } else if constexpr (mode == _bwt) {
            prev_sym = char_to_uchar(L[b]);
        } else if constexpr (mode == _bwt_file) {
            prev_sym = BWT_file_bufs[i_p][b];
            if (prev_sym == 2)
                prev_sym = 0;
        }

        // Iterate over the range L[b+1..e-1]
        for (pos_t i = b + 1; i < e; i++) {
            if constexpr (mode == _sa) {
                cur_sym = SA[i] == 0 ? 0 : T<i_sym_t>(SA[i] - 1);
            } else if constexpr (mode == _bwt) {
                cur_sym = char_to_uchar(L[i]);
            } else if constexpr (mode == _bwt_file) {
                cur_sym = BWT_file_bufs[i_p][i];
                if (cur_sym == 2)
                    cur_sym = 0;
            }

            // check if there is a run starting at L[i]
            if (cur_sym != prev_sym) {
                add_run(i_p, prev_sym, i - i_);
                if constexpr (byte_alphabet)
                    C[i_p][prev_sym] += i - i_;
                prev_sym = cur_sym;
                i_ = i;
            }
        }

        // add the run L[i'..e)
        add_run(i_p, prev_sym, e - i_);
        if constexpr (byte_alphabet)
            C[i_p][prev_sym] += e - i_;
        // Store in r_p[i_p] the number of runs starting in L[b..e).
        r_p[i_p] = RLBWT[i_p].size();
        RLBWT[i_p].shrink_to_fit();
    }

    // for i_p \in [1,p'-2], merge the last run in thread i_p's section with the first run in thread
    // i_p+1's section, if their characters are equal
    for (uint16_t i_p = 0; i_p < p_ - 1; i_p++) {
        i_sym_t c = run_sym(i_p, r_p[i_p] - 1);

        if (run_sym(i_p + 1, 0) == c) {
            pos_t l = run_len(i_p, r_p[i_p] - 1);
            set_run_len(i_p + 1, 0, run_len(i_p + 1, 0) + l);

            if constexpr (byte_alphabet) {
                C[i_p][c] -= l;
                C[i_p + 1][c] += l;
            }

            n_p[i_p + 1] -= l;
            r_p[i_p]--;
            RLBWT[i_p].resize(r_p[i_p]);
        }
    }

    if constexpr (mode == _bwt_file) {
        for (uint16_t i = 0; i < p; i++) {
            BWT_file_bufs[i].close();
        }

        BWT_file_bufs.clear();
        BWT_file_bufs.shrink_to_fit();
        std::filesystem::remove(prefix_tmp_files + ".bwt");
    }

    if (&L == &L_tmp) {
        L.clear();
        L.shrink_to_fit();
    }

    if (delete_T) {
        T_str.clear();
        T_str.shrink_to_fit();
        T_vec.clear();
        T_vec.shrink_to_fit();
    }

    if (support == _count && mode != _bwt) {
        SA.clear();
        SA.shrink_to_fit();
    }

    /* Now, r_p[i_p] stores the number of runs starting in the iteration range L[b..e] of thread
    i_p in [0..p'-1], and r_p[p'] = 0. We want r_p[i_p] to store the number of runs starting before
    the iteration range start position b of thread i_p in [0..p'-1]. Also, we want r_p[p'] to store
    the number of all runs. This way, we can build I_LF[r_p[i_p]..r_p[i_p+1]-1] with thread i_p in [0..p'-1]
    using the C-array in C[p'] while using and updating the rank-function in C[i_p] on-the-fly. */

    for (uint16_t i = p_; i > 0; i--) {
        r_p[i] = r_p[i - 1];
    }

    r_p[0] = 0;

    for (uint16_t i = 2; i <= p_; i++) {
        r_p[i] += r_p[i - 1];
    }

    /* Now, r_p[i_p] stores the number of runs starting before the iteration range start position b of
    thread i_p in [0..p'-1] and r_p[p'] stores the number of all runs, so we are done with r_p[0..p'] */

    r = r_p[p_];
    idx.r = r;

    if constexpr (int_alphabet) {
        C.emplace_back(std::vector<pos_t>(idx.sigma, 0));

        for (pos_t i = 0; i < r; i++) {
            C[0][run_sym(0, i)] += run_len(0, i);
        }
    }

    process_c();

    if (log) {
        if (mf_idx != NULL)
            *mf_idx << " time_build_rlbwt=" << time_diff_ns(time, now());
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::process_c()
{
    /* Now, C[i_p][c] is the number of occurrences of c in L[b..e], where [b..e] is the range of the
    thread i_p in [0..p'-1]. Also, we have C[p'][0..255] = 0. */

    /* We want to have C[i_p][c] = rank(L,c,b-1), where b is the iteration range start position of
    thread i_p in [0..p'-1]. Also, we want C[p'][0..255] to be the C-array, that is C[p'][c] stores
    the number of occurrences of all smaller characters c' < c in L[0..n-1], for c in [0..255]. */

    pos_t max_symbol = byte_alphabet ? 256 : idx.sigma;

    for (pos_t i = 1; i < p_; i++) {
        #pragma omp parallel for num_threads(p)
        for (uint64_t j = 0; j < max_symbol; j++) {
            C[i][j] += C[i - 1][j];
        }
    }

    /* Now, we have C[i_p][c] = rank(L,c,e), for each c in [0..255] and i_p in [0..p'-1],
    where e is the iteration range end position of thread i_p. */

    C.insert(C.begin(), std::vector<pos_t>());
    no_init_resize(C[0], max_symbol);

    #pragma omp parallel for num_threads(p)
    for (uint64_t j = 0; j < max_symbol; j++) {
        C[0][j] = 0;
    }

    /* Now, we have C[i_p][c] = rank(L,c,b-1), for each c in [0..255] and i_p in [0..p'-1],
    where b is the iteration range start position of thread i_p, so we are done with C[0..p'-1][0..255].
    Also, we have C[p'][c] = rank(L,c,n-1), for c in [0..255]. */

    pos_t cp_im1 = C[p_][0]; // stores at the start of the i-th iteration C[p'][i-1]
    C[p_][0] = 0;
    pos_t cp_im1_tmp;

    for (pos_t i = 1; i < max_symbol; i++) {
        cp_im1_tmp = C[p_][i];
        C[p_][i] = cp_im1 + C[p_][i - 1];
        cp_im1 = cp_im1_tmp;
    }

    // Now we are done with C, since C[p'] is the C-array.
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::build_ilf()
{
    if (log) {
        time = now();
        std::cout << "building I_LF" << std::flush;
    }

    no_init_resize(I_LF, r);

    #pragma omp parallel num_threads(p_)
    {
        // Index in [0..p'-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // Iteration range start position of thread i_p.
        pos_t b_r = r_p[i_p];

        // Number of BWT runs in thread i_p's section.
        pos_t rp_diff = r_p[i_p + 1] - r_p[i_p];

        // i', Start position of the last-seen run.
        pos_t i_ = n_p[i_p];

        // Build I_LF
        for (pos_t i = 0; i < rp_diff; i++) {
            /* Write the pair (i',LF(i')) to the next position i in I_LF, where
            LF(i') = C[L[i']] + rank(L,L[i'],i'-1) = C[p'][L[i']] + C[i_p][L[i']]. */
            I_LF[b_r + i] = std::make_pair(i_, C[p_][run_sym(i_p, i)] + C[i_p][run_sym(i_p, i)]);

            /* Update the rank-function in C[i_p] to store C[i_p][c] = rank(L,c,i'-1),
            for each c in [0..255] */
            C[i_p][run_sym(i_p, i)] += run_len(i_p, i);

            // Update the position of the last-seen run.
            i_ += run_len(i_p, i);
        }
    }

    C.clear();
    C.shrink_to_fit();

    if (log) {
        if (mf_idx != NULL)
            *mf_idx << " time_build_ilf=" << time_diff_ns(time, now());
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::build_mlf()
{
    if (log) {
        if (mf_mds != NULL) {
            *mf_mds << "RESULT"
                    << " type=build_mlf"
                    << " text=" << name_text_file
                    << " num_threads=" << p
                    << " a=" << idx.a;
        }
        time = now();
        std::cout << std::endl
                  << "building M_LF" << std::flush;
    }
    
    uint8_t width_l_ = byte_alphabet ? 8 :
        (uint8_t)(std::ceil(std::log2(idx.sigma + 1) / (double)8) * 8);

    idx._M_LF = move_data_structure_l_<pos_t, i_sym_t>(
        std::move(I_LF), n, {
            .num_threads = p,
            .a = idx.a,
            .log = log,
            .mf = mf_mds,
        },
        width_l_);

    r_ = idx._M_LF.num_intervals();
    idx.r_ = r_;

    if (log) {
        if (mf_mds != NULL)
            *mf_mds << std::endl;
        if (mf_idx != NULL) {
            *mf_idx << " time_build_mlf=" << time_diff_ns(time, now())
                    << " r_=" << r_;
        }
        std::cout << std::endl;
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
template <bool bigbwt, typename sa_sint_t>
void move_r<support, sym_t, pos_t>::construction::build_iphim1_sa()
{
    if (log) {
        time = now();
        std::cout << "building I_Phi^{-1}" << std::flush;
    }

    if constexpr (bigbwt) {
        for (uint16_t i = 0; i < p; i++) {
            SA_file_bufs.emplace_back(sdsl::int_vector_buffer<40>(
                prefix_tmp_files + ".sa", std::ios::in,
                128 * 1024, 40, true));
        }
    }

    no_init_resize(I_Phi_m1, r);

    I_Phi_m1[0] = std::make_pair(
        SA<bigbwt, sa_sint_t>(0, n - 1),
        SA<bigbwt, sa_sint_t>(0, 0));

    #pragma omp parallel num_threads(p_)
    {
        uint16_t i_p = omp_get_thread_num();

        // Iteration range start position of thread i_p.
        pos_t b_r = r_p[i_p];

        // Number of BWT runs within thread i_p's section.
        pos_t rp_diff = r_p[i_p + 1] - r_p[i_p];

        // start position of the current BWT run.
        pos_t j = n_p[i_p];

        if (rp_diff > 0) {
            if (r_p[i_p] != 0) {
                I_Phi_m1[b_r] = std::make_pair(
                    SA<bigbwt, sa_sint_t>(i_p, j - 1),
                    SA<bigbwt, sa_sint_t>(i_p, j));
            }

            j += run_len(i_p, 0);

            for (pos_t i = 1; i < rp_diff; i++) {
                I_Phi_m1[b_r + i] = std::make_pair(
                    SA<bigbwt, sa_sint_t>(i_p, j - 1),
                    SA<bigbwt, sa_sint_t>(i_p, j));

                j += run_len(i_p, i);
            }
        }
    }

    if (!build_sa_and_l && support != _locate_rlzdsa) {
        std::vector<sa_sint_t>& SA = get_sa<sa_sint_t>(); // [0..n-1] The suffix array

        SA.clear();
        SA.shrink_to_fit();
    }

    if (log) {
        if (mf_idx != NULL)
            *mf_idx << " time_build_iphi=" << time_diff_ns(time, now());
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
template <bool build_sas_>
void move_r<support, sym_t, pos_t>::construction::build_l__sas()
{
    if (log) {
        time = now();
        std::cout << "building L'" << (std::string)(build_sas_ ? " and SA_s" : "") << std::flush;
    }

    if constexpr (build_sas_) {
        if constexpr (support == _locate_move) {
            no_init_resize(SA_s, r_);
        } else {
            idx._SA_s.resize_no_init(r_);
        }
    }

    // Simultaneously iterate over the input intervals of M_LF nad the bwt runs to build L'
    #pragma omp parallel num_threads(p_)
    {
        // Index in [0..p'-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // Bwt range start position of thread i_p.
        pos_t b = n_p[i_p];

        // Iteration range start position of thread i_p.
        pos_t b_r = r_p[i_p];

        // Number of runs in thread i_p's section.
        pos_t rp_diff = r_p[i_p + 1] - r_p[i_p];

        // Index of the current input interval in M_LF, initially the index of the input interval of M_LF containing b
        pos_t j = bin_search_max_leq<pos_t>(b, 0, r_ - 1, [this](pos_t x) { return idx._M_LF.p(x); });
        // Starting position of the next bwt run.
        pos_t l_ = b;

        for (pos_t i = 0; i < rp_diff; i++) {
            idx._M_LF.template set_L_(j, run_sym(i_p, i));
            if constexpr (build_sas_) {
                if constexpr (support == _locate_move) {
                    SA_s[j] = I_Phi_m1[b_r + i].second;
                } else {
                    idx._SA_s.template set<0, pos_t>(j, I_Phi_m1[b_r + i].second);
                }
            }

            j++;

            // update l_ to the next run start position
            l_ += run_len(i_p, i);

            // iterate over all input intervals in M_LF within the i-th bwt run in thread i_p's section that have been
            // created by the balancing algorithm
            while (idx._M_LF.p(j) < l_) {
                if constexpr (build_sas_) {
                    if constexpr (support == _locate_move) {
                        SA_s[j] = n;
                    } else {
                        idx._SA_s.template set<0, pos_t>(j, n);
                    }
                }

                idx._M_LF.template set_L_(j, run_sym(i_p, i));
                j++;
            }
        }
    }

    n_p.clear();
    n_p.shrink_to_fit();

    r_p.clear();
    r_p.shrink_to_fit();

    RLBWT.clear();
    RLBWT.shrink_to_fit();

    if (log) {
        if (mf_idx != NULL)
            *mf_idx << " time_build_l__sas=" << time_diff_ns(time, now());
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::sort_iphim1()
{
    if (log) {
        time = now();
        std::cout << "sorting I_Phi^{-1}" << std::flush;
    }

    // Sort I_Phi^{-1} by the starting positions of its input intervals.
    auto comp_I_Phi = [](std::pair<pos_t, pos_t> p1, std::pair<pos_t, pos_t> p2) { return p1.first < p2.first; };

    // Choose the correct sorting algorithm.
    if (p > 1) {
        ips4o::parallel::sort(I_Phi_m1.begin(), I_Phi_m1.end(), comp_I_Phi);
    } else {
        ips4o::sort(I_Phi_m1.begin(), I_Phi_m1.end(), comp_I_Phi);
    }

    if (log) {
        if (mf_idx != NULL)
            *mf_idx << " time_sort_iphi=" << time_diff_ns(time, now());
        if (mf_mds != NULL) {
            *mf_mds << "RESULT"
                    << " type=build_mphi"
                    << " text=" << name_text_file
                    << " num_threads=" << p
                    << " a=" << idx.a;
        }
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::build_mphim1()
{
    if (log) {
        time = now();
        std::cout << std::endl
                  << "building M_Phi^{-1}" << std::flush;
    }

    idx._M_Phi_m1 = move_data_structure<pos_t>(
        std::move(I_Phi_m1), n, {
            .num_threads = p,
            .a = idx.a,
            .log = log,
            .mf = mf_mds,
        },
        &pi_mphi);

    r__ = idx._M_Phi_m1.num_intervals();
    idx.r__ = r__;

    if (log) {
        if (mf_mds != NULL)
            *mf_mds << std::endl;
        if (mf_idx != NULL) {
            *mf_idx << " time_build_mphi=" << time_diff_ns(time, now());
            *mf_idx << " r__=" << r__;
        }
        std::cout << std::endl;
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::build_saphim1()
{
    time = now();
    if (log)
        std::cout << "building SA_Phi^{-1}" << std::flush;

    no_init_resize(pi_, r_);

    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i < r_; i++) {
        pi_[i] = i;
    }

    auto comp_pi_ = [this](pos_t i, pos_t j) { return SA_s[i] < SA_s[j]; };
    if (p > 1) {
        ips4o::parallel::sort(pi_.begin(), pi_.end(), comp_pi_);
    } else {
        ips4o::sort(pi_.begin(), pi_.end(), comp_pi_);
    }

    idx.omega_idx = idx._M_Phi_m1.width_idx();
    idx._SA_Phi_m1 = interleaved_vectors<pos_t, pos_t>({ (uint8_t)(idx.omega_idx / 8) });
    idx._SA_Phi_m1.resize_no_init(r_);

    /* Now we will divide the range [0..n-1] up into p non-overlapping sub-ranges [s[i_p]..s[i_p+1]-1],
    for each i_p in [0..p-1], with 0 = s[0] < s[1] < ... < s[p] = n, where
    s[i_p] = min {s' in [0,n-1], s.t. x[i_p] + u[i_p] - 2 >= i_p * lfloor (r+r'')/p rfloor, where
                    x[i_p] = min {x' in [0,r''-1], s.t. M_Phi^{-1}.q(x') >= s'} and
                    u[i_p] = min {u' in [0,r-1], s.t. SA_s[u'] >= s'}
    }.
    By doing so, we ensure that the number of the output intervals of M_Phi^{-1} starting in the range
    [s[i_p]..s[i_p+1]-1] plus the number of suffix array samples in SA_s lying in the range
    [s[i_p]..s[i_p+1]-1] is lfloor (r+r'')/p rfloor +- 1. This property is useful, because it
    ensures that if with each thread i_p, we simultaneously iterate over those, then each thread
    iterates over almost exactly the same number lfloor (r+r'')/p rfloor +- 1 of entries in M_Phi^{-1}
    and SA_s combined. This way, we can acheive good load-balancing. Because we do not have to access
    s[0..p] later, we will not store those values in an array. */

    /* [0..p], x[i_p] = min {x' in [0,r''-1], s.t. M_Phi^{-1}.q(x') >= s'} stores the number of output
    intervals in M_Phi^{-1} starting before s[i_p]. */
    std::vector<pos_t> x(p + 1);
    x[0] = 0;
    x[p] = r__;

    /* [0..p], u[i_p] = min {u' in [0,r-1], s.t. SA_s[u'] >= s'} stores the number of suffix array
    samples in SA_s that are smaller than s[i_p]. */
    std::vector<pos_t> u(p + 1);
    u[0] = 0;
    u[p] = r;

    // Compute s[1..p-1], x[1..p-1] and u[1..p-1].
    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // The optimal value i_p * lfloor (r+r'')/p rfloor for s[i_p].
        pos_t o = i_p * ((r + r__) / p);

        // Left interval limit of the binary search for s[i_p].
        pos_t l_s;
        // Left interval limit of the binary search for x[i_p].
        pos_t l_x;
        // Left interval limit of the binary search for u[i_p].
        pos_t l_u;
        // Candidate position in the binary search for s[i_p].
        pos_t m_s;
        // Candidate position in the binary search for x[i_p].
        pos_t m_x;
        // Candidate position in the binary search for u[i_p].
        pos_t m_u;
        // Right interval limit of the binary search for s[i_p].
        pos_t r_s;
        // Right interval limit of the binary search for x[i_p].
        pos_t r_x;
        // Right interval limit of the binary search for u[i_p].
        pos_t r_u;

        // Initialize the search range for s[i_p] to [0..n-1].
        l_s = 0;
        r_s = n - 1;

        // Perform a binary search over [0..n-1], to find s[i_p].
        while (true) {
            /* Set the Candidate position for s[i_p] to the position in the middle
            between l_s and r_s. */
            m_s = l_s + (r_s - l_s) / 2;

            // Find the minimum x' in [0,r''-1], s.t. M_Phi^{-1}.q(x') >= m_s.
            l_x = 0;
            r_x = r__ - 1;
            while (l_x != r_x) {
                m_x = l_x + (r_x - l_x) / 2;
                if (idx._M_Phi_m1.q(pi_mphi[m_x]) < m_s) {
                    l_x = m_x + 1;
                } else {
                    r_x = m_x;
                }
            }

            // Find the minimum u' in [0,r-1], s.t. SA_s[pi'[u']] >= m_s.
            l_u = 0;
            r_u = r - 1;
            while (l_u != r_u) {
                m_u = l_u + (r_u - l_u) / 2;
                if (SA_s[pi_[m_u]] < m_s) {
                    l_u = m_u + 1;
                } else {
                    r_u = m_u;
                }
            }

            /* If l_s = r_s, then l_s is an optimal value for s[i_p] and l_x and
            l_u are valid values for x[i_p] and u[i_p], respectively. */
            if (l_s == r_s) {
                break;
            }

            // Else, adjust the range for the binary search over [0..n-1].
            if (l_x + l_u < o) {
                l_s = m_s + 1;
            } else {
                r_s = m_s;
            }
        }

        // Store l_x and l_u in x[i_p] and u[i_p], respectively.
        x[i_p] = l_x;
        u[i_p] = l_u;

        #pragma omp barrier

        // Iteration range start position in the output intervals of M_Phi^{-1}.
        pos_t i = x[i_p];
        // Iteration range start position in SA_s.
        pos_t j = u[i_p];
        // Iteration range end position + 1 in SA_s.
        pos_t j_ = u[i_p + 1];

        // Check if the range, over which the thread i_p has to iterate in SA_s, is empty
        if (j < j_) {
            while (idx._M_Phi_m1.q(pi_mphi[i]) != SA_s[pi_[j]]) {
                i++;
            }

            /* Iterate over SA_s[pi'[j]],SA_s[pi'[j+1]],...,SA_s[pi'[j'-1]] */
            while (j < j_) {
                // Skip the output intervals the balancing algorithm has added to I_Phi^{-1}
                while (idx._M_Phi_m1.q(pi_mphi[i]) != SA_s[pi_[j]]) {
                    i++;
                }

                idx.set_SA_Phi_m1(pi_[j], pi_mphi[i]);

                i++;
                j++;
            }
        }
    }

    /* Since we set SA_s[j] = n for each j-th input interval of M_LF, whiches starting position is not the starting position of a bwt run,
    * where j \in [0,r'), SA_s[pi'[r]],SA_s[pi'[r+1]],...,SA_s[pi'[r'-1]] = n holds, hence we set SA_Phi^{-1}[pi'[i]] = r'' for
    * i \in [r,r') to mark that we cannot recover SA_s[M_LF[x]] = M_Phi^{-1}.q[SA_Phi^{-1}[x]] for each x-th input interval of
    * M_LF, whiches starting position is not the starting position of a bwt run and where x \in [0,r'). */
    #pragma omp parallel for num_threads(p)
    for (uint64_t i = r; i < r_; i++) {
        idx.set_SA_Phi_m1(pi_[i], r__);
    }

    x.clear();
    x.shrink_to_fit();

    u.clear();
    u.shrink_to_fit();

    pi_mphi.clear();
    pi_mphi.shrink_to_fit();

    if (log) {
        if (mf_idx != NULL)
            *mf_idx << " time_build_saphim1=" << time_diff_ns(time, now());
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::build_de()
{
    if constexpr (supports_locate) {
        idx.p_r = std::min<pos_t>(256, std::max<pos_t>(1, r / 100));
    }

    idx._D_e.resize(idx.p_r - 1);

    #pragma omp parallel for num_threads(p)
    for (uint16_t i = 0; i < idx.p_r - 1; i++) {
        pos_t x = bin_search_min_geq<pos_t>((i + 1) * ((n - 1) / idx.p_r), 0, r - 1, [this](pos_t x) { return (((int64_t)SA_s[pi_[x]]) - 1) % n; });
        idx._D_e[i] = std::make_pair(pi_[x], (pos_t)((((int64_t)SA_s[pi_[x]]) - 1) % n));
    }

    SA_s.clear();
    SA_s.shrink_to_fit();

    pi_.clear();
    pi_.shrink_to_fit();
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::build_rsl_()
{
    if (log) {
        time = now();
        std::cout << "building RS_L'" << std::flush;
    }

    if constexpr (byte_alphabet) {
        idx._RS_L_ = rsl_t([this](pos_t i) { return idx.L_(i); }, 0, r_ - 1);
    } else {
        idx._RS_L_ = rsl_t([this](pos_t i) { return idx.L_(i); }, idx.sigma, 0, r_ - 1);
    }

    if (log) {
        if (mf_idx != NULL)
            *mf_idx << " time_build_rsl_=" << time_diff_ns(time, now());
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::store_rlbwt()
{
    if (log) {
        time = now();
        std::cout << "storing RLBWT to disk" << std::flush;
    }

    std::ofstream file_rlbwt(prefix_tmp_files + ".rlbwt");

    for (uint16_t i = 0; i < p_; i++) {
        RLBWT[i].serialize(file_rlbwt);
    }

    RLBWT.clear();
    RLBWT.shrink_to_fit();
    file_rlbwt.close();

    if (log) {
        if (mf_idx != NULL)
            *mf_idx << " time_store_rlbwt=" << time_diff_ns(time, now());
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::load_rlbwt()
{
    if (log) {
        time = now();
        std::cout << "loading RLBWT from disk" << std::flush;
    }

    std::ifstream file_rlbwt(prefix_tmp_files + ".rlbwt");
    RLBWT.resize(p_);

    for (uint16_t i = 0; i < p_; i++) {
        RLBWT[i].load(file_rlbwt);
    }

    file_rlbwt.close();
    std::filesystem::remove(prefix_tmp_files + ".rlbwt");

    if (log) {
        if (mf_idx != NULL)
            *mf_idx << " time_load_rlbwt=" << time_diff_ns(time, now());
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::store_mlf()
{
    if (log) {
        time = now();
        std::cout << "storing M_LF and L' to disk" << std::flush;
    }

    std::ofstream file_mlf(prefix_tmp_files + ".mlf");
    idx._M_LF.serialize(file_mlf);
    idx._M_LF = move_data_structure_l_<pos_t, i_sym_t>();
    file_mlf.close();

    if (log) {
        if (mf_idx != NULL)
            *mf_idx << " time_store_mlf=" << time_diff_ns(time, now());
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::load_mlf()
{
    if (log) {
        time = now();
        std::cout << "loading M_LF and L' from disk" << std::flush;
    }

    std::ifstream file_mlf(prefix_tmp_files + ".mlf");
    idx._M_LF.load(file_mlf);
    file_mlf.close();
    std::filesystem::remove(prefix_tmp_files + ".mlf");

    if (log) {
        if (mf_idx != NULL)
            *mf_idx << " time_load_mlf=" << time_diff_ns(time, now());
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::store_sas()
{
    if (log) {
        time = now();
        std::cout << "storing SA_s to disk" << std::flush;
    }

    std::ofstream file_sas(prefix_tmp_files + ".sas");
    write_to_file(file_sas, (char*)&SA_s[0], r_ * sizeof(pos_t));
    SA_s.clear();
    SA_s.shrink_to_fit();
    file_sas.close();

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::load_sas()
{
    if (log) {
        time = now();
        std::cout << "loading SA_s from disk" << std::flush;
    }

    std::ifstream file_sas(prefix_tmp_files + ".sas");
    no_init_resize(SA_s, r_);
    read_from_file(file_sas, (char*)&SA_s[0], r_ * sizeof(pos_t));
    file_sas.close();
    std::filesystem::remove(prefix_tmp_files + ".sas");

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::store_sas_idx()
{
    if (log) {
        time = now();
        std::cout << "storing SA_s to disk" << std::flush;
    }

    std::ofstream file_sas(prefix_tmp_files + ".sas");
    idx._SA_s.serialize(file_sas);
    idx._SA_s.clear();
    idx._SA_s.shrink_to_fit();
    file_sas.close();

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::load_sas_idx()
{
    if (log) {
        time = now();
        std::cout << "loading SA_s from disk" << std::flush;
    }

    std::ifstream file_sas(prefix_tmp_files + ".sas");
    idx._SA_s.load(file_sas);
    file_sas.close();
    std::filesystem::remove(prefix_tmp_files + ".sas");

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::store_rsl_()
{
    if (log) {
        time = now();
        std::cout << "storing RS_L' to disk" << std::flush;
    }

    std::ofstream file_rsl_(prefix_tmp_files + ".rsl_");
    idx._RS_L_.serialize(file_rsl_);
    file_rsl_.close();
    idx._RS_L_ = rsl_t();

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::load_rsl_()
{
    if (log) {
        time = now();
        std::cout << "loading RS_L' from disk" << std::flush;
    }

    std::ifstream file_rsl_(prefix_tmp_files + ".rsl_");
    idx._RS_L_.load(file_rsl_);
    file_rsl_.close();
    std::filesystem::remove(prefix_tmp_files + ".rsl_");

    if (log) {
        time = log_runtime(time);
    }
}