#pragma once

#include <vector>
#include <iostream>
#include <filesystem>
#include <omp.h>
#include <sdsl/wavelet_trees.hpp>

#include <move_r/misc/utils.hpp>
#include <move_r/data_structures/hybrid_bit_vector.hpp>
#include <move_r/data_structures/interleaved_vectors.hpp>

/**
 * @brief a rank-select data structure using hybrid bit vectors (either sd_array or plain bit vector), or the rs data structure
 * @tparam sym_t symbol type
 * @tparam pos_t position type
 * @tparam build_rank_support
 * @tparam build_select_support
 */
template <typename sym_t, typename pos_t = uint32_t, bool build_rank_support = true, bool build_select_support = true>
class rank_select_support {
    protected:

    static_assert(std::is_same_v<pos_t,uint32_t> || std::is_same_v<pos_t,uint64_t>);
    static_assert(build_rank_support || build_select_support);
    
    static_assert(
        std::is_same_v<sym_t,char> ||
        std::is_same_v<sym_t,uint8_t> ||
        std::is_same_v<sym_t,uint16_t> ||
        std::is_same_v<sym_t,uint32_t> ||
        std::is_same_v<sym_t,uint64_t>
    );

    static constexpr bool str_input = std::is_same_v<sym_t,char>; // true <=> the input is a string
    static constexpr bool byte_alphabet = sizeof(sym_t) == 1; // true <=> the input uses a byte alphabet
    static constexpr bool int_alphabet = !byte_alphabet; // true <=> the input uses an integer alphabet

    using i_sym_t = std::conditional_t<str_input,uint8_t,sym_t>; // internal (unsigned) symbol type
    using inp_t = std::conditional_t<str_input,std::string,std::vector<sym_t>>; // input container type
    using hybrid_bv_t = hybrid_bit_vector<pos_t,build_rank_support,false,build_select_support>; // hybrid bit vector type

    /* maximum number of occurrences to use scanning instead of binary search over the occurrences for answering rank */
    static constexpr pos_t max_occ_scan_rank = 16;
    /* minimum number of occurrences to use a bit vector for answering rank */
    static constexpr pos_t min_occ_vec_rank = 512;
    /* maximum relative number of occurrences to use an sd_array instead of a plain bitvector for answering rank (and select) */
    static constexpr double thrsh_sd_array = hybrid_bit_vector<pos_t>::compression_threshold;
    /* maximum relative number of occurrences to store all occurrences plainly to answer select */
    static constexpr double thrsh_plain_select = 0.01;

#ifndef BENCH_RANK_SELECT
    /* minimum number of occurrences to build a bit vector that marks all occurrences */
    static constexpr pos_t min_occ_vec = min_occ_vec_rank;
#else
    /* minimum number of occurrences to build a bit vector that marks all occurrences */
    static constexpr pos_t min_occ_vec = 3;
    /* maximum number of occurrences to answer select with a lookup */
    pos_t max_occ_lookup_select = 0;
#endif

    pos_t input_size = 0; // the size of the input
    pos_t sigma = 0; // the number of distinct symbols in the input
    pos_t num_vectors = 0; // the number of initialized vectors in hyb_bit_vecs

    // ############################# DATA STRUCTURES #############################

    /**
     * @brief hyb_bit_vecs[i] contains a hybrid bit vector [0..size-1] that marks (with ones) the occurrences of v in the
     *        input, for each value v occurring in the input(, where i = vec_idx[v], for int_alphabet = true and i = v, else)
     */
    std::vector<hybrid_bv_t> hyb_bit_vecs;

    // ############################# DATA STRUCTURES FOR byte_alphabet = true #############################

    /**
     * @brief [0..sigma-1] contains at position sym the frequency of sym in the input
     */
    std::vector<pos_t> freq;
    
    // ############################# DATA STRUCTURES FOR int_alphabet = true #############################
    
    /**
     * @brief [0..sigma-1] vec_idx[v] stores the position in hyb_bit_vecs of the the bit vector marking
     *        the occurrences of v in the input, if freq(v) > min_occ_vec; else vec_idx[v] = sigma
     */
    interleaved_vectors<pos_t,pos_t> vec_idx;

    /**
     * @brief [0..sigma-1] stores at position v the sum of the frequencies of all v' < v with freq(v') > min_occ_vec
     */
    interleaved_vectors<pos_t,pos_t> c_arr;
    
    /**
     * @brief stores at position i the position of the j-th occurrence of v
     *        in the input, where c_arr[v] <= i < c_arr[v+1] and j = i - c_arr[v] + 1
     */
    interleaved_vectors<pos_t,pos_t> occs;
    
    // ##########################################################

    pos_t symbol_idx(sym_t sym) const {
        if constexpr (str_input) {
            return char_to_uchar(sym);
        } else {
            return sym;
        }
    };

    /**
     * @brief builds the bit vectors
     * @param read function to read the input with; it is called with i in [l,r]
     * as a parameter and must return the value of the input at index i
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     */
    void build(const std::function<sym_t(pos_t)>& read, pos_t l, pos_t r) {
        input_size = r-l+1;
        uint8_t bytes_per_entry = 0;
        pos_t alphabet_range = byte_alphabet ? 256 : sigma;
        freq.resize(alphabet_range,0);

        for (pos_t i=l; i<=r; i++) {
            freq[symbol_idx(read(i))]++;
        }

        if constexpr (int_alphabet) {
            bytes_per_entry = (uint8_t)std::ceil(std::log2(input_size+1)/(double)8);
            c_arr = interleaved_vectors<pos_t,pos_t>({bytes_per_entry});
            c_arr.resize_no_init(alphabet_range+1);
            c_arr.template set<0,pos_t>(0,0);
        }

#ifndef BENCH_RANK_SELECT
        pos_t max_occ_plain = build_select_support ? input_size * thrsh_plain_select : min_occ_vec;
#else
        pos_t max_occ_plain = std::numeric_limits<pos_t>::max();
        max_occ_lookup_select = input_size * thrsh_plain_select;
#endif

        if constexpr (byte_alphabet) {
            for (pos_t v=0; v<256; v++) {
                if (freq[v] != 0) {
                    sigma++;
                }
            }
        } else {
            for (pos_t v=0; v<alphabet_range; v++) {
                if (freq[v] <= max_occ_plain) {
                    c_arr.template set<0,pos_t>(v+1,c_arr[v]+freq[v]);
                } else {
                    c_arr.template set<0,pos_t>(v+1,c_arr[v]);
                }
            }
        }
        
        std::vector<pos_t> occ_idx;
        std::vector<sdsl::sd_vector_builder> sdv_builders;
        std::vector<sdsl::bit_vector> plain_bvs;
        pos_t max_occ_sd_array = input_size * thrsh_sd_array;
        num_vectors = byte_alphabet ? 256 : 0;
        pos_t min_occ_vec_tmp = 0;

        if constexpr (int_alphabet) {
            min_occ_vec_tmp = std::min<pos_t>(min_occ_vec,max_occ_plain);
            vec_idx = interleaved_vectors<pos_t,pos_t>({
                (uint8_t)std::ceil(std::log2(sigma+1)/(double)8)});
            vec_idx.resize_no_init(sigma);
            no_init_resize(occ_idx,alphabet_range);
        } else {
            sdv_builders.resize(256);
            plain_bvs.resize(256);
        }

        for (pos_t v=0; v<alphabet_range; v++) {
            if constexpr (int_alphabet) {
                occ_idx[v] = c_arr[v];

                if (freq[v] > min_occ_vec_tmp) {
                    if (freq[v] <= max_occ_sd_array) {
                        sdv_builders.emplace_back(sdsl::sd_vector_builder(input_size,freq[v]));
                        plain_bvs.emplace_back(sdsl::bit_vector());
                    } else {
                        sdv_builders.emplace_back(sdsl::sd_vector_builder());
                        plain_bvs.emplace_back(sdsl::bit_vector(input_size));
                    }

                    vec_idx.template set<0,pos_t>(v,num_vectors);
                    num_vectors++;
                } else {
                    vec_idx.template set<0,pos_t>(v,sigma);
                }
            } else {
                if (freq[v] != 0) {
                    if (freq[v] <= max_occ_sd_array) {
                        sdv_builders[v] = sdsl::sd_vector_builder(input_size,freq[v]);
                    } else {
                        plain_bvs[v] = sdsl::bit_vector(input_size);
                    }
                }
            }
        }

        if constexpr (int_alphabet) {
            occs = interleaved_vectors<pos_t,pos_t>({bytes_per_entry});
            occs.resize_no_init(c_arr[alphabet_range]);
        }

        for (pos_t i=l; i<=r; i++) {
            pos_t v = symbol_idx(read(i));

            if constexpr (int_alphabet) {
                if (freq[v] <= max_occ_plain) {
                    occs.template set<0,pos_t>(occ_idx[v],i-l);
                    occ_idx[v]++;
                }

                if (freq[v] > min_occ_vec_tmp) {
                    if (freq[v] <= max_occ_sd_array) {
                        sdv_builders[vec_idx[v]].set(i-l);
                    } else {
                        plain_bvs[vec_idx[v]][i-l] = 1;
                    }
                }
            } else {
                if (freq[v] <= max_occ_sd_array) {
                    sdv_builders[v].set(i-l);
                } else {
                    plain_bvs[v][i-l] = 1;
                }
            }
        }

        if constexpr (int_alphabet) {
            freq.clear();
            freq.shrink_to_fit();
        }

        occ_idx.clear();
        occ_idx.shrink_to_fit();
        hyb_bit_vecs.resize(num_vectors);

        for (pos_t i=0; i<num_vectors; i++) {
            if (plain_bvs[i].size() != 0) {
                hyb_bit_vecs[i] = hybrid_bv_t(std::move(plain_bvs[i]));
            } else if (sdv_builders[i].size() != 0) {
                hyb_bit_vecs[i] = hybrid_bv_t(sdsl::sd_vector<>(sdv_builders[i]));
            }
        }
    }

    public:
    rank_select_support() = default;

    /**
     * @brief builds the data structure from the range [l,r] in input (0 <= l <= r < input.size()),
     * else if l > r, then the the data structure is built for the whole input
     * @param input the input
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     */
    rank_select_support(const std::string& input,pos_t l = 1, pos_t r = 0) requires(byte_alphabet) {
        if (l > r) {
            l = 0;
            r = input.size()-1;
        }
        
        r = std::min<pos_t>(r,input.size()-1);
        build([&input](pos_t i){return input[i];},l,r);
    }
    
    /**
     * @brief builds the data structure by reading the input using the function read
     * @param read function to read the input with; it is called with i in [l,r]
     * as a parameter and must return the value of the input at index i
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     */
    rank_select_support(const std::function<sym_t(pos_t)>& read, pos_t l = 1, pos_t r = 0) requires(byte_alphabet) {
        build(read,l,r);
    }

    /**
     * @brief builds the data structure from the range [l,r] in input (0 <= l <= r < input.size()),
     * else if l > r, then the the data structure is built for the whole input
     * @param input the input
     * @param alphabet_size maximum value in the input + 1 (should also be the number of distinct
     * values in the input)
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     */
    rank_select_support(const std::vector<sym_t>& input, pos_t alphabet_size, pos_t l = 1, pos_t r = 0) requires(int_alphabet) {
        if (l > r) {
            l = 0;
            r = input.size()-1;
        }
        
        r = std::min<pos_t>(r,input.size()-1);
        sigma = alphabet_size;
        build([&input](pos_t i){return input[i];},l,r);
    }
    
    /**
     * @brief builds the data structure by reading the input using the function read
     * @param read function to read the input with; it is called with i in [l,r]
     * as a parameter and must return the value of the input at index i
     * @param alphabet_size maximum value in the input + 1 (should also be the number of distinct
     * values in the input)
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     */
    rank_select_support(const std::function<sym_t(pos_t)>& read, pos_t alphabet_size, pos_t l = 1, pos_t r = 0) requires(int_alphabet) {
        sigma = alphabet_size;
        build(read,l,r);
    }

    /**
     * @brief returns the size of the input
     * @return the size of the input
     */
    inline pos_t size() const {
        return input_size;
    }

    /**
     * @brief returns the number of distinct values in the input
     * @return the number of distinct values in the input
     */
    inline pos_t alphabet_size() const {
        return sigma;
    }

    /**
     * @brief returns whether the input vector is empty
     * @return whether the input vector is empty
     */
    inline bool empty() const {
        return size() == 0;
    }

    /**
     * @brief returns the size of the data structure in bytes
     * @return size of the data structure in bytes
     */
    uint64_t size_in_bytes() const {
        uint64_t size = 0;
        
        size += vec_idx.size_in_bytes();
        size += c_arr.size_in_bytes();
        size += occs.size_in_bytes();

        for (pos_t i=0; i<num_vectors; i++) {
            size += hyb_bit_vecs[i].size_in_bytes();
        }

        return size;
    }

    /**
     * @brief returns whether v occurs in the input
     * @param v value
     * @return whether v occurs in the input
     */
    inline bool contains(sym_t v) const {
        if constexpr (byte_alphabet) {
            return frequency(v) != 0;
        } else {
            return v < sigma;
        }
    }

    /**
     * @brief returns the number of occurrences of the symbol v in the input (v must occur in the input)
     * @param v symbol that occurs in the input
     * @return the number of occurrences of v in the input
     */
    inline pos_t frequency(sym_t v) const {
        if constexpr (byte_alphabet) {
            return freq[symbol_idx(v)];
        } else {
            pos_t diff = c_arr[v+1]-c_arr[v];

            if (diff != 0) {
                return diff;
            } else {
                return hyb_bit_vecs[vec_idx[v]].num_ones();
            }
        }
    }
    
#ifdef BENCH_RANK_SELECT
    /**
     * @brief returns the number of occurrences of v before index i
     * @param v [0..alphabet_size-1] value that occurs in the input
     * @param i [1..input size] position in the input
     * @return number of occurrences of v before index i
     */
    inline pos_t rank_vec(sym_t v, pos_t i) const requires(int_alphabet) {
        static_assert(build_rank_support);
        return hyb_bit_vecs[vec_idx[v]].rank_1(i);
    }

    /**
     * @brief returns the number of occurrences of v before index i
     * @param v [0..alphabet_size-1] value that occurs in the input
     * @param i [1..input size] position in the input
     * @return number of occurrences of v before index i
     */
    inline pos_t rank_scan(sym_t v, pos_t i) const requires(int_alphabet) {
        static_assert(build_rank_support);
        
        pos_t v_s = c_arr[v];
        pos_t v_e = c_arr[v+1];
        pos_t pos = v_s;

        if (occs[v_s] >= i) {
            return 0;
        }

        pos = v_s;
        v_e--;

        while (pos < v_e && occs[pos+1] < i) {
            pos++;
        }

        return pos-v_s+1;
    }

    /**
     * @brief returns the number of occurrences of v before index i
     * @param v [0..alphabet_size-1] value that occurs in the input
     * @param i [1..input size] position in the input
     * @return number of occurrences of v before index i
     */
    inline pos_t rank_bin_search(sym_t v, pos_t i) const requires(int_alphabet) {
        static_assert(build_rank_support);
        pos_t v_s = c_arr[v];
        pos_t v_e = c_arr[v+1];
        pos_t pos = bin_search_max_lt<pos_t>(i,v_s,v_e-1,[this](pos_t x){return occs[x];});
        return occs[pos] >= i ? 0 : pos-v_s+1;
    }
#endif

    /**
     * @brief returns the number of occurrences of v before index i
     * @param v [0..alphabet_size-1] value that occurs in the input
     * @param i [1..input size] position in the input
     * @return number of occurrences of v before index i
     */
    inline pos_t rank(sym_t v, pos_t i) const {
        static_assert(build_rank_support);

        if constexpr (byte_alphabet) {
            return hyb_bit_vecs[symbol_idx(v)].rank_1(i);
        } else {
            pos_t v_s = c_arr[v];
            pos_t v_e = c_arr[v+1];

            if (v_e == v_s || v_e-v_s > min_occ_vec_rank) {
                return hyb_bit_vecs[vec_idx[v]].rank_1(i);
            } else {
                pos_t pos;

                if (v_e-v_s <= max_occ_scan_rank) {
                    if (occs[v_s] >= i) {
                        return 0;
                    }

                    pos = v_s;
                    v_e--;

                    while (pos < v_e && occs[pos+1] < i) {
                        pos++;
                    }

                    return pos-v_s+1;
                } else {
                    pos = bin_search_max_lt<pos_t>(i,v_s,v_e-1,[this](pos_t x){return occs[x];});
                    return occs[pos] >= i ? 0 : pos-v_s+1;
                }
            }
        }
    }

#ifdef BENCH_RANK_SELECT
    /**
     * @brief returns the index of the i-th occurrence of v
     * @param v a value that occurs in the input
     * @param i [1..number of occurrences of v in the input]
     * @return index of the i-th occurrence of v
     */
    inline pos_t select_lookup(sym_t v, pos_t i) const requires(int_alphabet) {
        static_assert(build_select_support);
        return occs[c_arr[v]+i-1];
    }

    /**
     * @brief returns the index of the i-th occurrence of v
     * @param v a value that occurs in the input
     * @param i [1..number of occurrences of v in the input]
     * @return index of the i-th occurrence of v
     */
    inline pos_t select_vec(sym_t v, pos_t i) const requires(int_alphabet) {
        static_assert(build_select_support);
        return hyb_bit_vecs[vec_idx[v]].select_1(i);
    }
#endif

    /**
     * @brief returns the index of the i-th occurrence of v
     * @param v a value that occurs in the input
     * @param i [1..number of occurrences of v in the input]
     * @return index of the i-th occurrence of v
     */
    inline pos_t select(sym_t v, pos_t i) const {
        static_assert(build_select_support);

        if constexpr (byte_alphabet) {
            return hyb_bit_vecs[symbol_idx(v)].select_1(i);
        } else {
            pos_t v_s = c_arr[v];
            pos_t v_e = c_arr[v+1];
            
#ifndef BENCH_RANK_SELECT
            if (v_e != v_s) {
#else
            if (v_e-v_s <= max_occ_lookup_select) {
#endif
                return occs[v_s+i-1];
            } else {
                return hyb_bit_vecs[vec_idx[v]].select_1(i);
            }
        }
    }

    /**
     * @brief serializes the rank_select_support to an output stream
     * @param out output stream
     */
    void serialize(std::ostream& out) const {
        out.write((char*)&input_size,sizeof(pos_t));

        if (input_size != 0) {
            out.write((char*)&sigma,sizeof(pos_t));
            out.write((char*)&num_vectors,sizeof(pos_t));
            vec_idx.serialize(out);
            c_arr.serialize(out);
            occs.serialize(out);

            if constexpr (byte_alphabet) {
                out.write((char*)&freq[0],256*sizeof(pos_t));
            }

            for (pos_t i=0; i<num_vectors; i++) {
                hyb_bit_vecs[i].serialize(out);
            }
        }
    }

    /**
     * @brief loads the rank_select_support from an input stream
     * @param in input stream
     */
    void load(std::istream& in) {
        in.read((char*)&input_size,sizeof(pos_t));

        if (input_size != 0) {
            in.read((char*)&sigma,sizeof(pos_t));
            in.read((char*)&num_vectors,sizeof(pos_t));
            vec_idx.load(in);
            c_arr.load(in);
            occs.load(in);
            hyb_bit_vecs.resize(num_vectors);

            if constexpr (byte_alphabet) {
                no_init_resize(freq,256);
                in.read((char*)&freq[0],256*sizeof(pos_t));
            }

            for (pos_t i=0; i<num_vectors; i++) {
                hyb_bit_vecs[i].load(in);
            }
        }
    }

    std::ostream& operator>>(std::ostream& os) const {
        serialize(os);
        return os;
    }

    std::istream& operator<<(std::istream& is) {
        load(is);
        return is;
    }
};