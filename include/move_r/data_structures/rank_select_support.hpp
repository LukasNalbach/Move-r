#pragma once

#include <vector>
#include <iostream>
#include <filesystem>
#include <sdsl/wt_gmr.hpp>
#include <move_r/data_structures/hybrid_bit_vector.hpp>
#include <move_r/data_structures/interleaved_vectors.hpp>

/**
 * @brief a rank-select data structure using hybrid bit vectors (either sd_array or plain bit vector), or the gmr data structure
 * @tparam sym_t symbol type type
 * @tparam pos_t position type
 */
template <typename sym_t, typename pos_t = uint32_t>
class rank_select_support {
    protected:

    static_assert(std::is_same_v<pos_t,uint32_t> || std::is_same_v<pos_t,uint64_t>);
    
    static_assert(
        std::is_same_v<sym_t,char> ||
        std::is_same_v<sym_t,uint8_t> ||
        std::is_same_v<sym_t,uint16_t> ||
        std::is_same_v<sym_t,uint32_t> ||
        std::is_same_v<sym_t,uint64_t>
    );

    static constexpr bool str_input = std::is_same_v<sym_t,char>; // true <=> the input is a string
    static constexpr bool byte_alphabet = sizeof(sym_t) == 1; // true <=> use bit vectors
    static constexpr bool int_alphabet = !byte_alphabet; // true <=> use gmr

    using i_sym_t = std::conditional_t<str_input,uint8_t,sym_t>; // internal (unsigned) symbol type
    using hybrid_bv_t = hybrid_bit_vector<pos_t,true,false,false>; // hybrid bit vector type

    // ############################# COMMON VARIABLES #############################

    /**
     * @brief [0..sigma-1] stores at position c the sum of the frequencies of all c' < c that occur in the input
     */
    interleaved_vectors<pos_t,pos_t> c_array;
    
    /**
     * @brief [0..input_size-1] stores at position i the position of the j-th occurrence of c
     *        in the input, where c_array[c] <= i < c_array[c+1] and j = i - c_array[c] + 1
     */
    interleaved_vectors<pos_t,pos_t> occs;

    // ############################# VARIABLES FOR byte_alphabet = true #############################

    pos_t input_size = 0; // the size of the input
    std::vector<i_sym_t> alphabet; // the symbols occurring in the input sorted descendingly by their frequency
    std::vector<pos_t> freq; // stores at position c the frequency of c in the input

    /**
     * @brief hyb_bit_vecs[c] contains a hybrid bit vector [0..size-1] that marks (with ones) the occurrences of c
     *        in the input, for each value c occurring in the input
     */
    std::vector<hybrid_bv_t> hyb_bit_vecs;

    // ############################# VARIABLES FOR byte_alphabet = false #############################

    pos_t sigma = 0; // the number of distinct symbols in the input

    /**
     * @brief rank-select data structure for integer alphabets
     */
    sdsl::wt_gmr_rs<> gmr;

    // ######################################################################################

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
     * @param num_threads the number of threads to use
     */
    void build_bvs(const std::function<sym_t(pos_t)>& read, pos_t l, pos_t r, uint16_t num_threads) {
        input_size = r-l+1;
        std::vector<std::vector<pos_t>> occurrences_thr(num_threads,std::vector<pos_t>(256,0));

        #pragma omp parallel for num_threads(num_threads)
        for (uint64_t i=l; i<=r; i++) {
            occurrences_thr[omp_get_thread_num()][symbol_idx(read(i))]++;
        }

        freq.resize(256);

        for (uint16_t cur_char=0; cur_char<256; cur_char++) {
            for (uint16_t i_p=0; i_p<num_threads; i_p++) {
                freq[cur_char] += occurrences_thr[i_p][cur_char];
            }

            if (freq[cur_char] != 0) {
                alphabet.emplace_back(cur_char);
            }
        }

        occurrences_thr.clear();
        occurrences_thr.shrink_to_fit();
        sigma = alphabet.size();
        uint16_t num_chars = sigma;
        std::sort(alphabet.begin(),alphabet.end(),[this](i_sym_t a, i_sym_t b){return freq[a] > freq[b];});
        std::vector<uint8_t> is_compressed(256,1);
        std::vector<sdsl::sd_vector_builder> sdv_builders(256);
        std::vector<sdsl::bit_vector> plain_bvs(256);

        #pragma omp parallel for num_threads(num_threads)
        for (uint16_t i=0; i<num_chars; i++) {
            uint8_t cur_char = alphabet[i];

            if (freq[cur_char] > input_size * hybrid_bit_vector<pos_t>::compression_threshold) {
                is_compressed[cur_char] = 0;
                plain_bvs[cur_char] = sdsl::bit_vector(input_size);
            } else {
                sdv_builders[cur_char] = sdsl::sd_vector_builder(input_size,freq[cur_char]);
            }
        }

        uint8_t cur_char;

        for (pos_t i=l; i<=r; i++) {
            cur_char = symbol_idx(read(i));

            if (is_compressed[cur_char]) {
                sdv_builders[cur_char].set(i-l);
            } else {
                plain_bvs[cur_char][i-l] = 1;
            }
        }

        hyb_bit_vecs.resize(256);

        #pragma omp parallel for num_threads(num_threads)
        for (uint16_t i=0; i<num_chars; i++) {
            uint8_t cur_char = symbol_idx(alphabet[i]);

            if (is_compressed[cur_char]) {
                hyb_bit_vecs[cur_char] = hybrid_bv_t(sdsl::sd_vector<>(sdv_builders[cur_char]));
            } else {
                hyb_bit_vecs[cur_char] = hybrid_bv_t(std::move(plain_bvs[cur_char]));
            }
        }
    }

    /**
     * @brief builds the gmr data structure
     * @param read function to read the input with; it is called with i in [l,r]
     * as a parameter and must return the value of the input at index i
     * @param alphabet_size number of distinct symbols in vector (= maximum value in the input)
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     */
    void build_gmr(const std::function<sym_t(pos_t)>& read, pos_t alphabet_size, pos_t l, pos_t r) {
        sigma = alphabet_size;
        std::string vector_buffer_filename = "vec_" + random_alphanumeric_string(10);
        sdsl::int_vector_buffer<> vector_buffer(vector_buffer_filename,std::ios::out);
        freq.resize(alphabet_size,0);

        for (uint64_t i=l; i<=r; i++) {
            sym_t c = read(i);
            vector_buffer[i-l] = c;
            freq[c]++;
        }

        gmr = sdsl::wt_gmr_rs<>(vector_buffer,r-l+1);
        vector_buffer.close();
        std::filesystem::remove(vector_buffer_filename);
    }

    /**
     * @brief builds the c-array and occurrences
     * @param read function to read the input with; it is called with i in [l,r]
     * as a parameter and must return the value of the input at index i
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     */
    void build_select(const std::function<sym_t(pos_t)>& read, pos_t l, pos_t r, uint16_t num_threads) {
        uint8_t bytes_occs = (uint8_t)std::ceil(std::log2(r-l+1)/(double)8);
        uint8_t bytes_c_arr = byte_alphabet ? sizeof(pos_t) : bytes_occs;
        pos_t size_c_arr = (byte_alphabet ? 256 : sigma)+1;
        c_array = interleaved_vectors<pos_t,pos_t>({bytes_c_arr});
        c_array.resize_no_init(size_c_arr);
        c_array.template set<0,pos_t>(0,0);

        for (pos_t c=1; c<size_c_arr; c++) {
            c_array.template set<0,pos_t>(c,c_array[c-1]+freq[c-1]);
        }
        
        if constexpr (int_alphabet) {
            freq.clear();
            freq.shrink_to_fit();
        }

        std::vector<pos_t> occ_idx;
        no_init_resize(occ_idx,size_c_arr-1);

        #pragma omp parallel for num_threads(num_threads)
        for (uint64_t c=0; c<size_c_arr-1; c++) {
            occ_idx[c] = c_array[c];
        }

        occs = interleaved_vectors<pos_t,pos_t>({bytes_occs});
        occs.resize_no_init(r-l+1);

        for (pos_t i=l; i<=r; i++) {
            pos_t c = symbol_idx(read(i));
            occs.template set<0,pos_t>(occ_idx[c],i-l);
            occ_idx[c]++;
        }

    }
    public:
    rank_select_support() = default;

    /**
     * @brief builds the data structure from the range [l,r] in an input (0 <= l <= r < string.size()),
     * else if l > r, then the the data structure is built for the whole input
     * @param string the input
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     * @param num_threads the number of threads to use
     */
    rank_select_support(
        const std::string& string,
        pos_t l = 1, pos_t r = 0,
        uint16_t num_threads = omp_get_max_threads()
    ) requires(byte_alphabet) {
        if (l > r) {
            l = 0;
            r = string.size()-1;
        }
        
        r = std::min<pos_t>(r,string.size()-1);
        build_bvs([&string](uint i){return string[i];},l,r,num_threads);
        build_select([&string](uint i){return string[i];},l,r,num_threads);
    }
    
    /**
     * @brief builds the data structure by reading the input using the function read
     * @param read function to read the input with; it is called with i in [l,r]
     * as a parameter and must return the value of the input at index i
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     * @param num_threads the number of threads to use
     */
    rank_select_support(
        const std::function<sym_t(pos_t)>& read,
        pos_t l = 1, pos_t r = 0,
        uint16_t num_threads = omp_get_max_threads()
    ) requires(byte_alphabet) {
        build_bvs(read,l,r,num_threads);
        build_select(read,l,r,num_threads);
    }

    /**
     * @brief builds the data structure by reading the input vector
     * @param vector input vector
     * @param alphabet_size number of distinct symbols in vector (= maximum value in vector)
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     * @param num_threads the number of threads to use
     */
    rank_select_support(
        const std::vector<sym_t>& vector,
        pos_t alphabet_size,
        pos_t l = 1, pos_t r = 0,
        uint16_t num_threads = omp_get_max_threads()
    ) requires(int_alphabet) {
        if (l > r) {
            l = 0;
            r = vector.size()-1;
        }
        
        r = std::min<pos_t>(r,vector.size()-1);
        build_gmr([&vector](uint i){return vector[i];},alphabet_size,l,r);
        build_select([&vector](uint i){return vector[i];},l,r,num_threads);
    }

    /**
     * @brief builds the data structure by reading the input using the function read
     * @param read function to read the input with; it is called with i in [l,r]
     * as a parameter and must return the value of the input at index i
     * @param alphabet_size number of distinct symbols in vector (= maximum value in vector)
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     */
    rank_select_support(
        const std::function<sym_t(pos_t)>& read,
        pos_t alphabet_size,
        pos_t l = 1, pos_t r = 0,
        uint16_t num_threads = omp_get_max_threads()
    ) requires(int_alphabet) {
        build_gmr(read,alphabet_size,l,r);
        build_select(read,l,r,num_threads);
    }

    /**
     * @brief returns the value at index i in the input
     * @param i [0..size-1]
     * @return the value at index i in the input
     */
    inline sym_t operator[](pos_t i) const {
        if constexpr (byte_alphabet) {
            for (i_sym_t c : alphabet) {
                if (hyb_bit_vecs[c][i] == 1) {
                    return c;
                }
            }

            return 0;
        } else {
            return gmr[i];
        }
    }

    /**
     * @brief returns the size of the input
     * @return the size of the input
     */
    inline pos_t size() const {
        if constexpr (byte_alphabet) {
            return input_size;
        } else {
            return gmr.size();
        }
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
        if constexpr (byte_alphabet) {
            uint64_t size = sizeof(pos_t)+alphabet.size()+sizeof(pos_t)*alphabet.size();
            size += c_array.size_in_bytes();
            size += occs.size_in_bytes();

            for (uint16_t i=0; i<hyb_bit_vecs.size(); i++) {
                size += hyb_bit_vecs[i].size_in_bytes();
            }

            return size;
        } else {
            return sdsl::size_in_bytes(gmr)+
                c_array.size_in_bytes()+
                occs.size_in_bytes();
        }
    }

    /**
     * @brief returns the number of distinct symbols in the input
     * @return the number of distinct symbols in the input
     */
    inline pos_t alphabet_size() const {
        return sigma;
    }

    inline bool contains(sym_t v) const {
        if constexpr (byte_alphabet) {
            return freq[symbol_idx(v)] != 0;
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
            return c_array[v+1]-c_array[v];
        }
    }

    /**
     * @brief returns the number of occurrences of c before index i
     * @param v a value that occurs in the input
     * @param i [0..string size]
     * @return number of occurrences of c before index i
     */
    inline pos_t rank(sym_t v, pos_t i) const {
        if constexpr (byte_alphabet) {
            return hyb_bit_vecs[symbol_idx(v)].rank_1(i);
        } else {
            return gmr.rank(i,v);
        }
    }
    
    /**
     * @brief returns the index of the i-th occurrence of c
     * @param v a value that occurs in the input
     * @param i [1..number of occurrences of c in the input]
     * @return index of the i-th occurrence of c
     */
    inline pos_t select(sym_t v, pos_t i) const {
        if constexpr (byte_alphabet) {
            return occs[c_array.template get_unsafe<0,pos_t>(v)+i-1];
        } else {
            return occs[c_array[v]+i-1];
        }
    }

    /**
     * @brief serializes the rank_select_support to an output stream
     * @param out output stream
     */
    void serialize(std::ostream& out) const {
        if constexpr (byte_alphabet) {
            out.write((char*)&input_size,sizeof(pos_t));
            
            if (input_size > 0) {
                uint8_t num_chars = alphabet.size();
                out.write((char*)&num_chars,1);

                out.write((char*)&alphabet[0],alphabet.size());
                out.write((char*)&freq[0],256*sizeof(pos_t));

                for (uint16_t i=0; i<alphabet.size(); i++) {
                    hyb_bit_vecs[alphabet[i]].serialize(out);
                }
            }
        } else {
            pos_t vector_size = gmr.size();
            out.write((char*)&vector_size,sizeof(pos_t));
            out.write((char*)&sigma,sizeof(pos_t));
            
            if (vector_size > 0) {
                gmr.serialize(out);
            }
        }

        c_array.serialize(out);
        occs.serialize(out);
    }

    /**
     * @brief loads the rank_select_support from an input stream
     * @param in input stream
     */
    void load(std::istream& in) {
        if constexpr (byte_alphabet) {
            in.read((char*)&input_size,sizeof(pos_t));

            if (input_size > 0) {
                uint8_t num_chars;
                in.read((char*)&num_chars,1);

                alphabet.resize(num_chars);
                in.read((char*)&alphabet[0],num_chars);

                freq.resize(256);
                in.read((char*)&freq[0],256*sizeof(pos_t));

                hyb_bit_vecs.resize(256);

                for (uint16_t i=0; i<alphabet.size(); i++) {
                    hyb_bit_vecs[alphabet[i]].load(in);
                }
            }
        } else {
            pos_t vector_size;
            in.read((char*)&vector_size,sizeof(pos_t));
            in.read((char*)&sigma,sizeof(pos_t));

            if (vector_size > 0) {
                gmr.load(in);
            }
        }

        c_array.load(in);
        occs.load(in);
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