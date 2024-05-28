#pragma once

#include <vector>
#include <iostream>
#include <filesystem>
#include <sdsl/wt_gmr.hpp>
#include <move_r/data_structures/hybrid_bit_vector.hpp>

/**
 * @brief a string rank-select data structure using hybrid bit vectors (either sd_array or plain bit vector)
 * @tparam pos_t unsigned integer type
 */
template <typename sym_t, typename pos_t = uint32_t>
class rank_select_support {
    static_assert(std::is_same<pos_t,uint32_t>::value || std::is_same<pos_t,uint64_t>::value);
    
    static_assert(
        std::is_same<sym_t,char>::value ||
        std::is_same<sym_t,uint8_t>::value ||
        std::is_same<sym_t,uint16_t>::value ||
        std::is_same<sym_t,uint32_t>::value ||
        std::is_same<sym_t,uint64_t>::value
    );
    
    using hybrid_bv_t = hybrid_bit_vector<pos_t,true,false,true>;

    protected:
    pos_t string_size = 0; // the size of the input
    pos_t sigma = 0; // the number of distinct symbols in the input
    std::vector<char> characters; // vector containing exactly the characters occurring in the input (only used if sym_t = char)
    std::vector<pos_t> frequencies; // [0..255/sigma-1] (for sym_t = char/.) stores at position c the frequency of c in the input

    /**
     * @brief hyb_bit_vecs[c] contains a hybrid bit vector [0..size-1] that marks (with ones) the occurrences of c
     *        in the input, for each value c occurring in the input
     */
    std::vector<hybrid_bv_t> hyb_bit_vecs;

    /**
     * @brief rank-select data structure for integer alphabets (for sym_t != char)
     */
    sdsl::wt_gmr_rs<> gmr;

    pos_t symbol_idx(sym_t sym) const {
        if constexpr (std::is_same<sym_t,char>::value) {
            return char_to_uchar(sym);
        } else {
            return sym;
        }
    };

    public:
    rank_select_support() = default;

    /**
     * @brief builds the data structure from the range [l,r] in an input (0 <= l <= r < string.size()),
     * else if l > r, then the the data structure is built for the whole input
     * @param string the input
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     * @param p the number of threads to use
     */
    rank_select_support(const std::string& string, pos_t l = 1, pos_t r = 0, uint16_t p = omp_get_max_threads()) requires(std::is_same<sym_t,char>::value) {
        if (l > r) {
            l = 0;
            r = string.size()-1;
        }
        
        r = std::min(r,string.size()-1);
        *this = std::move(rank_select_support([&string](uint i){return string[i];},l,r,p));
    }
    
    /**
     * @brief builds the data structure by reading the input using the function read
     * @param read function to read the input with; it is called with i in [l,r]
     * as a parameter and must return the value of the input at index i
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     * @param p the number of threads to use
     */
    rank_select_support(const std::function<char(pos_t)>& read, pos_t l = 1, pos_t r = 0, uint16_t p = omp_get_max_threads()) requires(std::is_same<sym_t,char>::value) {
        string_size = r-l+1;
        std::vector<std::vector<pos_t>> occurrences_thr(p,std::vector<pos_t>(256,0));

        #pragma omp parallel for num_threads(p)
        for (uint64_t i=l; i<=r; i++) {
            occurrences_thr[omp_get_thread_num()][char_to_uchar(read(i))]++;
        }

        frequencies.resize(256);

        for (uint16_t cur_char=0; cur_char<256; cur_char++) {
            for (uint16_t i_p=0; i_p<p; i_p++) {
                frequencies[cur_char] += occurrences_thr[i_p][cur_char];
            }

            if (frequencies[cur_char] != 0) {
                characters.emplace_back(uchar_to_char((uint8_t)cur_char));
            }
        }

        occurrences_thr.clear();
        occurrences_thr.shrink_to_fit();
        sigma = characters.size();
        uint16_t num_chars = sigma;
        std::sort(characters.begin(),characters.end(),[this](char a, char b){return frequencies[char_to_uchar(a)] > frequencies[char_to_uchar(b)];});
        std::vector<uint8_t> is_compressed(256,1);
        std::vector<sdsl::sd_vector_builder> sdv_builders(256);
        std::vector<sdsl::bit_vector> plain_bvs(256);

        #pragma omp parallel for num_threads(p)
        for (uint16_t i=0; i<num_chars; i++) {
            uint8_t cur_char = char_to_uchar(characters[i]);

            if (frequencies[cur_char] > string_size * hybrid_bit_vector<pos_t>::compression_threshold) {
                is_compressed[cur_char] = 0;
                plain_bvs[cur_char] = std::move(sdsl::bit_vector(string_size));
            } else {
                sdv_builders[cur_char] = std::move(sdsl::sd_vector_builder(string_size,frequencies[cur_char]));
            }
        }

        uint8_t cur_char;

        for (uint64_t i=l; i<=r; i++) {
            cur_char = char_to_uchar(read(i));

            if (is_compressed[cur_char]) {
                sdv_builders[cur_char].set(i-l);
            } else {
                plain_bvs[cur_char][i-l] = 1;
            }
        }

        hyb_bit_vecs.resize(256);

        #pragma omp parallel for num_threads(p)
        for (uint16_t i=0; i<num_chars; i++) {
            uint8_t cur_char = char_to_uchar(characters[i]);

            if (is_compressed[cur_char]) {
                hyb_bit_vecs[cur_char] = std::move(hybrid_bv_t(std::move(sdsl::sd_vector<>(sdv_builders[cur_char]))));
            } else {
                hyb_bit_vecs[cur_char] = std::move(hybrid_bv_t(std::move(plain_bvs[cur_char])));
            }
        }
    }

    /**
     * @brief builds the data structure by reading the input vector
     * @param vector input vector
     * @param alphabet_size number of distinct symbols in vector (= maximum value in vector)
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     * @param p the number of threads to use
     */
    rank_select_support(const std::vector<sym_t>& vector, pos_t alphabet_size, pos_t l = 1, pos_t r = 0) requires(!std::is_same<sym_t,char>::value) {
        if (l > r) {
            l = 0;
            r = vector.size()-1;
        }
        
        r = std::min(r,vector.size()-1);
        *this = std::move(rank_select_support([&vector](uint i){return vector[i];},alphabet_size,l,r));
    }

    /**
     * @brief builds the data structure by reading the input using the function read
     * @param read function to read the input with; it is called with i in [l,r]
     * as a parameter and must return the value of the input at index i
     * @param alphabet_size number of distinct symbols in vector (= maximum value in vector)
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     */
    rank_select_support(const std::function<sym_t(pos_t)>& read, pos_t alphabet_size, pos_t l = 1, pos_t r = 0) requires(!std::is_same<sym_t,char>::value) {
        sigma = alphabet_size;
        std::string vector_buffer_filename = "vec_" + random_alphanumeric_string(10);
        sdsl::int_vector_buffer<> vector_buffer(vector_buffer_filename,std::ios::out);
        frequencies.resize(alphabet_size,0);

        for (uint64_t i=l; i<=r; i++) {
            vector_buffer[i-l] = read(i);
            frequencies[symbol_idx(read(i-l))]++;
        }

        gmr = std::move(sdsl::wt_gmr_rs<>(vector_buffer,r-l+1));
        vector_buffer.close();
        std::filesystem::remove(vector_buffer_filename);
    }

    /**
     * @brief returns the value at index i in the input
     * @param i [0..size-1]
     * @return the value at index i in the input
     */
    inline sym_t operator[](pos_t i) const {
        if constexpr (std::is_same<sym_t,char>::value) {
            for (char c : characters) {
                if (hyb_bit_vecs[char_to_uchar(c)][i] == 1) {
                    return c;
                }
            }

            return uchar_to_char(0);
        } else {
            return gmr[i];
        }
    }

    /**
     * @brief returns the size of the input
     * @return the size of the input
     */
    inline pos_t size() const {
        if constexpr (std::is_same<sym_t,char>::value) {
            return string_size;
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
        if constexpr (std::is_same<sym_t,char>::value) {
            uint64_t size = sizeof(pos_t)+characters.size()+sizeof(pos_t)*characters.size();

            for (uint16_t i=0; i<hyb_bit_vecs.size(); i++) {
                size += hyb_bit_vecs[i].size_in_bytes();
            }

            return size;
        } else {
            return sdsl::size_in_bytes(gmr)+sigma*sizeof(pos_t);
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
        if constexpr (std::is_same<sym_t,char>::value) {
            return frequencies[symbol_idx(v)] > 0;
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
        return frequencies[symbol_idx(v)];
    }

    /**
     * @brief returns the number of occurrences of c before index i
     * @param v a value that occurs in the input
     * @param i [0..string size]
     * @return number of occurrences of c before index i
     */
    inline pos_t rank(sym_t v, pos_t i) const {
        if constexpr (std::is_same<sym_t,char>::value) {
            return hyb_bit_vecs[char_to_uchar(v)].rank_1(i);
        } else {
            return gmr.rank(i,v);
        }
    }
    
    /**
     * @brief returns the index of the i-th occurrence of c
     * @param v a value that occurs in the input
     * @param i [0..number of occurrences of c in the input - 1]
     * @return index of the i-th occurrence of c
     */
    inline pos_t select(sym_t v, pos_t i) const {
        if constexpr (std::is_same<sym_t,char>::value) {
            return hyb_bit_vecs[char_to_uchar(v)].select_1(i);
        } else {
            return gmr.select(i,v);
        }
    }

    /**
     * @brief serializes the rank_select_support to an output stream
     * @param out output stream
     */
    void serialize(std::ostream& out) const {
        if constexpr (std::is_same<sym_t,char>::value) {
            out.write((char*)&string_size,sizeof(pos_t));
            
            if (string_size > 0) {
                uint8_t num_chars = characters.size();
                out.write((char*)&num_chars,1);

                out.write((char*)&characters[0],characters.size());
                out.write((char*)&frequencies[0],256*sizeof(pos_t));

                for (uint16_t i=0; i<characters.size(); i++) {
                    hyb_bit_vecs[char_to_uchar(characters[i])].serialize(out);
                }
            }
        } else {
            pos_t vector_size = gmr.size();
            out.write((char*)&vector_size,sizeof(pos_t));
            out.write((char*)&sigma,sizeof(pos_t));
            out.write((char*)&frequencies[0],sigma*sizeof(pos_t));
            
            if (vector_size > 0) {
                gmr.serialize(out);
            }
        }
    }

    /**
     * @brief loads the rank_select_support from an input stream
     * @param in input stream
     */
    void load(std::istream& in) {
        if constexpr (std::is_same<sym_t,char>::value) {
            in.read((char*)&string_size,sizeof(pos_t));

            if (string_size > 0) {
                uint8_t num_chars;
                in.read((char*)&num_chars,1);

                characters.resize(num_chars);
                in.read((char*)&characters[0],num_chars);

                frequencies.resize(256);
                in.read((char*)&frequencies[0],256*sizeof(pos_t));

                hyb_bit_vecs.resize(256);

                for (uint16_t i=0; i<characters.size(); i++) {
                    hyb_bit_vecs[char_to_uchar(characters[i])].load(in);
                }
            }
        } else {
            pos_t vector_size;
            in.read((char*)&vector_size,sizeof(pos_t));
            in.read((char*)&sigma,sizeof(pos_t));
            in.read((char*)&frequencies[0],sigma*sizeof(pos_t));

            if (vector_size > 0) {
                gmr.load(in);
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