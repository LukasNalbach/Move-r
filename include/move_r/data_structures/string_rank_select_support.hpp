#pragma once

#include <vector>
#include <iostream>
#include <move_r/data_structures/hybrid_bit_vector.hpp>

/**
 * @brief a string rank-select data structure using hybrid bit vectors (either sd_array or plain bit vector)
 * @tparam uint_t unsigned integer type
 */
template <typename uint_t = uint32_t, bool build_rank_support = true, bool build_select_support = true, bool build_select_other_support = false>
class string_rank_select_support {
    static_assert(std::is_same<uint_t,uint32_t>::value || std::is_same<uint_t,uint64_t>::value);
    
    using hybrid_bv_t = hybrid_bit_vector<uint_t,build_rank_support,build_select_other_support,build_select_support>;

    protected:
    uint_t size = 0; // the size of the input string    
    std::vector<char> characters; // vector containing exactly the characters occurring in the input string
    std::vector<uint_t> occurrences; // [0..255] stores at position c the number of occurrences of c in the input string

    /**
     * @brief hyb_bit_vecs[c] contains a hybrid bit vector [0..size-1] that marks (with ones) the occurrences of c
     *        in the input string, for each character c occurring in the input string
     */
    std::vector<hybrid_bv_t> hyb_bit_vecs;

    public:
    string_rank_select_support() = default;

    /**
     * @brief builds the data structure from the range [l,r] in an input string (0 <= l <= r < string.size()),
     * else if l > r, then the the data structure is built for the whole input string
     * @param string the input string
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     * @param p the number of threads to use
     */
    string_rank_select_support(
        const std::string& string,
        uint_t l = 1,
        uint_t r = 0,
        uint16_t p = omp_get_max_threads()
    ) {
        r = std::max(r,string.size()-1);

        if (l > r) {
            l = 0;
            r = string.size()-1;
        }

        *this = std::move(string_rank_select_support([&string](uint i){return string[i];},l,r,p));
    }
    
    /**
     * @brief builds the data structure by reading the input using the function read
     * @param read function to read the input string with; it is called with i in [l,r]
     * as a parameter and must return the character of the input string at index i
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     * @param p the number of threads to use
     */
    string_rank_select_support(
        const std::function<char(uint_t)>& read,
        uint_t l = 1,
        uint_t r = 0,
        uint16_t p = omp_get_max_threads()
    ) {
        size = r-l+1;
        std::vector<std::vector<uint_t>> occurrences_thr(p,std::vector<uint_t>(256,0));

        #pragma omp parallel for num_threads(p)
        for (uint64_t i=l; i<=r; i++) {
            occurrences_thr[omp_get_thread_num()][char_to_uchar(read(i))]++;
        }

        occurrences.resize(256);

        for (uint16_t cur_char=0; cur_char<256; cur_char++) {
            for (uint16_t i_p=0; i_p<p; i_p++) {
                occurrences[cur_char] += occurrences_thr[i_p][cur_char];
            }

            if (occurrences[cur_char] != 0) {
                characters.emplace_back(uchar_to_char((uint8_t)cur_char));
            }
        }

        occurrences_thr.clear();
        occurrences_thr.shrink_to_fit();
        uint16_t num_chars = characters.size();
        std::vector<uint8_t> is_compressed(256,1);
        std::vector<sdsl::sd_vector_builder> sdv_builders(256);
        std::vector<sdsl::bit_vector> plain_bvs(256);

        #pragma omp parallel for num_threads(p)
        for (uint16_t i=0; i<num_chars; i++) {
            uint8_t cur_char = char_to_uchar(characters[i]);

            if (occurrences[cur_char] > size * hybrid_bit_vector<uint_t>::compression_threshold) {
                is_compressed[cur_char] = 0;
                plain_bvs[cur_char] = std::move(sdsl::bit_vector(size));
            } else {
                sdv_builders[cur_char] = std::move(sdsl::sd_vector_builder(size,occurrences[cur_char]));
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
     * @brief returns the number of distinct characters in the input string (alphabet size)
     * @return alphabet_size 
     */
    inline uint8_t alphabet_size() {
        return characters.size();
    }

    /**
     * @brief returns the number of occurrences of c in the input string
     * @param c character
     * @return the number of occurrences of c in the input string
     */
    inline uint_t num_occurrences(char c) {
        return occurrences[char_to_uchar(c)];
    }

    /**
     * @brief returns the size of the data structure in bytes
     * @return size of the data structure in bytes
     */
    uint64_t size_in_bytes() {
        uint64_t size = sizeof(uint_t)+characters.size()+sizeof(uint_t)*characters.size();

        for (uint16_t i=0; i<hyb_bit_vecs.size(); i++) {
            size += hyb_bit_vecs[i].size_in_bytes();
        }

        return size;
    }

    /**
     * @brief returns the number of occurrences of c before index i
     * @param c a character that occurs in the input string
     * @param i [0..size]
     * @return number of occurrences of c before index i
     */
    inline uint_t rank(char c, uint_t i) {
        static_assert(build_rank_support);
        return hyb_bit_vecs[char_to_uchar(c)].rank_1(i);
    }

    /**
     * @brief returns the number of occurrences of characters other than c before index i
     * @param c a character that occurs in the input string
     * @param i [0..size]
     * @return number of occurrences of characters other than c before index i
     */
    inline uint_t rank_other(char c, uint_t i) {
        static_assert(build_rank_support);
        return hyb_bit_vecs[char_to_uchar(c)].rank_0(i);
    }
    
    /**
     * @brief returns the index of the i-th occurrence of c
     * @param c a character that occurs in the input string
     * @param i [0..number of occurrences of c in the input string - 1]
     * @return index of the i-th occurrence of c
     */
    inline uint_t select(char c, uint_t i) {
        static_assert(build_select_support);
        return hyb_bit_vecs[char_to_uchar(c)].select_1(i);
    }

    /**
     * @brief returns the i-th index, at which there is a character other than c
     * @param c a character that occurs in the input string
     * @param i [0..sum of the numbers of occurrences of all characters other than c in the input string - 1]
     * @return i-th index, at which there is a character other than c
     */
    inline uint_t select_other(char c, uint_t i) {
        static_assert(build_select_other_support);
        return hyb_bit_vecs[char_to_uchar(c)].select_0(i);
    }

    /**
     * @brief returns the index of the next occurrence of c after index i
     * @param c a character that occurs after index i in the input string
     * @param i [0..size-1]
     * @return the index of the next occurrence of c after index i
     */
    inline uint_t next(char c, uint_t i) {
        return hyb_bit_vecs[char_to_uchar(c)].next_1(i);
    }

    /**
     * @brief returns the index of the previous occurrence of c before index i
     * @param c a character that occurs before index i in the input string
     * @param i [0..size-1]
     * @return the index of the next occurrence of c before index i
     */
    inline uint_t previous(char c, uint_t i) {
        return hyb_bit_vecs[char_to_uchar(c)].previous_1(i);
    }

    /**
     * @brief returns the index of the next character other than c after index i
     * @param c a character that occurs in the input string
     * @param i [0..size-1] an index in the input string, after which
     * there occurs a character c' != c in the input string
     * @return the index of the next character other than c after index i 
     */
    inline uint_t next_other(char c, uint_t i) {
        return hyb_bit_vecs[char_to_uchar(c)].next_0(i);
    }

    /**
     * @brief returns the index of the previous character other than c before index i
     * @param c a character that occurs in the input string
     * @param i [0..size-1] an index in the input string, before which
     * there occurs a character c' != c in the input string
     * @return the index of the previous character other than c before index i
     */
    inline uint_t previous_other(char c, uint_t i) {
        return hyb_bit_vecs[char_to_uchar(c)].previous_0(i);
    }

    /**
     * @brief returns whether there is a c at index i
     * @param i [0..size-1]
     * @param c a character that occurs in the input string
     * @return whether there is a c at index i
     */
    inline bool equals_at(uint_t i, char c) {
        return hyb_bit_vecs[char_to_uchar(c)][i];
    }

    /**
     * @brief returns whether the input string contains c
     * @param c a character
     * @return whether the input string contains c
     */
    inline bool contains_character(char c) {
        return occurrences[char_to_uchar(c)] != 0;
    }

    /**
     * @brief returns whether the character c occurs before index i in the input string
     * @param c character that occurs in the input string
     * @param i [0..size-1]
     * @return whether c occurs before index i
     */
    inline bool occurs_before(char c, uint_t i) {
        return rank(char_to_uchar(c),i) != 0;
    }

    /**
     * @brief returns whether the character c occurs after index i in the input string
     * @param c character that occurs in the input string
     * @param i [0..size-2]
     * @return whether c occurs after index i
     */
    inline bool occurs_after(char c, uint_t i) {
        return rank(char_to_uchar(c),i+1) < occurrences[c];
    }

    /**
     * @brief returns a vector containing all characters in the input string (the alphabet)
     * @return the alphabet
     */
    std::vector<char> alphabet() {
        return characters;
    }

    /**
     * @brief serializes the string_rank_select_support to an output stream
     * @param out output stream
     */
    void serialize(std::ostream& out) {
        out.write((char*)&size,sizeof(uint_t));
        
        if (size > 0) {
            uint8_t num_chars = characters.size();
            out.write((char*)&num_chars,1);

            out.write((char*)&characters[0],characters.size());
            out.write((char*)&occurrences[0],256*sizeof(uint_t));

            for (uint16_t i=0; i<characters.size(); i++) {
                hyb_bit_vecs[char_to_uchar(characters[i])].serialize(out);
            }
        }
    }

    /**
     * @brief loads the string_rank_select_support from an input stream
     * @param in input stream
     */
    void load(std::istream& in) {
        in.read((char*)&size,sizeof(uint_t));

        if (size > 0) {
            uint8_t num_chars;
            in.read((char*)&num_chars,1);

            characters.resize(num_chars);
            in.read((char*)&characters[0],num_chars);

            occurrences.resize(256);
            in.read((char*)&occurrences[0],256*sizeof(uint_t));

            hyb_bit_vecs.resize(256);

            for (uint16_t i=0; i<characters.size(); i++) {
                hyb_bit_vecs[char_to_uchar(characters[i])].load(in);
            }
        }
    }
};