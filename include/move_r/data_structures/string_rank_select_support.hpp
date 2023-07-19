#pragma once

#include <vector>
#include <iostream>
#include <move_r/data_structures/sd_array.hpp>

/**
 * @brief a string rank-select data structure using sd_arrays
 * @tparam uint_t unsigned integer type
 */
template <typename uint_t = uint32_t>
class string_rank_select_support {
	static_assert(std::is_same<uint_t,uint32_t>::value || std::is_same<uint_t,uint64_t>::value);
    
    protected:
    uint_t size = 0;
    
    std::vector<uint8_t> uchars = {}; // vector containing exactly the characters occurring in the input string
    std::vector<uint8_t> contains_uchar = {}; // [0..256] lookup table, contains_uchar[c] = 1 <=> c occurs in the input string

    /**
     * @brief sd_arrays with rank_1- and select_1 support of each character occurring in the input string;
     *        sd_arrays[c] contains an sd_vector of size [0..size-1] that marks (with ones) the occurrences of c
     *        in the input string
     */
    std::vector<sd_array<uint_t>> sd_arrays;

    public:
    string_rank_select_support() = default;
    string_rank_select_support(string_rank_select_support&& other) = default;
    string_rank_select_support(const string_rank_select_support& other) = default;
    string_rank_select_support& operator=(string_rank_select_support&& other) = default;
    string_rank_select_support& operator=(const string_rank_select_support& other) = default;

    /**
     * @brief builds the data structure from the range [l,r] in an input string (0 <= l <= r < string.size()),
     * else if l > r, then the the data structure is built for the whole input string
     * @param string the input string
	 * @param l left range limit (l <= r)
	 * @param r right range limit (l <= r)
     * @param chars vector that contains exactly the characters occurring in the input string in the range [l,r]
     * (optional, saves one scan over the input string)
     * @param p the number of threads to use
     */
    string_rank_select_support(
        const std::string& string,
        uint_t l = 1,
        uint_t r = 0,
        const std::vector<char>& chars = {},
        uint16_t p = omp_get_max_threads()
    ) {
        r = std::max(r,string.size()-1);

		if (l > r) {
			l = 0;
			r = string.size()-1;
		}

        *this = std::move(string_rank_select_support([&string](uint i){return string[i];},l,r,chars,p));
    }
    
    /**
     * @brief builds the data structure by reading the input using the function read
     * @param read function to read the input string with; it is called with i in [l,r]
     * as a parameter and must return the character of the input string at index i
	 * @param l left range limit (l <= r)
	 * @param r right range limit (l <= r)
     * @param chars vector that contains exactly the characters occurring in the input
     * string (optional, saves one scan over the input string)
     * @param p the number of threads to use
     */
    string_rank_select_support(
        const std::function<char(uint_t)>& read,
        uint_t l = 1,
        uint_t r = 0,
        const std::vector<char>& chars = {},
        uint16_t p = omp_get_max_threads()
    ) {
        size = r-l+1;

        if (chars.empty()) {
            std::vector<std::vector<uint8_t>> contains_uchar_thr(p,std::vector<uint8_t>(256,0));

            #pragma omp parallel for num_threads(p)
            for (uint64_t i=l; i<=r; i++) {
                contains_uchar_thr[omp_get_thread_num()][char_to_uchar(read(i))] = 1;
            }

            contains_uchar.resize(256,0);

            for (uint16_t i=0; i<256; i++) {
                for (uint16_t i_p=0; i_p<p; i_p++) {
                    if (contains_uchar_thr[i_p][i] == 1) {
                        contains_uchar[i] = 1;
                        break;
                    }
                }
            }

            for (uint16_t i=0; i<256; i++) {
                if (contains_uchar[i] == 1) {
                    uchars.emplace_back(i);
                }
            }
        } else {
            for (char c : chars) {
                uchars.push_back(char_to_uchar(c));
            }

            contains_uchar.resize(256,0);

            for (uint16_t i=0; i<uchars.size(); i++) {
                contains_uchar[uchars[i]] = 1;
            }
        }

        std::vector<sdsl::bit_vector> plain_bvs(256);
        uint16_t num_chars = uchars.size();

        #pragma omp parallel for num_threads(p)
        for (uint16_t i=0; i<num_chars; i++) {
            plain_bvs[uchars[i]] = std::move(sdsl::bit_vector(size));
        }

        #pragma omp parallel for num_threads(p)
        for (uint64_t i=l; i<=r; i++) {
            plain_bvs[char_to_uchar(read(i))][i] = 1;
        }

        sd_arrays.resize(256);

        #pragma omp parallel for num_threads(p)
        for (uint16_t i=0; i<num_chars; i++) {
            sd_arrays[uchars[i]] = std::move(sd_array<uint_t>(plain_bvs[uchars[i]]));
            plain_bvs[uchars[i]] = sdsl::bit_vector();
        }
    }

    /**
     * @brief returns the number of distinct characters in the input string (alphabet size)
     * @return alphabet_size 
     */
    inline uint8_t alphabet_size() {
        return uchars.size();
    }

    /**
     * @brief returns the size of the data structure in bytes
     * @return size of the data structure in bytes
     */
    uint64_t size_in_bytes() {
        uint64_t size = sizeof(uint_t)+uchars.size(); // variables

        for (uint16_t i=0; i<sd_arrays.size(); i++) {
            size += sd_arrays[i].size_in_bytes();
        }

        return size;
    }

    /**
     * @brief returns the number of occurrences of c before index i
     * @param c a character
     * @param i [0..size-1]
     * @return number of occurrences of c before index i
     */
    inline uint_t rank(uint8_t c, uint_t i) {
        return sd_arrays[c].rank_1(i);
    }

    /**
     * @brief returns the number of occurrences of c before index i
     * @param c a character
     * @param i [0..size-1]
     * @return number of occurrences of c before index i
     */
    inline uint_t rank(char c, uint_t i) {
        return rank(char_to_uchar(c),i);
    }

    /**
     * @brief returns the number of occurrences of characters other than c before index i
     * @param c a character
     * @param i [0..size-1]
     * @return number of coccurrences of characters other than c before index i
     */
    inline uint_t rank_other(uint8_t c, uint_t i) {
        return sd_arrays[c].rank_0(i);
    }

    /**
     * @brief returns the number of occurrences of characters other than c before index i
     * @param c a character
     * @param i [0..size-1]
     * @return number of coccurrences of characters other than c before index i
     */
    inline uint_t rank_other(char c, uint_t i) {
        return rank_other(char_to_uchar(c),i);
    }

    /**
     * @brief returns the index of the i-th occurrence of c
     * @param c a character that occurs in the input string
     * @param i [0..number of occurrences of c in the input string - 1]
     * @return index of the i-th occurrence of c
     */
    inline uint_t select(uint8_t c, uint_t i) {
        return sd_arrays[c].select_1(i);
    }

    /**
     * @brief returns the index of the i-th occurrence of c
     * @param c a character that occurs in the input string
     * @param i [0..number of occurrences of c in the input string - 1]
     * @return index of the i-th occurrence of c
     */
    inline uint_t select(char c, uint_t i) {
        return select(char_to_uchar(c),i);
    }

    /**
     * @brief returns the i-th index, at which there is a character other than c
     * @param c a character that occurs in the input string
     * @param i [0..sum of the numbers of occurences of all characters other than c in the input string - 1]
     * @return i-th index, at which there is a character other than c
     */
    inline uint_t select_other(uint8_t c, uint_t i) {
        return sd_arrays[c].select_1(i);
    }

    /**
     * @brief returns the i-th index, at which there is a character other than c
     * @param c a character that occurs in the input string
     * @param i [0..sum of the numbers of occurences of all characters other than c in the input string - 1]
     * @return i-th index, at which there is a character other than c
     */
    inline uint_t select_other(char c, uint_t i) {
        return select_other(char_to_uchar(c),i);
    }

    /**
     * @brief returns the index of the next occurence of c after index i
     * @param c a character that occurs in the input string
     * @param i [0..size-1]
     * @return the index of the next occurence of c after index i
     */
    inline uint_t next(uint8_t c, uint_t i) {
        return sd_arrays[c].next_1(i);
    }

    /**
     * @brief returns the index of the next occurence of c after index i
     * @param c a character that occurs in the input string
     * @param i [0..size-1]
     * @return the index of the next occurence of c after index i
     */
    inline uint_t next(char c, uint_t i) {
        return next(char_to_uchar(c),i);
    }

    /**
     * @brief returns the index of the previous occurence of c before index i
     * @param c a character that occurs in the input string
     * @param i [0..size-1]
     * @return the index of the next occurence of c before index i
     */
    inline uint_t previous(uint8_t c, uint_t i) {
        return sd_arrays[c].previous_1(i);
    }

    /**
     * @brief returns the index of the previous occurence of c before index i
     * @param c a character that occurs in the input string
     * @param i [0..size-1]
     * @return the index of the next occurence of c before index i
     */
    inline uint_t previous(char c, uint_t i) {
        return previous(char_to_uchar(c),i);
    }

    /**
     * @brief returns the index of the next character other than c after index i
     * @param c a character that occurs in the input string
     * @param i [0..size-1]
     * @return the index of the next character other than c after index i 
     */
    inline uint_t next_other(uint8_t c, uint_t i) {
        return sd_arrays[c].next_0(i);
    }

    /**
     * @brief returns the index of the next character other than c after index i
     * @param c a character that occurs in the input string
     * @param i [0..size-1]
     * @return the index of the next character other than c after index i 
     */
    inline uint_t next_other(char c, uint_t i) {
        return next_other(char_to_uchar(c),i);
    }

    /**
     * @brief returns the index of the previous character other than c before index i
     * @param c a character that occurs in the input string
     * @param i [0..size-1]
     * @return the index of the previous character other than c before index i
     */
    inline uint_t previous_other(uint8_t c, uint_t i) {
        return sd_arrays[c].previous_0(i);
    }

    /**
     * @brief returns the index of the previous character other than c before index i
     * @param c a character that occurs in the input string
     * @param i [0..size-1]
     * @return the index of the previous character other than c before index i
     */
    inline uint_t previous_other(char c, uint_t i) {
        return previous_other(char_to_uchar(c),i);
    }

    /**
     * @brief returns whether there is a c at index i
     * @param i [0..size-1]
     * @param c a character that occurs in the input string
     * @return whether there is a c at index i
     */
    inline bool equals_at(uint_t i, uint8_t c) {
        return sd_arrays[c][i];
    }

    /**
     * @brief returns whether there is a c at index i
     * @param i [0..size-1]
     * @param c a character that occurs in the input string
     * @return whether there is a c at index i
     */
    inline bool equals_at(uint_t i, char c) {
        return equals_at(i,char_to_uchar(c));
    }

    /**
     * @brief whether the input string contains c
     * @param c a character
     * @return whether the input string contains c
     */
    inline bool contains_character(uint8_t c) {
        return contains_uchar[c] != 0;
    }

    /**
     * @brief whether the input string contains c
     * @param c a character
     * @return whether the input string contains c
     */
    inline bool contains_character(char c) {
        return contains_uchar[char_to_uchar(c)] != 0;
    }

    /**
     * @brief returns a vector containing all characters in the input string (the alphabet)
     * @tparam char_t character type
     * @return the alphabet
     */
    template <typename char_t>
    std::vector<char_t> alphabet() {
	    static_assert(std::is_same<char_t,uint8_t>::value || std::is_same<char_t,char>::value);

        if constexpr (std::is_same<char_t,uint8_t>::value) {
            return uchars;
        } else {
            std::vector<char> chars(0);

            for (uint8_t c : uchars) {
                chars.emplace_back(uchar_to_char(c));
            }
            
            return chars;
        }
    }

    /**
     * @brief serializes the string_rank_select_support to an output stream
     * @param out output stream
     */
    void serialize(std::ostream& out) {
        out.write((char*)&size,sizeof(uint_t));
        
        if (size > 0) {
            uint8_t num_chars = uchars.size();
            out.write((char*)&num_chars,1);

            out.write((char*)&uchars[0],uchars.size());
            out.write((char*)&contains_uchar[0],256);

            for (uint16_t i=0; i<uchars.size(); i++) {
                sd_arrays[uchars[i]].serialize(out);
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

            uchars.resize(num_chars);
            in.read((char*)&uchars[0],num_chars);

            contains_uchar.resize(256);
            in.read((char*)&contains_uchar[0],256);

            sd_arrays.resize(256);

            for (uint16_t i=0; i<uchars.size(); i++) {
                sd_arrays[uchars[i]].load(in);
            }
        }
    }
};