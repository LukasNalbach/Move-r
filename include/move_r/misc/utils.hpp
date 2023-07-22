#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <climits>
#include <string>
#include <unistd.h>
#include <malloc_count.h>

uint64_t ram_size() {
    uint64_t pages = sysconf(_SC_PHYS_PAGES);
    uint64_t page_size = sysconf(_SC_PAGE_SIZE);
    return pages*page_size;
}

std::chrono::steady_clock::time_point now() {
    return std::chrono::steady_clock::now();
}

std::string format_time(uint64_t ns) {
    std::string time_str;

    if (ns > 10000000000) {
        time_str = std::to_string(ns/1000000000) + " s";
    } else if (ns > 10000000) {
        time_str = std::to_string(ns/1000000) + " ms";
    } else if (ns > 10000) {
        time_str = std::to_string(ns/1000) + " us";
    } else {
        time_str = std::to_string(ns) + " ns";
    }

    return time_str;
}

std::string format_query_throughput(uint64_t num_queries, uint64_t ns) {
    std::string str;
    double queries_per_ns = num_queries/(double)ns;

    if (queries_per_ns < 0.000001) {
        str = std::to_string(queries_per_ns*1000000000) + " queries/s";
    } else if (queries_per_ns < 0.001) {
        str = std::to_string(queries_per_ns*1000000) + " queries/ms";
    } else if (queries_per_ns < 1) {
        str = std::to_string(queries_per_ns*1000) + " queries/us";
    } else {
        str = std::to_string(queries_per_ns) + " queries/ns";
    }

    return str;
}

std::string format_size(uint64_t B) {
    std::string size_str;

    if (B > 10000000000) {
        size_str = std::to_string(B/1000000000) + " GB";
    } else if (B > 10000000) {
        size_str = std::to_string(B/1000000) + " MB";
    } else if (B > 10000) {
        size_str = std::to_string(B/1000) + " KB";
    } else {
        size_str = std::to_string(B) + " B";
    }
    
    return size_str;
}

std::string format_threads(uint16_t p) {
    if (p == 1) {
        return "1 thread";
    } else {
        return std::to_string(p) + " threads";
    }
}

uint64_t time_diff_ns(std::chrono::steady_clock::time_point t1, std::chrono::steady_clock::time_point t2) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
}

uint64_t time_diff_ns(std::chrono::steady_clock::time_point t) {
    return time_diff_ns(t,std::chrono::steady_clock::now());
}

std::chrono::steady_clock::time_point log_runtime(std::chrono::steady_clock::time_point t1, std::chrono::steady_clock::time_point t2) {
    std::cout << ", in ~ " << format_time(time_diff_ns(t1,t2)) << std::endl;
    return std::chrono::steady_clock::now();
}

std::chrono::steady_clock::time_point log_runtime(std::chrono::steady_clock::time_point t) {
    return log_runtime(t,std::chrono::steady_clock::now());
}

void log_peak_mem_usage(uint64_t baseline_memory_allocation) {
    std::cout << "peak memory allocation until now: " << format_size(malloc_count_peak() - baseline_memory_allocation) << std::endl;
}

void log_mem_usage(uint64_t baseline_memory_allocation) {
    std::cout << "current memory allocation: " << format_size(malloc_count_current() - baseline_memory_allocation) << std::endl;
}

void log_message(std::string message) {
    std::cout << message << std::flush;
}

void print_header_error() {
    std::cout << "Error: malformed header in patterns file" << std::endl;
    std::cout << "Take a look here for more info on the file format: http://pizzachili.dcc.uchile.cl/experiments.html" << std::endl;
    exit(0);
}

uint64_t number_of_patterns(std::string header) {
    uint64_t start_pos = header.find("number=");

    if (start_pos == std::string::npos or start_pos + 7 >= header.size()) {
        print_header_error();
    }

    start_pos += 7;
    uint64_t end_pos = header.substr(start_pos).find(" ");

    if (end_pos == std::string::npos) {
        print_header_error();
    }

    return std::atoi(header.substr(start_pos).substr(0, end_pos).c_str());
}

uint64_t patterns_length(std::string header) {
    uint64_t start_pos = header.find("length=");

    if (start_pos == std::string::npos or start_pos + 7 >= header.size()) {
        print_header_error();
    }

    start_pos += 7;
    uint64_t end_pos = header.substr(start_pos).find(" ");

    if (end_pos == std::string::npos) {
        print_header_error();
    }

    return std::atoi(header.substr(start_pos).substr(0, end_pos).c_str());
}

void read_from_file(std::istream& in, const char* data, uint64_t size) {
    uint64_t size_left = size;
    uint64_t bytes_to_read;

    while (size_left > 0) {
        bytes_to_read = std::min(size_left,(uint64_t)INT_MAX);
        in.read((char*)&data[size-size_left],bytes_to_read);
        size_left -= bytes_to_read;
    }
}

void write_to_file(std::ostream& out, const char* data, uint64_t size) {
    uint64_t size_left = size;
    uint64_t bytes_to_write;

    while (size_left > 0) {
        bytes_to_write = std::min(size_left,(uint64_t)INT_MAX);
        out.write((char*)&data[size-size_left],bytes_to_write);
        size_left -= bytes_to_write;
    }
}

/**
 * @brief interprets a uint8_t as a char
 * @param c a uint8_t character
 * @return the same character interpreted as a char 
 */
inline char uchar_to_char(uint8_t c) {
    return *reinterpret_cast<char*>(&c);
}

/**
 * @brief interprets a char as a char
 * @param c a char character
 * @return the same character interpreted as a uint8_t 
 */
inline uint8_t char_to_uchar(char c) {
    return *reinterpret_cast<uint8_t*>(&c);
}

template <typename T>
bool contains(std::vector<T> vec, T val) {
    return std::find(vec.begin(),vec.end(),val) != vec.end();
}

template <typename T>
bool is_subset_of(std::vector<T> first, std::vector<T> second) {
    std::sort(first.begin(), first.end());
    std::sort(second.begin(), second.end());
    return std::includes(second.begin(),second.end(),first.begin(),first.end());
}

std::string random_alphanumeric_string(uint64_t length) {
    static std::string possible_chars = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

    std::string str_rand;
    str_rand.reserve(length);
    
    for (uint64_t i=0; i<length; i++) {
        str_rand.push_back(possible_chars[std::rand()%possible_chars.size()]);
    }

    return str_rand;
}

template <typename uint_t = uint32_t>
std::string to_string(std::pair<uint_t,uint_t> pair) {
    std::string str = "(";
    str.append(std::to_string(pair.first));
    str.push_back(',');
    str.append(std::to_string(pair.second));
    str.push_back(')');
    return str;
}

template<typename T>
class no_init {
    static_assert(std::is_fundamental<T>::value);

    private:
    T v_;

    public:
    no_init() noexcept {}
    constexpr no_init(T value) noexcept: v_{value} {}
    constexpr operator T() const noexcept {return v_;}
};

template <typename T, typename Alloc = std::allocator<T>>
class default_init_allocator : public Alloc {
    using a_t = std::allocator_traits<Alloc>;
    using Alloc::Alloc;

    public:
    template<typename U> struct rebind {};
    template<typename U> void construct (U* ptr) noexcept(std::is_nothrow_default_constructible<U>::value) {::new(static_cast<void*>(ptr)) U;}
    template<typename U, typename... Args> void construct (U* ptr, Args&&... args) {}
};

void no_init_resize(std::string& str, size_t size) {
    (*reinterpret_cast<std::basic_string<char,std::char_traits<char>,default_init_allocator<char>>*>(&str)).resize(size);
}