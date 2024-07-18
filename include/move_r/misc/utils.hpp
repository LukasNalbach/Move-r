#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <climits>
#include <string>
#include <functional>
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

uint64_t time_diff_min(std::chrono::steady_clock::time_point t1, std::chrono::steady_clock::time_point t2) {
    return std::chrono::duration_cast<std::chrono::minutes>(t2-t1).count();
}

uint64_t time_diff_ns(std::chrono::steady_clock::time_point t1, std::chrono::steady_clock::time_point t2) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
}

uint64_t time_diff_min(std::chrono::steady_clock::time_point t) {
    return time_diff_min(t,std::chrono::steady_clock::now());
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

inline char uchar_to_char(uint8_t c) {
    return *reinterpret_cast<char*>(&c);
}

inline uint8_t char_to_uchar(char c) {
    return *reinterpret_cast<uint8_t*>(&c);
}

template <typename T>
bool contains(std::vector<T> vec, T val) {
    return std::find(vec.begin(),vec.end(),val) != vec.end();
}

template <typename T>
void remove(std::vector<T> vec, T val) {
    while (contains(vec,val)) {
        vec.erase(std::find(vec.begin(),vec.end(),val));
    }
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

template <typename pos_t = uint32_t>
std::string to_string(std::pair<pos_t,pos_t> pair) {
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

template <typename T>
void no_init_resize(std::vector<T>& vec, size_t size) {
    (*reinterpret_cast<std::vector<no_init<T>>*>(&vec)).resize(size);
}

template <typename T1, typename T2>
void no_init_resize(std::vector<std::pair<T1,T2>>& vec, size_t size) {
    (*reinterpret_cast<std::vector<std::pair<no_init<T1>,no_init<T2>>>*>(&vec)).resize(size);
}

template <typename T>
void no_init_resize(std::vector<std::tuple<T,T,T>>& vec, size_t size) {
    (*reinterpret_cast<std::vector<std::tuple<no_init<T>,no_init<T>,no_init<T>>>*>(&vec)).resize(size);
}

template <typename pos_t>
pos_t bin_search_max_leq(pos_t value, pos_t left, pos_t right, std::function<pos_t(pos_t)> value_at) {
    pos_t middle;

    while (left != right) {
        middle = left+(right-left)/2+1;

        if (value_at(middle) <= value) {
            left = middle;
        } else {
            right = middle-1;
        }
    }

    return left;
}

template <typename pos_t>
pos_t bin_search_min_geq(pos_t value, pos_t left, pos_t right, std::function<pos_t(pos_t)> value_at) {
    pos_t middle;

    while (left != right) {
        middle = left+(right-left)/2;

        if (value <= value_at(middle)) {
            right = middle;
        } else {
            left = middle+1;
        }
    }

    return left;
}

template <typename pos_t>
pos_t bin_search_max_lt(pos_t value, pos_t left, pos_t right, std::function<pos_t(pos_t)> value_at) {
    pos_t middle;

    while (left != right) {
        middle = left+(right-left)/2+1;

        if (value_at(middle) < value) {
            left = middle;
        } else {
            right = middle-1;
        }
    }

    return left;
}

template <typename pos_t>
pos_t bin_search_min_gt(pos_t value, pos_t left, pos_t right, std::function<pos_t(pos_t)> value_at) {
    pos_t middle;

    while (left != right) {
        middle = left+(right-left)/2;

        if (value < value_at(middle)) {
            right = middle;
        } else {
            left = middle+1;
        }
    }

    return left;
}

enum exp_search_dir {LEFT,RIGHT};

template <typename pos_t, exp_search_dir search_dir>
pos_t exp_search_max_leq(pos_t value, pos_t left, pos_t right, std::function<pos_t(pos_t)> value_at) {
    pos_t cur_step_size = 1;

    if constexpr (search_dir == LEFT) {
        right -= cur_step_size;

        while (value < value_at(right)) {
            cur_step_size *= 2;

            if (right < left+cur_step_size) {
                cur_step_size = right;
                right = 0;
                break;
            }

            right -= cur_step_size;
        }

        return bin_search_max_leq<pos_t>(value,right,right+cur_step_size-1,value_at);
    } else {
        left += cur_step_size;

        while (value_at(left) < value) {
            cur_step_size *= 2;

            if (right < cur_step_size || right-cur_step_size < left) {
                cur_step_size = right-left;
                left = right;
                break;
            }

            left += cur_step_size;
        }

        return bin_search_max_leq<pos_t>(value,left-cur_step_size+1,left,value_at);
    }
}

uint64_t malloc_count_peak_memory_usage(std::ifstream& log_file) {
    std::string log_file_content;
    log_file.seekg(0,std::ios::end);
    no_init_resize(log_file_content,log_file.tellg());
    log_file.seekg(0,std::ios::beg);
    log_file.read((char*)&log_file_content[0],log_file_content.size());
    int32_t pos = 0;
    uint64_t cur_peak = 0;
    std::string str_cur_peak;
    
    while ((pos = log_file_content.find(", peak",pos)) != -1) {
        while (!('0' <= log_file_content[pos] && log_file_content[pos] <= '9')) {
            pos++;
        }

        while (('0' <= log_file_content[pos] && log_file_content[pos] <= '9') || log_file_content[pos] == '.') {
            if (log_file_content[pos] != '.') {
                str_cur_peak.push_back(log_file_content[pos]);
            }
            
            pos++;
        }

        cur_peak = std::max(cur_peak,(uint64_t)stol(str_cur_peak));
        str_cur_peak.clear();
    }

    return cur_peak;
}

template <auto start, auto end, auto inc, class T>
constexpr void for_constexpr(T&& f) {
    if constexpr (start < end) {
        f(std::integral_constant<decltype(start),start>());
        for_constexpr<start+inc,end,inc>(f);
    }
}

template<bool B, typename T>
struct constexpr_case {
    static constexpr bool value = B;
    using type = T;
};

template <bool B, typename TrueF, typename FalseF>
struct eval_if {
    using type = typename TrueF::type;
};

template <typename TrueF, typename FalseF>
struct eval_if<false, TrueF, FalseF> {
    using type = typename FalseF::type;
};    

template <bool B, typename T, typename F>
using eval_if_t = typename eval_if<B, T, F>::type;

template<typename Head, typename... Tail>
struct constexpr_switch {
    using type = eval_if_t<Head::value, Head, constexpr_switch<Tail...>>;
};

template <typename T>
struct constexpr_switch<T> {
    using type = T;
};

template <bool B, typename T>
struct constexpr_switch<constexpr_case<B, T>> {
    static_assert(B, "!");
    using type = T;
};

template<typename Head, typename... Tail>
using constexpr_switch_t = typename constexpr_switch<Head, Tail...>::type;