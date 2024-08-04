#pragma once

#define MOVE_R_BENCH 1

#include <thread>
#include <omp.h>
#include <libsais.h>
#include <libsais64.h>
#include <move_r/move_r.hpp>

void preprocess_input() {
    std::vector<std::vector<uint8_t>> contains_uchar_thr(max_num_threads,std::vector<uint8_t>(256,0));

    #pragma omp parallel for num_threads(max_num_threads)
    for (uint64_t pos=0; pos<input_size-1; pos++) {
        contains_uchar_thr[omp_get_thread_num()][char_to_uchar(input[pos])] = 1;
    }

    std::vector<uint8_t> contains_uchar(256,0);

    for (uint16_t cur_char=0; cur_char<256; cur_char++) {
        for (uint16_t cur_thread=0; cur_thread<max_num_threads; cur_thread++) {
            if (contains_uchar_thr[cur_thread][cur_char] == 1) {
                contains_uchar[cur_char] = 1;
                break;
            }
        }
    }

    contains_uchar_thr.clear();
    contains_uchar_thr.shrink_to_fit();
    uint8_t alphabet_size = 1;

    for (uint16_t cur_char=0; cur_char<256; cur_char++) {
        if (contains_uchar[cur_char] == 1) {
            alphabet_size++;
        }
    }
    
    bool contains_invalid_char = false;
    uint8_t min_valid_uchar = 3;

    for (uint8_t cur_char=0; cur_char<min_valid_uchar; cur_char++) {
        if (contains_uchar[cur_char] == 1) {
            contains_invalid_char = true;
            break;
        }
    }

    // make sure the input does not contain '\n' = 10
    if (contains_uchar[10] == 1) {
        contains_invalid_char = true;
    }

    if (contains_invalid_char) {
        if (alphabet_size > 251) {
            std::cout << "Error: the input contains more than 251 distinct characters" << std::endl;
        }

        chars_remapped = true;

        map_char.resize(256,0);
        uint16_t next_valid_char = min_valid_uchar;

        for (uint16_t cur_char=0; cur_char<256; cur_char++) {
            if (contains_uchar[cur_char] == 1) {
                map_char[cur_char] = next_valid_char;
                next_valid_char++;

                // make sure the input does not contain '\n' = 10
                if (next_valid_char == 10) {
                    next_valid_char++;
                }
            }
        }

        #pragma omp parallel for num_threads(max_num_threads)
        for (uint64_t pos=0; pos<input_size-1; pos++) {
            input[pos] = uchar_to_char(map_char[char_to_uchar(input[pos])]);
        }
    }

    contains_uchar.clear();
    contains_uchar.shrink_to_fit();
}

void update_peak_memory_usage(std::ifstream& log_file) {
    external_peak_memory_usage = std::max(external_peak_memory_usage,malloc_count_peak_memory_usage(log_file));
}

template <typename idx_t>
void build_index(idx_t& index, uint16_t num_threads);

template <typename sa_sint_t, typename idx_t>
uint64_t build_index_sa_and_bwt(idx_t& index, uint16_t num_threads);

template <typename idx_t>
void destroy_index(idx_t& index);

template <typename pos_t, typename idx_t>
pos_t count_pattern(idx_t& index, std::string& pattern);

template <typename pos_t, typename idx_t>
void locate_pattern(idx_t& index, std::string& pattern, std::vector<pos_t>& occurrences);

struct build_result {
    uint16_t num_threads;
    uint64_t time_build;
    uint64_t peak_memory_usage;
    uint64_t size_index;
};

struct query_result {
    uint64_t num_queries;
    uint64_t pattern_length;
    uint64_t num_occurrences;
    uint64_t time_query;
};

void map_string(std::string& str) {
    for (uint64_t pos=0; pos<str.size(); pos++) {
        str[pos] = uchar_to_char(map_char[char_to_uchar(str[pos])]);
    }
}

template <typename pos_t, typename idx_t>
query_result count_patterns(idx_t& index) {
    patterns_file_1.seekg(0);
    std::string header;
    std::getline(patterns_file_1,header);
    uint64_t num_queries = number_of_patterns(header);
    uint64_t pattern_length = patterns_length(header);
    uint64_t num_occurrences = 0;
    uint64_t time_query = 0;
    std::string pattern;
    no_init_resize(pattern,pattern_length);
    std::chrono::steady_clock::time_point t1,t2;

    for (uint8_t i=1; i<=2; i++) {
        for (uint64_t cur_query=0; cur_query<num_queries; cur_query++) {
            patterns_file_1.read((char*)&pattern[0],pattern_length);
            if (chars_remapped) map_string(pattern);
            count_pattern<pos_t,idx_t>(index,pattern);
        }

        patterns_file_1.seekg(0);
        std::getline(patterns_file_1,header);
    }

    for (uint64_t cur_query=0; cur_query<num_queries; cur_query++) {
        patterns_file_1.read((char*)&pattern[0],pattern_length);
        if (chars_remapped) map_string(pattern);
        t1 = now();
        num_occurrences += count_pattern<pos_t,idx_t>(index,pattern);
        t2 = now();
        time_query += time_diff_ns(t1,t2);
    }

    return query_result{num_queries,pattern_length,num_occurrences,time_query};
}

template <typename pos_t, typename idx_t>
query_result locate_patterns(idx_t& index, std::ifstream& patterns_file) {
    patterns_file.seekg(0);
    std::string header;
    std::getline(patterns_file,header);
    uint64_t num_queries = number_of_patterns(header);
    uint64_t pattern_length = patterns_length(header);
    uint64_t num_occurrences = 0;
    std::vector<pos_t> occurrences;
    uint64_t time_query = 0;
    std::string pattern;
    no_init_resize(pattern,pattern_length);
    std::chrono::steady_clock::time_point t1,t2;
    bool correct = true;

    for (uint8_t i=1; i<=2; i++) {
        for (uint64_t cur_query=0; cur_query<num_queries; cur_query++) {
            patterns_file.read((char*)&pattern[0],pattern_length);
            if (chars_remapped) map_string(pattern);
            locate_pattern<pos_t,idx_t>(index,pattern,occurrences);
            occurrences.clear();
        }

        patterns_file.seekg(0);
        std::getline(patterns_file,header);
    }

    for (uint64_t cur_query=0; cur_query<num_queries; cur_query++) {
        patterns_file.read((char*)&pattern[0],pattern_length);
        if (chars_remapped) map_string(pattern);
        t1 = now();
        locate_pattern<pos_t,idx_t>(index,pattern,occurrences);
        t2 = now();
        time_query += time_diff_ns(t1,t2);
        num_occurrences += occurrences.size();

        if (check_correctness) {
            for (pos_t occurrence : occurrences) {
                for (pos_t pos=0; pos<pattern_length; pos++) {
                    if (input[occurrence+pos] != pattern[pos]) {
                        correct = false;
                        break;
                    }
                }

                if (!correct) break;
            }

            if (!correct || occurrences.size() != count_pattern<pos_t,idx_t>(index,pattern)) {
                correct = false;
                break;
            }
        }

        occurrences.clear();
    }

    if (check_correctness) {
        if (correct) std::cout << " (no wrong occurrences)";
        else std::cout << " (wrong occurrences)";
    }

    return query_result{num_queries,pattern_length,num_occurrences,time_query};
}

void write_measurement_data(
    bool bench_count,
    bool bench_locate,
    std::string index_log_name,
    uint64_t size_index,
    build_result result_build = build_result{},
    query_result result_count = query_result{},
    query_result result_locate_1 = query_result{},
    query_result result_locate_2 = query_result{}
) {
    mf << "RESULT"
        << " type=comparison_build"
        << " implementation=" << index_log_name
        << " text=" << name_text_file
        << " text_length=" << input_size
        << " num_threads=" << result_build.num_threads
        << " time_build=" << result_build.time_build
        << " peak_memory_usage=" << result_build.peak_memory_usage
        << " size_index=" << size_index
        << std::endl;

    if (bench_count) {
        mf << "RESULT"
            << " type=comparison_count"
            << " implementation=" << index_log_name
            << " text=" << name_text_file
            << " text_length=" << input_size
            << " num_queries=" << result_count.num_queries
            << " pattern_length=" << result_count.pattern_length
            << " num_occurrences=" << result_count.num_occurrences
            << " time_query=" << result_count.time_query
            << " size_index=" << size_index
            << std::endl;
    }

    if (bench_locate) {
        std::vector<query_result> locate_results = {result_locate_1,result_locate_2};

        for (query_result& res : locate_results) {
            mf << "RESULT"
                << " type=comparison_locate"
                << " implementation=" << index_log_name
                << " text=" << name_text_file
                << " text_length=" << input_size
                << " num_queries=" << res.num_queries
                << " pattern_length=" << res.pattern_length
                << " num_occurrences=" << result_count.num_occurrences
                << " time_query=" << res.time_query
                << " size_index=" << result_build.size_index
                << std::endl;
        }
    }
}

template <typename pos_t, typename idx_t, bool bench_count, bool bench_locate>
void bench_index(std::string index_name, std::string index_log_name) {
    idx_t index;
    std::chrono::steady_clock::time_point t1,t2;
    uint64_t m1,m2;
    build_result result_build;
    query_result result_count,result_locate_1,result_locate_2;
    external_peak_memory_usage = 0;

    std::cout << "############## benchmarking " << index_name << " ##############" << std::endl << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(1));
    std::cout << "building " << index_name << std::flush;

    malloc_count_reset_peak();
    m1 = malloc_count_current();
    t1 = now();
    build_index<idx_t>(index,1);
    t2 = now();
    m2 = malloc_count_current();
    
    result_build = build_result{1,time_diff_ns(t1,t2),std::max({malloc_count_peak()-m1,external_peak_memory_usage,m2-m1}),m2-m1};

    std::cout << std::endl;
    std::cout << "build time: " << format_time(result_build.time_build) << std::endl;
    std::cout << "peak memory usage: " << format_size(result_build.peak_memory_usage) << std::endl;
    std::cout << "index size: " << format_size(result_build.size_index) << std::endl << std::endl;

    std::this_thread::sleep_for(std::chrono::seconds(1));
    
    if constexpr (bench_count) {
        std::cout << "counting the first set of patterns" << std::flush;
        result_count = count_patterns<pos_t>(index);
        if (!check_correctness) std::cout << ": " << format_query_throughput(result_count.num_queries,result_count.time_query);
        std::cout << std::endl << "total number of occurrences: " << result_count.num_occurrences << std::endl;

        std::this_thread::sleep_for(std::chrono::seconds(1));
    }

    if constexpr (bench_locate) {
        std::cout << "locating the first set of patterns" << std::flush;
        result_locate_1 = locate_patterns<pos_t,idx_t>(index,patterns_file_1);
        if (!check_correctness) std::cout << ": " << format_query_throughput(result_locate_1.num_queries,result_locate_1.time_query);
        std::cout << std::endl << "total number of occurrences: " << result_locate_1.num_occurrences << std::endl;

        std::this_thread::sleep_for(std::chrono::seconds(1));
        
        std::cout << "locating the second set of patterns" << std::flush;
        result_locate_2 = locate_patterns<pos_t,idx_t>(index,patterns_file_2);
        if (!check_correctness) std::cout << ": " << format_query_throughput(result_locate_2.num_queries,result_locate_2.time_query);
        std::cout << std::endl << "total number of occurrences: " << result_locate_2.num_occurrences << std::endl;
    }

    if constexpr (bench_count || bench_locate) std::cout << std::endl;

    if (mf.is_open()) {
        write_measurement_data(
            bench_count,
            bench_locate,
            index_log_name,
            result_build.size_index,
            result_build,
            result_count,
            result_locate_1,
            result_locate_2
        );
    }
}

template <typename pos_t, typename sa_sint_t>
void bench_a() {
    std::cout << std::endl << "building the suffix array" << std::flush;
    std::vector<sa_sint_t>& SA = get_sa<sa_sint_t>();
    no_init_resize(SA,input_size);

    if constexpr (std::is_same_v<sa_sint_t,int32_t>) {
        if (max_num_threads == 1) {
            libsais((uint8_t*)&input[0],&SA[0],input_size,0,NULL);
        } else {
            libsais_omp((uint8_t*)&input[0],&SA[0],input_size,0,NULL,max_num_threads);
        }
    } else {
        if (max_num_threads == 1) {
            libsais64((uint8_t*)&input[0],&SA[0],input_size,0,NULL);
        } else {
            libsais64_omp((uint8_t*)&input[0],&SA[0],input_size,0,NULL,max_num_threads);
        }
    }

    std::cout << std::endl << "building the BWT" << std::flush;
    no_init_resize(BWT,input_size);

    #pragma omp parallel for num_threads(max_num_threads)
    for (uint64_t pos=0; pos<input_size; pos++) {
        BWT[pos] = input[SA[pos] == 0 ? input_size-1 : SA[pos]-1];
    }

    std::cout << std::endl << std::endl;

    move_r<_locate_move,char,pos_t> index;
    uint64_t m1,m2,size_index;
    std::chrono::steady_clock::time_point t1,t2;
    query_result result_count,result_locate_1,result_locate_2;

    for (uint16_t a=2; a<=8192; a*=2) {
        std::cout << "############## benchmarking move-r with a = " << std::to_string(a) << " ##############" << std::endl;

        std::this_thread::sleep_for(std::chrono::seconds(1));

        std::cout << std::endl << "building move_r using " << format_threads(max_num_threads) << std::flush;
        index = move_r<_locate_move,char,pos_t>();
        malloc_count_reset_peak();
        m1 = malloc_count_current();
        t1 = now();
        index = move_r<_locate_move,char,pos_t>(get_sa<sa_sint_t>(),BWT,{.num_threads=max_num_threads,.a=a});
        t2 = now();
        m2 = malloc_count_current();
        size_index = m2-m1;

        std::cout << std::endl;
        std::cout << "build time: " << format_time(time_diff_ns(t1,t2)) << std::endl;
        std::cout << "peak memory usage: " << format_size(std::max({malloc_count_peak()-m1,m2-m1})) << std::endl;
        std::cout << "index size: " << format_size(size_index) << std::endl;

        std::this_thread::sleep_for(std::chrono::seconds(1));
        
        std::cout << std::endl << "counting the first set of patterns" << std::flush;
        result_count = count_patterns<pos_t>(index);
        std::cout << ": " << format_query_throughput(result_count.num_queries,result_count.time_query);
        std::cout << std::endl << "total number of occurrences: " << result_count.num_occurrences << std::endl;

        std::this_thread::sleep_for(std::chrono::seconds(1));

        std::cout << "locating the first set of patterns" << std::flush;
        result_locate_1 = locate_patterns<pos_t,move_r<_locate_move,char,pos_t>>(index,patterns_file_1);
        std::cout << ": " << format_query_throughput(result_locate_1.num_queries,result_locate_1.time_query);
        std::cout << std::endl << "total number of occurrences: " << result_locate_1.num_occurrences << std::endl;

        std::this_thread::sleep_for(std::chrono::seconds(1));
        
        std::cout << "locating the second set of patterns" << std::flush;
        result_locate_2 = locate_patterns<pos_t,move_r<_locate_move,char,pos_t>>(index,patterns_file_2);
        std::cout << ": " << format_query_throughput(result_locate_2.num_queries,result_locate_2.time_query);
        std::cout << std::endl << "total number of occurrences: " << result_locate_2.num_occurrences << std::endl;

        std::cout << std::endl;

        if (mf.is_open()) {
            mf << "RESULT"
                << " type=comparison_a_count"
                << " text=" << name_text_file
                << " text_length=" << input_size
                << " a=" << std::to_string(a)
                << " num_queries=" << result_count.num_queries
                << " pattern_length=" << result_count.pattern_length
                << " num_occurrences=" << result_count.num_occurrences
                << " time_query=" << result_count.time_query
                << " size_index=" << size_index
                << std::endl;

            std::vector<query_result> locate_results = {result_locate_1,result_locate_2};

            for (query_result& res : locate_results) {
                mf << "RESULT"
                    << " type=comparison_a_locate"
                    << " text=" << name_text_file
                    << " text_length=" << input_size
                    << " a=" << std::to_string(a)
                    << " num_queries=" << res.num_queries
                    << " pattern_length=" << res.pattern_length
                    << " num_occurrences=" << res.num_occurrences
                    << " time_query=" << res.time_query
                    << " size_index=" << size_index
                    << std::endl;
            }
        }
    }
}