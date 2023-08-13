#include <iostream>
#include <filesystem>
#include <string_view>

#include <move_r/move_r.hpp> // move-r
#include <rindex_types.hpp> // rcomp
#include "../external/r-index-prezza/internal/r_index.hpp" // r-index
#include "../external/r-index-mun/internal/r_index.hpp" // r-index-bigbwt
#include <OnlineRindex.hpp> // OnlineRlbwt
#include <DynRleForRlbwt.hpp>
#include <DynSuccForRindex.hpp>
#include <BitsUtil.cpp>
#include <dynamic/dynamic.hpp> // rle_bwt
#include <r_index_f.hpp> // r-index-f

// rcomp
using rcomp_lfig = rcomp::rindex_types::lfig_naive<7>::type;
using rcomp_glfig_16 = rcomp::rindex_types::glfig_serialized<16>::type;

// rle_bwt
using rle_bwt = dyn::rle_bwt;

// OnlineRlbwt
using BTreeNodeT = itmmti::BTreeNode<16>; // BTree arity = {16, 32, 64, 128}
using BtmNodeMT = itmmti::BtmNodeM_StepCode<BTreeNodeT, 32>; // BtmNode arity in {16, 32, 64, 128}.
using BtmMInfoT = itmmti::BtmMInfo_BlockVec<BtmNodeMT, 512>; // Each block has 512 btmNodeM.
using BtmNodeST = itmmti::BtmNodeS<BTreeNodeT, uint32_t, 8>; // CharT = uint32_t. BtmNode arity = {4, 8, 16, 32, 64, 128}.
using BtmSInfoT = itmmti::BtmSInfo_BlockVec<BtmNodeST, 1024>; // Each block has 1024 btmNodeS.
using DynRleT = itmmti::DynRleForRlbwt<itmmti::WBitsBlockVec<1024>, itmmti::Samples_WBitsBlockVec<1024>, BtmMInfoT, BtmSInfoT>;
using BtmNodeInSucc = itmmti::BtmNodeForPSumWithVal<16>; // BtmNode arity = {16, 32, 64, 128}.
using DynSuccT = itmmti::DynSuccForRindex<BTreeNodeT, BtmNodeInSucc>;
using OnlineRlbwt = itmmti::OnlineRlbwtIndex<DynRleT, DynSuccT>;

bool chars_remapped = false;
std::vector<uint8_t> map_char;
std::vector<uint8_t> unmap_char;
uint16_t max_num_threads;
std::string input;
uint64_t input_size;
std::vector<int32_t> SA_32;
std::vector<int64_t> SA_64;
std::string BWT;
std::string path_mf;
bool check_correctness = false;
bool bench_a_bigbwt = false;
uint64_t external_peak_memory_usage;
std::ofstream mf;
std::ifstream input_file;
std::ifstream patterns_file_1;
std::ifstream patterns_file_2;
std::string path_input_file;
std::string path_patterns_file_1;
std::string path_patterns_file_2;
std::string name_text_file;
int ptr = 1;

template <typename sa_sint_t>
constexpr std::vector<sa_sint_t>& get_sa() {
    if constexpr (std::is_same<sa_sint_t,int32_t>::value) {
        return SA_32;
    } else {
        return SA_64;
    }
}

void help(std::string msg) {
    if (msg != "") std::cout << msg << std::endl;
    std::cout << "move-r-bench: benchmarks construction- and query performance of move-r, move-r-bigbwt, r-index, r-index-bigbwt," << std::endl;
    std::cout << "              r-index-f, rcomp-lfig, rcomp-glfig, OnlineRLBWT and rle_bwt; has to be executed from the base folder." << std::endl << std::endl;
    std::cout << "usage 1: move-r-bench [options] <input_file> <patterns_file_1> <patterns_file_2> <num_threads>" << std::endl;
    std::cout << "   -c                 check for correctnes if possible; disables the -m option; will not print" << std::endl;
    std::cout << "                      runtime data if the runtime could be affected by checking for correctness" << std::endl;
    std::cout << "   -m <m_file>        writes measurement data to m_file" << std::endl;
    std::cout << "   <input_file>       input file" << std::endl;
    std::cout << "   <patterns_file_1>  file containing patterns (pattern length ~ number of occurrences) from <input_file>" << std::endl;
    std::cout << "                      to count and locate" << std::endl;
    std::cout << "   <patterns_file_2>  file containing patterns (pattern length << number of occurrences) from <input_file>" << std::endl;
    std::cout << "                      to locate" << std::endl;
    std::cout << "   <num_threads>      maximum number of threads to use" << std::endl;
    std::cout << std::endl;
    std::cout << "usage 2: move-r-bench -sa [options] <input_file> <num_threads>" << std::endl;
    std::cout << "                   builds the suffix array and bwt once using libsais and constructs only the" << std::endl;
    std::cout << "                   static indexes from the suffix array and the bwt (move-r and r-index)" << std::endl;
    std::cout << "                   or from the output of  prefix-free parsing (r-index-bigbwt and r-index-f)." << std::endl;
    std::cout << "   -m <m_file>     writes measurement data to m_file" << std::endl;
    std::cout << "   <input_file>    input file" << std::endl;
    std::cout << "   <num_threads>   maximum number of threads to use" << std::endl;
    std::cout << std::endl;
    std::cout << "usage 3: move-r-bench -a [options] <input_file> <patterns_file_1> <patterns_file_2> <num_threads>" << std::endl;
    std::cout << "                      constructs move_r using <num_threads> threads and measures count- and locate" << std::endl;
    std::cout << "                      performance of move_r for a=2, a=4, ..., a=8192." << std::endl;
    std::cout << "   -m <m_file>        writes measurement data to m_file" << std::endl;
    std::cout << "   <input_file>       input file" << std::endl;
    std::cout << "   <patterns_file_1>  file containing patterns (pattern length ~ number of occurrences) from <input_file>" << std::endl;
    std::cout << "                      to count and locate" << std::endl;
    std::cout << "   <patterns_file_2>  file containing patterns (pattern length << number of occurrences) from <input_file>" << std::endl;
    std::cout << "                      to locate" << std::endl;
    std::cout << "   <num_threads>   maximum number of threads to use" << std::endl;
    exit(0);
}

void parse_args(char** argv, int argc, int &ptr) {
    std::string s = argv[ptr];
    ptr++;

    if (s == "-m") {
        if (ptr >= argc-1) help("error: missing parameter after -m option");
        path_mf = argv[ptr++];
        mf.open(path_mf,std::filesystem::exists(path_mf) ? std::ios::app : std::ios::out);
        if (!mf.good()) help("error: cannot open or create measurement file");
    } else if (s == "-c") {
        if (ptr >= argc-1) help("error: missing parameter after -c option");
        check_correctness = true;
    } else {
        help("error: unrecognized '" + s +  "' option");
    }
}

// ############################# benchmarking framework #############################

uint64_t peak_memory_usage(std::ifstream& log_file) {
    std::string log_file_content;
    log_file.seekg(0,std::ios::end);
    no_init_resize(log_file_content,log_file.tellg());
    log_file.seekg(0,std::ios::beg);
    log_file.read((char*)&log_file_content[0],log_file_content.size());
    int32_t pos = 0;
    uint64_t cur_peak = 0;
    std::string str_cur_peak;

    while ((pos = log_file_content.find(", peak: ",pos)) != -1) {
        pos += 8;

        while (log_file_content[pos] != ',') {
            str_cur_peak.push_back(log_file_content[pos]);
            pos++;
        }

        cur_peak = std::max(cur_peak,(uint64_t)stol(str_cur_peak));
        str_cur_peak.clear();
    }

    return cur_peak;
}

template <typename idx_t, bool alternative_build_mode>
void build_index(idx_t& index, uint16_t num_threads);

template <typename sa_sint_t, typename idx_t>
uint64_t build_index_from_sa_and_bwt(idx_t& index, uint16_t num_threads);

template <typename idx_t>
void destroy_index(idx_t& index);

template <typename uint_t, typename idx_t>
uint_t count_pattern(idx_t& index, std::string& pattern);

template <typename uint_t, typename idx_t>
void locate_pattern(idx_t& index, std::string& pattern, std::vector<uint_t>& occurrences);

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
    for (uint_t cur_pos=0; cur_pos<str.size(); cur_pos++) {
        str[cur_pos] = uchar_to_char(map_char[char_to_uchar(str[cur_pos])]);
    }
}

void unmap_string(std::string& str) {
    for (uint_t cur_pos=0; cur_pos<str.size(); cur_pos++) {
        str[cur_pos] = uchar_to_char(unmap_char[char_to_uchar(str[cur_pos])]);
    }
}

template <typename uint_t, typename idx_t>
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
            count_pattern<uint_t,idx_t>(index,pattern);
        }

        patterns_file_1.seekg(0);
        std::getline(patterns_file_1,header);
    }

    for (uint64_t cur_query=0; cur_query<num_queries; cur_query++) {
        patterns_file_1.read((char*)&pattern[0],pattern_length);
        if (chars_remapped) map_string(pattern);
        t1 = now();
        num_occurrences += count_pattern<uint_t,idx_t>(index,pattern);
        t2 = now();
        time_query += time_diff_ns(t1,t2);
    }

    return query_result{num_queries,pattern_length,num_occurrences,time_query};
}

template <typename uint_t, typename idx_t>
query_result locate_patterns(idx_t& index, std::ifstream& patterns_file) {
    patterns_file.seekg(0);
    std::string header;
    std::getline(patterns_file,header);
    uint64_t num_queries = number_of_patterns(header);
    uint64_t pattern_length = patterns_length(header);
    uint64_t num_occurrences = 0;
    std::vector<uint_t> occurrences;
    uint64_t time_query = 0;
    std::string pattern;
    no_init_resize(pattern,pattern_length);
    std::chrono::steady_clock::time_point t1,t2;
    bool correct = true;

    for (uint8_t i=1; i<=2; i++) {
        for (uint64_t cur_query=0; cur_query<num_queries; cur_query++) {
            patterns_file.read((char*)&pattern[0],pattern_length);
            if (chars_remapped) map_string(pattern);
            locate_pattern<uint_t,idx_t>(index,pattern,occurrences);
            occurrences.clear();
        }

        patterns_file.seekg(0);
        std::getline(patterns_file,header);
    }

    for (uint64_t cur_query=0; cur_query<num_queries; cur_query++) {
        patterns_file.read((char*)&pattern[0],pattern_length);
        if (chars_remapped) map_string(pattern);
        t1 = now();
        locate_pattern<uint_t,idx_t>(index,pattern,occurrences);
        t2 = now();
        time_query += time_diff_ns(t1,t2);
        num_occurrences += occurrences.size();

        if (check_correctness) {
            for (uint_t occurrence : occurrences) {
                for (uint_t pos=0; pos<pattern_length; pos++) {
                    if (input[occurrence+pos] != pattern[pos]) {
                        correct = false;
                        break;
                    }
                }

                if (!correct) break;
            }

            if (!correct || occurrences.size() != count_pattern<uint_t,idx_t>(index,pattern)) {
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

template <typename uint_t, typename idx_t, bool alternative_build_mode, bool measure_count, bool measure_locate>
void measure(std::string index_name, std::string index_log_name, uint16_t max_build_threads) {
    idx_t index;
    std::chrono::steady_clock::time_point t1,t2;
    uint64_t m1,m2;
    std::vector<build_result> results_build;
    query_result result_count,result_locate_1,result_locate_2;
    std::cout << "############## benchmarking " << index_name << " ##############" << std::endl << std::endl;

    std::this_thread::sleep_for(1s);

    for (uint16_t cur_num_threads=1; cur_num_threads<=max_build_threads; cur_num_threads*=2) {
        std::cout << "building " << index_name << " using " << format_threads(cur_num_threads) << std::flush;

        destroy_index<idx_t>(index);
        external_peak_memory_usage = 0;
        malloc_count_reset_peak();
        m1 = malloc_count_current();
        t1 = now();
        build_index<idx_t,alternative_build_mode>(index,cur_num_threads);
        t2 = now();
        m2 = malloc_count_current();
        
        build_result res{
            .num_threads = cur_num_threads,
            .time_build = time_diff_ns(t1,t2),
            .peak_memory_usage = std::max(malloc_count_peak()-m1,external_peak_memory_usage),
            .size_index = m2-m1
        };

        results_build.emplace_back(res);
        std::cout << std::endl;
        std::cout << "build time: " << format_time(res.time_build) << std::endl;
        std::cout << "peak memory usage: " << format_size(res.peak_memory_usage) << std::endl;
        std::cout << "index size: " << format_size(res.size_index) << std::endl << std::endl;
    }
    
    if constexpr (measure_count) {
        std::cout << "counting the first set of patterns" << std::flush;
        result_count = count_patterns<uint_t>(index);
        if (!check_correctness) std::cout << ": " << format_query_throughput(result_count.num_queries,result_count.time_query);
        std::cout << std::endl << "total number of occurrences: " << result_count.num_occurrences << std::endl;

        std::this_thread::sleep_for(1s);
    }

    if constexpr (measure_locate) {
        std::cout << "locating the first set of patterns" << std::flush;
        result_locate_1 = locate_patterns<uint_t,idx_t>(index,patterns_file_1);
        if (!check_correctness) std::cout << ": " << format_query_throughput(result_locate_1.num_queries,result_locate_1.time_query);
        std::cout << std::endl << "total number of occurrences: " << result_locate_1.num_occurrences << std::endl;

        std::this_thread::sleep_for(1s);
        
        std::cout << "locating the second set of patterns" << std::flush;
        result_locate_2 = locate_patterns<uint_t,idx_t>(index,patterns_file_2);
        if (!check_correctness) std::cout << ": " << format_query_throughput(result_locate_2.num_queries,result_locate_2.time_query);
        std::cout << std::endl << "total number of occurrences: " << result_locate_2.num_occurrences << std::endl;
    }

    if constexpr (measure_count || measure_locate) std::cout << std::endl;

    if (mf.is_open()) {
        uint64_t size_index = results_build.front().size_index;

        for (build_result& res : results_build) {
            mf << "RESULT"
                << " type=comparison_build"
                << " implementation=" << index_log_name
                << " text=" << name_text_file
                << " num_threads=" << res.num_threads
                << " time_build=" << res.time_build
                << " peak_memory_usage=" << res.peak_memory_usage
                << " size_index=" << res.size_index
                << std::endl;
        }

        if constexpr (measure_count) {
            mf << "RESULT"
                << " type=comparison_count"
                << " implementation=" << index_log_name
                << " text=" << name_text_file
                << " num_queries=" << result_count.num_queries
                << " pattern_length=" << result_count.pattern_length
                << " num_occurrences=" << result_count.num_occurrences
                << " time_query=" << result_count.time_query
                << " size_index=" << size_index
                << std::endl;
        }

        if constexpr (measure_locate) {
            std::vector<query_result> locate_results = {result_locate_1,result_locate_2};

            for (query_result& res : locate_results) {
                mf << "RESULT"
                    << " type=comparison_locate"
                    << " implementation=" << index_log_name
                    << " text=" << name_text_file
                    << " num_queries=" << res.num_queries
                    << " pattern_length=" << res.pattern_length
                    << " num_occurrences=" << result_count.num_occurrences
                    << " time_query=" << res.time_query
                    << " size_index=" << results_build.back().size_index
                    << std::endl;
            }
        }
    }
}

template <typename sa_sint_t, typename idx_t>
void measure_construct_from_sa_and_bwt(std::string index_name, std::string index_log_name, uint16_t max_build_threads) {
    idx_t index;
    std::chrono::steady_clock::time_point t1,t2;
    uint64_t m1,m2;
    std::vector<build_result> results_build;
    uint64_t time_build_override;
    std::cout << "############## benchmarking " << index_name << " ##############" << std::endl << std::endl;

    for (uint16_t cur_num_threads=1; cur_num_threads<=max_build_threads; cur_num_threads*=2) {
        std::cout << "building " << index_name << " using " << format_threads(cur_num_threads) << std::flush;
        destroy_index<idx_t>(index);

        std::this_thread::sleep_for(1s);

        external_peak_memory_usage = 0;
        malloc_count_reset_peak();
        m1 = malloc_count_current();
        t1 = now();
        time_build_override = build_index_from_sa_and_bwt<sa_sint_t,idx_t>(index,cur_num_threads);
        t2 = now();
        m2 = malloc_count_current();
        
        build_result res{
            .num_threads = cur_num_threads,
            .time_build = time_build_override != 0 ? time_build_override : time_diff_ns(t1,t2),
            .peak_memory_usage = std::max(malloc_count_peak()-m1,external_peak_memory_usage),
            .size_index = m2-m1
        };

        results_build.emplace_back(res);
        std::cout << std::endl;
        std::cout << "build time: " << format_time(res.time_build) << std::endl;
        std::cout << "peak memory usage: " << format_size(res.peak_memory_usage) << std::endl;
        std::cout << "index size: " << format_size(res.size_index) << std::endl << std::endl;
    }

    if (mf.is_open()) {
        for (build_result& res : results_build) {
            mf << "RESULT"
                << " type=comparison_build_from_sa_and_bwt"
                << " implementation=" << index_log_name
                << " text=" << name_text_file
                << " num_threads=" << res.num_threads
                << " time_build=" << res.time_build
                << " peak_memory_usage=" << res.peak_memory_usage
                << " size_index=" << res.size_index
                << std::endl;
        }
    }
}

template <typename uint_t, typename sa_sint_t>
void benchmark_a() {
    move_r<uint_t> index;
    uint64_t m1,m2,size_index;
    std::chrono::steady_clock::time_point t1,t2;
    query_result result_count,result_locate_1,result_locate_2;

    for (uint16_t a=2; a<=8192; a*=2) {
        std::cout << "############## benchmarking move-r with a = " << std::to_string(a) << " ##############" << std::endl;

        std::this_thread::sleep_for(1s);

        std::cout << std::endl << "building move_r using " << format_threads(max_num_threads) << std::flush;
        index = move_r<uint_t>();
        malloc_count_reset_peak();
        m1 = malloc_count_current();
        t1 = now();
        index = move_r<uint_t>(get_sa<sa_sint_t>(),BWT,full_support,max_num_threads,a);
        t2 = now();
        m2 = malloc_count_current();
        size_index = m2-m1;

        std::cout << std::endl;
        std::cout << "build time: " << format_time(time_diff_ns(t1,t2)) << std::endl;
        std::cout << "peak memory usage: " << format_size(malloc_count_peak()-m1) << std::endl;
        std::cout << "index size: " << format_size(size_index) << std::endl;

        std::this_thread::sleep_for(1s);
        
        std::cout << std::endl << "counting the first set of patterns" << std::flush;
        result_count = count_patterns<uint_t>(index);
        std::cout << ": " << format_query_throughput(result_count.num_queries,result_count.time_query);
        std::cout << std::endl << "total number of occurrences: " << result_count.num_occurrences << std::endl;

        std::this_thread::sleep_for(1s);

        std::cout << "locating the first set of patterns" << std::flush;
        result_locate_1 = locate_patterns<uint_t,move_r<uint_t>>(index,patterns_file_1);
        std::cout << ": " << format_query_throughput(result_locate_1.num_queries,result_locate_1.time_query);
        std::cout << std::endl << "total number of occurrences: " << result_locate_1.num_occurrences << std::endl;

        std::this_thread::sleep_for(1s);
        
        std::cout << "locating the second set of patterns" << std::flush;
        result_locate_2 = locate_patterns<uint_t,move_r<uint_t>>(index,patterns_file_2);
        std::cout << ": " << format_query_throughput(result_locate_2.num_queries,result_locate_2.time_query);
        std::cout << std::endl << "total number of occurrences: " << result_locate_2.num_occurrences << std::endl;

        std::cout << std::endl;

        if (mf.is_open()) {
            mf << "RESULT"
                << " type=comparison_a_count"
                << " text=" << name_text_file
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

template <typename uint_t>
void measure_all() {
    measure<uint_t,move_r<uint_t>,false,true,true>("move-r","move_r",max_num_threads);
    measure<uint_t,move_r<uint_t>,true,false,false>("move-r-bigbwt","move_r_bigbwt",max_num_threads);
    measure<uint_t,r_index_f<>,false,true,false>("r-index-f","r_index_f",1);
    measure<uint_t,rcomp_lfig,false,true,true>("rcomp-lfig","rcomp_lfig",1);
    measure<uint_t,rcomp_glfig_16,false,true,true>("rcomp-glfig","rcomp_glfig_g16",1);
    measure<uint_t,ri::r_index<>,false,true,true>("r-index","r_index",1);
    measure<uint_t,ri_mun::r_index<>,false,false,false>("r-index-bigbwt","r_index_bigbwt",1);
    measure<uint_t,OnlineRlbwt,false,true,true>("OnlineRlbwt","online_rlbwt",1);
    measure<uint_t,rle_bwt,false,true,false>("rle_bwt","rle_bwt",1);
}

void preprocess_input() {
    std::vector<std::vector<uint8_t>> contains_uchar_thr(max_num_threads,std::vector<uint8_t>(256,0));

    #pragma omp parallel for num_threads(max_num_threads)
    for (uint64_t i=0; i<input_size-1; i++) {
        contains_uchar_thr[omp_get_thread_num()][char_to_uchar(input[i])] = 1;
    }

    std::vector<uint8_t> contains_uchar(256,0);

    for (uint16_t i=0; i<256; i++) {
        for (uint16_t j=0; j<max_num_threads; j++) {
            if (contains_uchar_thr[j][i] == 1) {
                contains_uchar[i] = 1;
                break;
            }
        }
    }

    contains_uchar_thr.clear();
    contains_uchar_thr.shrink_to_fit();
    uint8_t alphabet_size = 1;

    for (uint16_t i=0; i<256; i++) {
        if (contains_uchar[i] == 1) {
            alphabet_size++;
        }
    }
    
    bool contains_invalid_char = false;
    uint8_t min_valid_char = 3;

    for (uint8_t i=0; i<min_valid_char; i++) {
        if (contains_uchar[i] == 1) {
            contains_invalid_char = true;
            break;
        }
    }

    if (contains_invalid_char) {
        if (alphabet_size > 252) {
            std::cout << "Error: the input contains more than 252 distinct characters" << std::endl;
        }

        chars_remapped = true;

        map_char.resize(256,0);
        unmap_char.resize(256,0);
        uint16_t j = min_valid_char;

        for (uint16_t i=0; i<256; i++) {
            if (contains_uchar[i] == 1) {
                map_char[i] = j;
                unmap_char[j] = i;
                j++;
            }
        }

        #pragma omp parallel for num_threads(max_num_threads)
        for (uint64_t i=0; i<input_size-1; i++) {
            input[i] = uchar_to_char(map_char[char_to_uchar(input[i])]);
        }
    }

    contains_uchar.clear();
    contains_uchar.shrink_to_fit();
}

template <bool bench_a, typename uint_t, typename sa_sint_t>
void measure_all_construct_from_sa_and_bwt() {
    std::cout << std::endl << "building the suffix array" << std::flush;
    std::vector<sa_sint_t>& SA = get_sa<sa_sint_t>();
    no_init_resize(SA,input_size);

    if constexpr (std::is_same<sa_sint_t,int32_t>::value) {
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

    if constexpr (bench_a) {
        benchmark_a<uint_t,sa_sint_t>();
    } else {
        measure_construct_from_sa_and_bwt<sa_sint_t,move_r<uint_t>>("move-r","move_r",max_num_threads);
        measure_construct_from_sa_and_bwt<sa_sint_t,r_index_f<>>("r-index-f","r_index_f",1);
        measure_construct_from_sa_and_bwt<sa_sint_t,ri::r_index<>>("r-index","r_index",1);
        measure_construct_from_sa_and_bwt<sa_sint_t,ri_mun::r_index<>>("r-index-bigbwt","r_index_bigbwt",1);
    }
}

// ############################# move-r #############################

template <>
void build_index<move_r<uint32_t>,false>(move_r<uint32_t>& index, uint16_t num_threads) {
    index = std::move(move_r<uint32_t>(input,full_support,runtime,num_threads));
}

template <>
void build_index<move_r<uint64_t>,false>(move_r<uint64_t>& index, uint16_t num_threads) {
    index = std::move(move_r<uint64_t>(input,full_support,runtime,num_threads));
}

template <>
void build_index<move_r<uint32_t>,true>(move_r<uint32_t>& index, uint16_t num_threads) {
    index = std::move(move_r<uint32_t>(input,full_support,space,num_threads));
}

template <>
void build_index<move_r<uint64_t>,true>(move_r<uint64_t>& index, uint16_t num_threads) {
    index = std::move(move_r<uint64_t>(input,full_support,space,num_threads));
}

template <>
uint64_t build_index_from_sa_and_bwt<int32_t,move_r<uint32_t>>(move_r<uint32_t>& index, uint16_t num_threads) {
    index = std::move(move_r<uint32_t>(get_sa<int32_t>(),BWT,full_support,num_threads));
    return 0;
}

template <>
uint64_t build_index_from_sa_and_bwt<int64_t,move_r<uint32_t>>(move_r<uint32_t>& index, uint16_t num_threads) {
    index = std::move(move_r<uint32_t>(get_sa<int64_t>(),BWT,full_support,num_threads));
    return 0;
}

template <>
uint64_t build_index_from_sa_and_bwt<int64_t,move_r<uint64_t>>(move_r<uint64_t>& index, uint16_t num_threads) {
    index = std::move(move_r<uint64_t>(get_sa<int64_t>(),BWT,full_support,num_threads));
    return 0;
}

template <>
void destroy_index<move_r<uint32_t>>(move_r<uint32_t>& index) {
    index = std::move(move_r<uint32_t>());
}

template <>
void destroy_index<move_r<uint64_t>>(move_r<uint64_t>& index) {
    index = std::move(move_r<uint64_t>());
}

template <>
uint32_t count_pattern<uint32_t,move_r<uint32_t>>(move_r<uint32_t>& index, std::string& pattern) {
    return index.count(pattern);
}

template <>
uint64_t count_pattern<uint64_t,move_r<uint64_t>>(move_r<uint64_t>& index, std::string& pattern) {
    return index.count(pattern);
}

template <>
void locate_pattern<uint32_t,move_r<uint32_t>>(move_r<uint32_t>& index, std::string& pattern, std::vector<uint32_t>& occurrences) {
    index.locate(pattern,occurrences);
}

template <>
void locate_pattern<uint64_t,move_r<uint64_t>>(move_r<uint64_t>& index, std::string& pattern, std::vector<uint64_t>& occurrences) {
    index.locate(pattern,occurrences);
}

// ############################# rcomp #############################

template <>
void build_index<rcomp_lfig,false>(rcomp_lfig& index, uint16_t) {
    for (uint64_t pos=1; pos<input_size; pos++) {
        index.extend(input[input_size-1-pos]);
    }
}

template <>
void build_index<rcomp_glfig_16,false>(rcomp_glfig_16& index, uint16_t) {
    for (uint64_t pos=1; pos<input_size; pos++) {
        index.extend(input[input_size-1-pos]);
    }
}

template <>
void destroy_index<rcomp_lfig>(rcomp_lfig&) {}

template <>
void destroy_index<rcomp_glfig_16>(rcomp_glfig_16&) {}

template <>
uint32_t count_pattern<uint32_t,rcomp_lfig>(rcomp_lfig& index, std::string& pattern) {
    std::string reversed_pattern(pattern.rbegin(),pattern.rend());
    return index.count(rcomp::make_range(reversed_pattern));
}

template <>
uint64_t count_pattern<uint64_t,rcomp_lfig>(rcomp_lfig& index, std::string& pattern) {
    std::string reversed_pattern(pattern.rbegin(),pattern.rend());
    return index.count(rcomp::make_range(reversed_pattern));
}

template <>
uint32_t count_pattern<uint32_t,rcomp_glfig_16>(rcomp_glfig_16& index, std::string& pattern) {
    std::string reversed_pattern(pattern.rbegin(),pattern.rend());
    return index.count(rcomp::make_range(reversed_pattern));
}

template <>
uint64_t count_pattern<uint64_t,rcomp_glfig_16>(rcomp_glfig_16& index, std::string& pattern) {
    std::string reversed_pattern(pattern.rbegin(),pattern.rend());
    return index.count(rcomp::make_range(reversed_pattern));
}

template <>
void locate_pattern<uint32_t,rcomp_lfig>(rcomp_lfig& index, std::string& pattern, std::vector<uint32_t>& occurrences) {
    std::string reversed_pattern(pattern.rbegin(),pattern.rend());
    index.locate(rcomp::make_range(reversed_pattern),[&occurrences](uint32_t occurrence){occurrences.emplace_back(input_size-1-occurrence);});
}

template <>
void locate_pattern<uint64_t,rcomp_lfig>(rcomp_lfig& index, std::string& pattern, std::vector<uint64_t>& occurrences) {
    std::string reversed_pattern(pattern.rbegin(),pattern.rend());
    index.locate(rcomp::make_range(reversed_pattern),[&occurrences](uint64_t occurrence){occurrences.emplace_back(input_size-1-occurrence);});
}

template <>
void locate_pattern<uint32_t,rcomp_glfig_16>(rcomp_glfig_16& index, std::string& pattern, std::vector<uint32_t>& occurrences) {
    std::string reversed_pattern(pattern.rbegin(),pattern.rend());
    index.locate(rcomp::make_range(reversed_pattern),[&occurrences](uint32_t occurrence){occurrences.emplace_back(input_size-1-occurrence);});
}

template <>
void locate_pattern<uint64_t,rcomp_glfig_16>(rcomp_glfig_16& index, std::string& pattern, std::vector<uint64_t>& occurrences) {
    std::string reversed_pattern(pattern.rbegin(),pattern.rend());
    index.locate(rcomp::make_range(reversed_pattern),[&occurrences](uint64_t occurrence){occurrences.emplace_back(input_size-1-occurrence);});
}

// ############################# r-index #############################

template <>
void build_index<ri::r_index<>,false>(ri::r_index<>& index, uint16_t) {
    std::streambuf* cout_rfbuf = cout.rdbuf();
    std::cout.rdbuf(NULL);
    index = ri::r_index<>(input,true);
    std::cout.rdbuf(cout_rfbuf);
}

template <>
uint64_t build_index_from_sa_and_bwt<int32_t,ri::r_index<>>(ri::r_index<>& index, uint16_t) {
    index = std::move(ri::r_index<>(get_sa<int32_t>(),BWT));
    return 0;
}

template <>
uint64_t build_index_from_sa_and_bwt<int64_t,ri::r_index<>>(ri::r_index<>& index, uint16_t) {
    index = std::move(ri::r_index<>(get_sa<int64_t>(),BWT));
    return 0;
}

template <>
void destroy_index<ri::r_index<>>(ri::r_index<>& index) {
    index = std::move(ri::r_index<>());
}

template <>
uint32_t count_pattern<uint32_t,ri::r_index<>>(ri::r_index<>& index, std::string& pattern) {
    auto range = index.count(pattern);
    return range.second - range.first + 1;
}

template <>
uint64_t count_pattern<uint64_t,ri::r_index<>>(ri::r_index<>& index, std::string& pattern) {
    auto range = index.count(pattern);
    return range.second - range.first + 1;
}

template <>
void locate_pattern<uint32_t,ri::r_index<>>(ri::r_index<>& index, std::string& pattern, std::vector<uint32_t>& occurrences) {
    index.locate_all<uint32_t>(pattern,occurrences);
}

template <>
void locate_pattern<uint64_t,ri::r_index<>>(ri::r_index<>& index, std::string& pattern, std::vector<uint64_t>& occurrences) {
    index.locate_all<uint64_t>(pattern,occurrences);
}

// ############################# r-index-bigbwt #############################

template <>
void build_index<ri_mun::r_index<>,false>(ri_mun::r_index<>& index, uint16_t num_threads) {
    std::streambuf* cout_rfbuf = cout.rdbuf();
    std::cout.rdbuf(NULL);
    std::string name_text_file = "r-index-" + random_alphanumeric_string(10);
    std::ofstream text_file(name_text_file);
    write_to_file(text_file,input.c_str(),input_size-1);
    text_file.close();
    index = ri_mun::r_index<>(name_text_file,"bigbwt",num_threads);
    std::filesystem::remove(name_text_file);
    std::ifstream log_file(name_text_file + ".log");
    external_peak_memory_usage = peak_memory_usage(log_file);
    log_file.close();
    std::filesystem::remove(name_text_file + ".log");
    system("rm -f nul");
    std::cout.rdbuf(cout_rfbuf);
}

uint64_t build_r_index_bigbwt_from_sa_and_bwt(ri_mun::r_index<>& index) {
    std::streambuf* cout_rfbuf = cout.rdbuf();
    std::cout.rdbuf(NULL);
    std::string name_text_file = "r-index-" + random_alphanumeric_string(10);
    std::ofstream text_file(name_text_file);
    write_to_file(text_file,input.c_str(),input_size-1);
    text_file.close();
    
    uint64_t time_build_override;
    index = ri_mun::r_index<>(name_text_file,"bigbwt",1,&time_build_override);
    std::filesystem::remove(name_text_file);
    std::filesystem::remove(name_text_file + ".log");
    system("rm -f nul");
    std::cout.rdbuf(cout_rfbuf);

    return time_build_override;
}

template <>
uint64_t build_index_from_sa_and_bwt<int32_t,ri_mun::r_index<>>(ri_mun::r_index<>& index, uint16_t num_threads) {
    return build_r_index_bigbwt_from_sa_and_bwt(index);
}

template <>
uint64_t build_index_from_sa_and_bwt<int64_t,ri_mun::r_index<>>(ri_mun::r_index<>& index, uint16_t num_threads) {
    return build_r_index_bigbwt_from_sa_and_bwt(index);
}

template <>
void destroy_index<ri_mun::r_index<>>(ri_mun::r_index<>& index) {
    index = std::move(ri_mun::r_index<>());
}

template <>
uint32_t count_pattern<uint32_t,ri_mun::r_index<>>(ri_mun::r_index<>& index, std::string& pattern) {
    auto range = index.count(pattern);
    return range.second - range.first + 1;
}

template <>
uint64_t count_pattern<uint64_t,ri_mun::r_index<>>(ri_mun::r_index<>& index, std::string& pattern) {
    auto range = index.count(pattern);
    return range.second - range.first + 1;
}

template <>
void locate_pattern<uint32_t,ri_mun::r_index<>>(ri_mun::r_index<>& index, std::string& pattern, std::vector<uint32_t>& occurrences) {
    index.locate_all<uint32_t>(pattern,occurrences);
}

template <>
void locate_pattern<uint64_t,ri_mun::r_index<>>(ri_mun::r_index<>& index, std::string& pattern, std::vector<uint64_t>& occurrences) {
    index.locate_all<uint64_t>(pattern,occurrences);
}

// ############################# OnlineRlbwt #############################

template <>
void build_index<OnlineRlbwt,false>(OnlineRlbwt& index, uint16_t) {
    for (uint64_t pos=0; pos<input_size-1; pos++) {
        index.extend(input[pos]);
    }
}

template <>
void destroy_index<OnlineRlbwt>(OnlineRlbwt&) {}

template <>
uint32_t count_pattern<uint32_t,OnlineRlbwt>(OnlineRlbwt& index, std::string& pattern) {
    OnlineRlbwt::PatTracker tracker = index.getInitialPatTracker();

    for (uint32_t pos=0; pos<pattern.size(); pos++) {
        index.lfMap(tracker,pattern[pos]);
    }

    return index.getNumOcc(tracker);
}

template <>
uint64_t count_pattern<uint64_t,OnlineRlbwt>(OnlineRlbwt& index, std::string& pattern) {
    OnlineRlbwt::PatTracker tracker = index.getInitialPatTracker();

    for (uint64_t pos=0; pos<pattern.size(); pos++) {
        index.lfMap(tracker,pattern[pos]);
    }

    return index.getNumOcc(tracker);
}

template <>
void locate_pattern<uint32_t,OnlineRlbwt>(OnlineRlbwt& index, std::string& pattern, std::vector<uint32_t>& occurrences) {
       OnlineRlbwt::PatTracker tracker = index.getInitialPatTracker();

    for (uint32_t pos=0; pos<pattern.size(); pos++) {
        index.lfMap(tracker,pattern[pos]);
    }

    uint64_t numOcc = index.getNumOcc(tracker);
    uint64_t curPos = index.calcFstOcc(tracker);
    occurrences.emplace_back(curPos-pattern.size());

    for (uint32_t pos=1; pos<numOcc; pos++) {
        curPos = index.calcNextPos(curPos);
          occurrences.emplace_back(curPos-pattern.size());
    }
}

template <>
void locate_pattern<uint64_t,OnlineRlbwt>(OnlineRlbwt& index, std::string& pattern, std::vector<uint64_t>& occurrences) {
       OnlineRlbwt::PatTracker tracker = index.getInitialPatTracker();

    for (uint64_t pos=0; pos<pattern.size(); pos++) {
        index.lfMap(tracker,pattern[pos]);
    }

    uint64_t numOcc = index.getNumOcc(tracker);
    uint64_t curPos = index.calcFstOcc(tracker);
    occurrences.emplace_back(curPos-pattern.size());

    for (uint64_t pos=1; pos<numOcc; pos++) {
        curPos = index.calcNextPos(curPos);
          occurrences.emplace_back(curPos-pattern.size());
    }
}

// ############################# rle_bwt #############################

template <>
void build_index<rle_bwt,false>(rle_bwt& index, uint16_t) {
    for (uint64_t pos=0; pos<input_size-1; pos++) {
        index.extend(input[pos]);
    }
}

template <>
void destroy_index<rle_bwt>(rle_bwt&) {}

template <>
uint32_t count_pattern<uint32_t,rle_bwt>(rle_bwt& index, std::string& pattern) {
    auto range = index.count(pattern);
    return range.second - range.first;
}

template <>
uint64_t count_pattern<uint64_t,rle_bwt>(rle_bwt& index, std::string& pattern) {
    auto range = index.count(pattern);
    return range.second - range.first;
}

// ############################# r-index-f #############################

template <>
void build_index<r_index_f<>,false>(r_index_f<>& index, uint16_t) {
    std::streambuf* cout_rfbuf = cout.rdbuf();
    std::cout.rdbuf(NULL);
    std::string name_text_file = "r-index-f-" + random_alphanumeric_string(10);
    std::ofstream text_file(name_text_file);
    write_to_file(text_file,input.c_str(),input_size-1);
    text_file.close();
    std::string cmd_newscan = "build/external/Big-BWT/newscanNT.x " + name_text_file + " >nul 2>nul";
    system(cmd_newscan.c_str());
    std::string cmd_pfp_thresholds = "build/external/pfp-thresholds/pfp-thresholds " +
    name_text_file + " -r > " + name_text_file + ".log 2>" + name_text_file + ".log";
    system(cmd_pfp_thresholds.c_str());
    index = std::move(r_index_f<>(name_text_file));
    std::ifstream log_file(name_text_file + ".log");
    external_peak_memory_usage = peak_memory_usage(log_file);
    log_file.close();
    std::string cmd_rm = "rm -f " + name_text_file + "* >nul 2>nul";
    system(cmd_rm.c_str());
    system("rm -f nul");
    std::cout.rdbuf(cout_rfbuf);
}

uint64_t build_r_index_f_from_sa_and_bwt(r_index_f<>& index) {
    std::streambuf* cout_rfbuf = cout.rdbuf();
    std::cout.rdbuf(NULL);
    std::string name_text_file = "r-index-f-" + random_alphanumeric_string(10);
    std::ofstream text_file(name_text_file);
    write_to_file(text_file,input.c_str(),input_size-1);
    text_file.close();
    std::string cmd_newscan = "build/external/Big-BWT/newscanNT.x " + name_text_file + " >nul 2>nul";
    system(cmd_newscan.c_str());
    std::string cmd_pfp_thresholds = "build/external/pfp-thresholds/pfp-thresholds " +
    name_text_file + " -r > " + name_text_file + ".log 2>" + name_text_file + ".log";
    system(cmd_pfp_thresholds.c_str());

    auto t1 = now();
    index = std::move(r_index_f<>(name_text_file));
    auto t2 = now();
    
    std::string cmd_rm = "rm -f " + name_text_file + "* >nul 2>nul";
    system(cmd_rm.c_str());
    system("rm -f nul");
    std::cout.rdbuf(cout_rfbuf);
    
    return time_diff_ns(t1,t2);
}

template <>
uint64_t build_index_from_sa_and_bwt<int32_t,r_index_f<>>(r_index_f<>& index, uint16_t num_threads) {
    return build_r_index_f_from_sa_and_bwt(index);
}

template <>
uint64_t build_index_from_sa_and_bwt<int64_t,r_index_f<>>(r_index_f<>& index, uint16_t) {
    return build_r_index_f_from_sa_and_bwt(index);
}

template <>
void destroy_index<r_index_f<>>(r_index_f<>& index) {
    index = std::move(r_index_f<>());
}

template <>
uint32_t count_pattern<uint32_t,r_index_f<>>(r_index_f<>& index, std::string& pattern) {
    return index.count(pattern);
}

template <>
uint64_t count_pattern<uint64_t,r_index_f<>>(r_index_f<>& index, std::string& pattern) {
    return index.count(pattern);
}

// ############################# benchmarking framework #############################

template <bool bench_a>
int main_construct_from_sa_and_bwt(int argc, char** argv) {
    if constexpr (bench_a) {
        if (argc == 8) {
            if (std::string(argv[2]) != "-m") help("");

            path_mf = argv[3];
            mf.open(path_mf,std::filesystem::exists(path_mf) ? std::ios::app : std::ios::out);
            if (!mf.good()) help("error: cannot open or create measurement file");

            path_input_file = argv[4];
            path_patterns_file_1 = argv[5];
            path_patterns_file_2 = argv[6];
            max_num_threads = atoi(argv[7]);
        } else if (argc == 6) {
            path_input_file = argv[2];
            path_patterns_file_1 = argv[3];
            path_patterns_file_2 = argv[4];
            max_num_threads = atoi(argv[5]);
        } else {
            help("");
        }

        patterns_file_1.open(path_patterns_file_1);
        patterns_file_2.open(path_patterns_file_2);
        if (!patterns_file_1.good()) help("error: invalid input, could read <patterns_file_1>");
        if (!patterns_file_2.good()) help("error: invalid input, could read <patterns_file_2>");
    } else {
        if (argc == 6) {
            if (std::string(argv[2]) != "-m") help("");

            path_mf = argv[3];
            path_input_file = argv[4];
            max_num_threads = atoi(argv[5]);

            mf.open(path_mf,std::filesystem::exists(path_mf) ? std::ios::app : std::ios::out);
            if (!mf.good()) help("error: cannot open or create measurement file");
        } else if (argc == 4) {
            path_input_file = argv[2];
            max_num_threads = atoi(argv[3]);
        } else {
            help("");
        }
    }
    
    input_file.open(path_input_file);
    if (!input_file.good()) help("error: invalid input, could not read <input_file>");
    if (max_num_threads == 0 || max_num_threads > omp_get_max_threads()) help("error: invalid number of threads");

    system("chmod +x build/external/Big-BWT/*");
    system("chmod +x build/external/pfp-thresholds/*");
    std::cout << std::setprecision(4);
    name_text_file = path_input_file.substr(path_input_file.find_last_of("/\\") + 1);
    std::cout << "benchmarking " << path_input_file << std::flush;

    input_file.seekg(0,std::ios::end);
    input_size = input_file.tellg()+(std::streamsize)+1;
    input_file.seekg(0,std::ios::beg);
    no_init_resize(input,input_size);
    read_from_file(input_file,input.c_str(),input_size-1);
    input[input_size-1] = 1;
    input_file.close();
    preprocess_input();
    
    std::cout << " (" << format_size(input_size-1) << ") using up to " << format_threads(max_num_threads) << std::endl;

    if (input_size <= UINT_MAX) {
        if (input_size <= INT_MAX) {
            measure_all_construct_from_sa_and_bwt<bench_a,uint32_t,int32_t>();
        } else {
            measure_all_construct_from_sa_and_bwt<bench_a,uint32_t,int64_t>();
        }
    } else {
        measure_all_construct_from_sa_and_bwt<bench_a,uint64_t,int64_t>();
    }

    if (mf.is_open()) mf.close();
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 1) help("");
    if (std::string(argv[1]) == "-sa") return main_construct_from_sa_and_bwt<false>(argc,argv);
    if (std::string(argv[1]) == "-a") return main_construct_from_sa_and_bwt<true>(argc,argv);
    if (argc < 5) help("");
    while (ptr < argc - 4) parse_args(argv, argc, ptr);

    path_input_file = argv[ptr];
    path_patterns_file_1 = argv[ptr+1];
    path_patterns_file_2 = argv[ptr+2];
    max_num_threads = atoi(argv[ptr+3]);

    input_file.open(path_input_file);
    patterns_file_1.open(path_patterns_file_1);
    patterns_file_2.open(path_patterns_file_2);

    if (!input_file.good()) help("error: invalid input, could not read <input_file>");
    if (!patterns_file_1.good()) help("error: invalid input, could not read <patterns_file_1>");
    if (!patterns_file_2.good()) help("error: invalid input, could not read <patterns_file_2>");
    if (mf.is_open() && check_correctness) help("error: cannot output measurement data when checking for correctness");
    if (max_num_threads == 0 || max_num_threads > omp_get_max_threads()) help("error: invalid number of threads");

    system("chmod +x build/external/Big-BWT/*");
    system("chmod +x build/external/pfp-thresholds/*");
    std::cout << std::setprecision(4);
    name_text_file = path_input_file.substr(path_input_file.find_last_of("/\\") + 1);
    std::cout << "benchmarking " << path_input_file << std::flush;

    input_file.seekg(0,std::ios::end);
    input_size = input_file.tellg()+(std::streamsize)+1;
    input_file.seekg(0,std::ios::beg);
    input.reserve(input_size);
    no_init_resize(input,input_size-1);
    read_from_file(input_file,input.c_str(),input_size-1);
    input_file.close();
    preprocess_input();

	std::cout << " (" << format_size(input_size-1) << ") using up to " << format_threads(max_num_threads) << std::endl;
    if (check_correctness) std::cout << "correctnes will be checked if possible" << std::endl;
    std::cout << std::endl;

    if (input_size <= UINT_MAX) {
        measure_all<uint32_t>();
    } else {
        measure_all<uint64_t>();
    }

    patterns_file_1.close();
    patterns_file_2.close();

    if (mf.is_open()) mf.close();
    return 0;
}