#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>
#include <string_view>

bool chars_remapped = false;
std::vector<uint8_t> map_char;
uint16_t max_num_threads;
std::string input;
uint64_t input_size;
std::vector<int32_t> SA_32;
std::vector<int64_t> SA_64;
std::string BWT;
std::string path_mf;
bool check_correctness = false;
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

#include "framework.hpp"
#include "indexes.hpp"

void help(std::string msg) {
    if (msg != "") std::cout << msg << std::endl;
    std::cout << "move-r-bench: benchmarks construction- and query performance of move-r, block-rlbwt-2, block-rlbwt-v," << std::endl;
    std::cout << "              block-rlbwt-r, r-index, r-index-f, rcomp-glfig and online-rlbwt;" << std::endl;
    std::cout << "              has to be executed from the base folder." << std::endl << std::endl;
    std::cout << "usage 1: move-r-bench [options] <input_file> <patterns_file_1> <patterns_file_2>" << std::endl;
    std::cout << "   -c                 check for correctnes if possible; disables the -m option; will not print" << std::endl;
    std::cout << "                      runtime data if the runtime could be affected by checking for correctness" << std::endl;
    std::cout << "   -m <m_file>        writes bench_indexment data to m_file" << std::endl;
    std::cout << "   <input_file>       input file" << std::endl;
    std::cout << "   <patterns_file_1>  file containing patterns (pattern length ~ number of occurrences) from <input_file>" << std::endl;
    std::cout << "                      to count and locate" << std::endl;
    std::cout << "   <patterns_file_2>  file containing patterns (pattern length << number of occurrences) from <input_file>" << std::endl;
    std::cout << "                      to locate" << std::endl;
    std::cout << std::endl;
    std::cout << "usage 2: move-r-bench -a [options] <input_file> <patterns_file_1> <patterns_file_2> <num_threads>" << std::endl;
    std::cout << "                      constructs move_r using <num_threads> threads and bench_indexs count- and locate" << std::endl;
    std::cout << "                      performance of move_r for a=2, a=4, ..., a=8192." << std::endl;
    std::cout << "   -m <m_file>        writes bench_indexment data to m_file" << std::endl;
    std::cout << "   <input_file>       input file" << std::endl;
    std::cout << "   <patterns_file_1>  file containing patterns (pattern length ~ number of occurrences) from <input_file>" << std::endl;
    std::cout << "                      to count and locate" << std::endl;
    std::cout << "   <patterns_file_2>  file containing patterns (pattern length << number of occurrences) from <input_file>" << std::endl;
    std::cout << "                      to locate" << std::endl;
    std::cout << "   <num_threads>   maximum number of threads to use" << std::endl;
    exit(0);
}

template <typename pos_t>
void bench_indexes() {
    bench_index<pos_t,move_r<_phi,char,pos_t>,true,true>("move-r","move_r");
    
    /*
    block_rlbwt_data bd = prepare_blockrlbwt();
    measure_blockrlbwt("block_rlbwt_2",bd);
    measure_blockrlbwt("block_rlbwt_v",bd);
    measure_blockrlbwt("block_rlbwt_r",bd);
    cleanup_grlbwt(bd);
    */

    bench_index<pos_t,r_index,true,true>("r-index","r_index");
    //bench_index<pos_t,r_index_f<>,true,false>("r-index-f","r_index_f");
    //bench_index<pos_t,rcomp_glfig,true,true>("rcomp-glfig","rcomp_glfig");
    //bench_index<pos_t,online_rlbwt,true,true>("online-rlbwt","online_rlbwt");
}

void parse_args(char** argv, int argc, int &ptr) {
    std::string s = argv[ptr];
    ptr++;

    if (s == "-m") {
        if (ptr >= argc-1) help("error: missing parameter after -m option");
        path_mf = argv[ptr++];
        mf.open(path_mf,std::filesystem::exists(path_mf) ? std::ios::app : std::ios::out);
        if (!mf.good()) help("error: cannot open or create bench_indexment file");
    } else if (s == "-c") {
        if (ptr >= argc-1) help("error: missing parameter after -c option");
        check_correctness = true;
    } else {
        help("error: unrecognized '" + s +  "' option");
    }
}

int main_bench_a(int argc, char** argv) {
    if (argc == 8) {
        if (std::string(argv[2]) != "-m") help("");

        path_mf = argv[3];
        mf.open(path_mf,std::filesystem::exists(path_mf) ? std::ios::app : std::ios::out);
        if (!mf.good()) help("error: cannot open or create bench_indexment file");

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
    
    input_file.open(path_input_file);
    if (!input_file.good()) help("error: invalid input, could not read <input_file>");
    if (max_num_threads == 0 || max_num_threads > omp_get_max_threads()) help("error: invalid number of threads");

    std::cout << std::setprecision(4);
    name_text_file = path_input_file.substr(path_input_file.find_last_of("/\\") + 1);
    std::cout << "benchmarking " << path_input_file << std::flush;

    input_file.seekg(0,std::ios::end);
    input_size = input_file.tellg()+(std::streamsize)+1;
    input_file.seekg(0,std::ios::beg);
    no_init_resize(input,input_size);
    read_from_file(input_file,input.c_str(),input_size-1);
    preprocess_input();
    input[input_size-1] = 1;
    input_file.close();
    
    std::cout << " (" << format_size(input_size-1) << ") using up to " << format_threads(max_num_threads) << std::endl;

    if (input_size <= UINT_MAX) {
        if (input_size <= INT_MAX) {
            bench_a<uint32_t,int32_t>();
        } else {
            bench_a<uint32_t,int64_t>();
        }
    } else {
        bench_a<uint64_t,int64_t>();
    }

    if (mf.is_open()) mf.close();
    return 0;
}

int main_bench_indexes(int argc, char** argv) {
    if (argc < 4) help("");
    while (ptr < argc - 3) parse_args(argv, argc, ptr);

    path_input_file = argv[ptr];
    path_patterns_file_1 = argv[ptr+1];
    path_patterns_file_2 = argv[ptr+2];
    max_num_threads = 1;

    input_file.open(path_input_file);
    patterns_file_1.open(path_patterns_file_1);
    patterns_file_2.open(path_patterns_file_2);

    if (!input_file.good()) help("error: invalid input, could not read <input_file>");
    if (!patterns_file_1.good()) help("error: invalid input, could not read <patterns_file_1>");
    if (!patterns_file_2.good()) help("error: invalid input, could not read <patterns_file_2>");
    if (mf.is_open() && check_correctness) help("error: cannot output bench_indexment data when checking for correctness");

    system("chmod +x external/Big-BWT/*");
    system("chmod +x external/r-index/Big-BWT/*");
    system("chmod +x build/external/pfp-thresholds/*");
    system("chmod +x build/external/block_RLBWT/*");
    system("chmod +x build/external/grlBWT/*");
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

    std::cout << " (" << format_size(input_size-1) << ")" << std::endl;
    if (check_correctness) std::cout << "correctnes will be checked if possible" << std::endl;
    std::cout << std::endl;

    if (input_size <= UINT_MAX) {
        bench_indexes<uint32_t>();
    } else {
        bench_indexes<uint64_t>();
    }

    patterns_file_1.close();
    patterns_file_2.close();

    if (mf.is_open()) mf.close();
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 1) help("");
    if (std::string(argv[1]) == "-a") return main_bench_a(argc,argv);
    else return main_bench_indexes(argc,argv);
}