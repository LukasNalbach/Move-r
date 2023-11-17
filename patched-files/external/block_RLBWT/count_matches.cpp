//#define VERB

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <tuple>

#include "include/reader.hpp"
#include "include/types.hpp"

void help() {
    std::cout << "count matches in RLBWT data structure.\n\n";
    std::cout << "Usage: count_matches [options] <bwt_file> <patterns>\n";
    std::cout << "   bwt_file   Path to block rlbwt element.\n";
    std::cout << "   patterns   Path to file containing patterns.\n";
    std::cout << "   -s         Block rlbwt is space optimized.\n";
    std::cout << "   -c         Blocks contains a constant number of runs.\n";
    std::cout << "   -t         Don't include query times in std::cout\n";
    std::cout << "Bwt and pattern files are required.\n\n";
    std::cout << "Example: count_matches bwt.bin Einstein.txt >> /dev/null" << std::endl;
    exit(0);
}

template <class bwt_type>
std::tuple<uint64_t, uint64_t, uint64_t> bench(const std::string& in_file_path, std::ifstream& patterns, bool o_t, double& bps) {
    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::nanoseconds;

    bwt_type bwt(in_file_path);
    bps = 8 * double(bwt.bytes()) / bwt.size();
    std::string p;
    size_t i = 0;
    uint64_t time_query = 0;
    uint64_t num_queries = 0;
    uint64_t num_occurrenes = 0;
    while (std::getline(patterns, p)) {
        auto start = high_resolution_clock::now();
        num_occurrenes += bwt.count(p);
        auto end = high_resolution_clock::now();
        time_query += duration_cast<nanoseconds>(end - start).count();
        /*
        if (o_t) {
            std::cout << p << "\t" << count << "\t" << time << std::endl;
        } else {
            std::cout << p << "\t" << count << std::endl;
        }
        */
        num_queries++;
    }
    return {time_query, num_queries, num_occurrenes};
}

int main(int argc, char const* argv[]) {
    if (argc < 3) {
        std::cerr << "Input and pattern files are required\n" << std::endl;
        help();
    }
    std::string in_file_path = "";
    std::string patterns = "";
    bool space_op = false;
    bool run_block = false;
    bool output_time = true;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-s") == 0) {
            space_op = true;
        } else if (strcmp(argv[i], "-c") == 0) {
            run_block = true;
        } else if (strcmp(argv[i], "-t") == 0) {
            output_time = false;
        } else if (in_file_path.size() == 0) {
            in_file_path = argv[i];
        } else {
            patterns = argv[i];
        }
    }
    std::ifstream p(patterns);
    std::cerr << "looking for patterns from " << patterns << " in " << in_file_path << std::endl;
    //std::cout << "Pattern\tcount\ttime" << std::endl;
    std::tuple<uint64_t, uint64_t, uint64_t> res;
    double bps = 0;
    if (run_block) {
        res = bench<bbwt::run<>>(in_file_path, p, output_time, bps);
    } else if (space_op) {
        res = bench<bbwt::vbyte<>>(in_file_path, p, output_time, bps);
    } else {
        res = bench<bbwt::two_byte<>>(in_file_path, p, output_time, bps);
    }
    std::cerr << "Mean query time: " << std::get<0>(res) << " / " << std::get<1>(res) << " = " << (std::get<0>(res) / std::get<1>(res)) << "ns\n" 
              << " with " << bps << " bits per symbol" << std::endl;
    std::cerr << "total number of occurrences: " << std::get<2>(res) << std::endl;
    std::ofstream results_file("blockrlbwt_count_result.txt");
    results_file.write((char*)&std::get<0>(res),sizeof(uint64_t));
    results_file.write((char*)&std::get<1>(res),sizeof(uint64_t));
    results_file.write((char*)&std::get<2>(res),sizeof(uint64_t));
    results_file.close();
}