#include <iostream>
#include <filesystem>
#include <move_r/misc/utils.hpp>
#include <move_r/move_r.hpp>

int ptr = 1;
std::ofstream measurement_file;
std::string path_index_file;
std::string path_patterns_file;
std::ifstream index_file;
std::ifstream patterns_file;
std::string name_textfile;

void help(std::string msg) {
    if (msg != "") std::cout << msg << std::endl;
    std::cout << "move-r-count: count all occurrences of the input patterns." << std::endl << std::endl;
    std::cout << "usage: move-r-count <index_file> <patterns_file>" << std::endl;
    std::cout << "   -m <m_file> <text_name>    m_file is the file to write measurement data to," << std::endl;
    std::cout << "                              text_name should be the name of the original file" << std::endl;
    std::cout << "   <index_file>               index file (with extension .move-r)" << std::endl;
    std::cout << "   <patterns_file>            file in pizza&chili format containing the patterns." << std::endl;
    exit(0);
}

void parse_args(char **argv, int argc, int &ptr) {
    std::string s = argv[ptr];
    ptr++;

    if (s == "-m") {
        if (ptr >= argc - 1) help("error: missing parameter after -o option.");
        std::string path_m_file = argv[ptr++];
        measurement_file.open(path_m_file,std::filesystem::exists(path_m_file) ? std::ios::app : std::ios::out);
        if (!measurement_file.good()) help("error: cannot open measurement file");
        name_textfile = argv[ptr++];
    } else {
        help("error: unrecognized '" + s + "' option");
    }
}

template <typename uint_t>
void measure_count() {
    std::cout << std::setprecision(4);
    std::cout << "loading the index" << std::flush;
    auto t1 = now();
    move_r<uint_t> index;
    index.load(index_file,{move_r_support::count});
    log_runtime(t1);
    index_file.close();
    std::cout << std::endl;
    index.log_data_structure_sizes();
    std::cout << std::endl << "searching patterns ... " << std::endl;
    std::string header;
    std::getline(patterns_file,header);
    uint64_t num_patterns = number_of_patterns(header);
    uint64_t pattern_length = patterns_length(header);
    uint64_t perc;
    uint64_t last_perc = 0;
    uint64_t num_occurrences = 0;
    uint64_t time_count = 0;
    std::chrono::steady_clock::time_point t2,t3;
    std::string pattern;
    no_init_resize(pattern,pattern_length);

    for (uint64_t i=0; i<num_patterns; i++) {
        perc = (100*i) / num_patterns;

        if (perc > last_perc) {
            std::cout << perc << "% done .." << std::endl;
            last_perc = perc;
        }

        patterns_file.read((char*)&pattern[0],pattern_length);
        t2 = now();
        num_occurrences += index.count(pattern);
        t3 = now();
        time_count += time_diff_ns(t2,t3);
    }

    patterns_file.close();

    if (num_occurrences == 0) {
        std::cout << "found no occurrences" << std::endl;
    } else {
        std::cout << "average occurrences per pattern: " << (num_occurrences/num_patterns) << std::endl;
        std::cout << "number of patterns: " << num_patterns << std::endl;
        std::cout << "pattern length: " << pattern_length << std::endl;
        std::cout << "total number of occurrences: " << num_occurrences << std::endl;
        std::cout << "count time: " << format_time(time_count) << std::endl;
        std::cout << "            " << format_time(time_count/num_patterns) << "/pattern" << std::endl;
        std::cout << "            " << format_time(time_count/num_occurrences) << "/occurrence" << std::endl;
    }

    if (measurement_file.is_open()) {
        measurement_file << "RESULT";
        measurement_file << " type=count";
        measurement_file << " text=" << name_textfile;
        measurement_file << " index_impl=move_r";
        measurement_file << " a=" << index.balancing_parameter();
        measurement_file << " n=" << index.input_size();
        measurement_file << " sigma=" << std::to_string(index.alphabet_size());
        measurement_file << " r=" << index.num_bwt_runs();
        measurement_file << " r_=" << index.num_intervals_m_lf();
        measurement_file << " r__=" << index.num_intervals_m_phi();
        measurement_file << " pattern_length=" << pattern_length;
        index.log_data_structure_sizes(measurement_file);
        measurement_file << " num_patterns=" << num_patterns;
        measurement_file << " num_occurrences=" << num_occurrences;
        measurement_file << " time_count=" << time_count;
        measurement_file << std::endl;
        measurement_file.close();
    }
}

int main(int argc, char **argv) {
    if (argc < 3) help("");
    while (ptr < argc - 2) parse_args(argv, argc, ptr);

    path_index_file = argv[ptr];
    path_patterns_file = argv[ptr+1];

    index_file.open(path_index_file);
    patterns_file.open(path_patterns_file);

    if (!index_file.good()) help("error: could not read <index_file>");
    if (!patterns_file.good()) help("error: could not read <patterns_file>");

    bool is_64_bit;
    index_file.read((char*)&is_64_bit,1);
    index_file.seekg(0,std::ios::beg);

    if (is_64_bit) {
        measure_count<uint64_t>();
    } else {
        measure_count<uint32_t>();
    }
}