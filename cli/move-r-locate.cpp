#include <filesystem>
#include <iostream>
#include <move_r/move_r.hpp>

int ptr = 1;
bool output_occurrences = false;
bool check_correctness = false;
std::string input;
std::ofstream mf;
std::string path_index_file;
std::string path_patterns_file;
std::string path_text_file;
std::string path_outputfile;
std::ifstream index_file;
std::ifstream patterns_file;
std::ifstream input_file;
std::ofstream output_file;
std::string name_text_file;

void help(std::string msg)
{
    if (msg != "") std::cout << msg << std::endl;
    std::cout << "move-r-locate: locate all occurrences of the input patterns." << std::endl << std::endl;
    std::cout << "usage: move-r-locate [options] <index_file> <patterns>" << std::endl;
    std::cout << "   -c <input_file>            check correctness of each pattern occurrence on" << std::endl;
    std::cout << "                              this input file (must be the indexed input file)" << std::endl;
    std::cout << "   -m <m_file> <text_name>    m_file is the file to write measurement data to," << std::endl;
    std::cout << "                              text_name should be the name of the original file" << std::endl;
    std::cout << "   -o <output_file>           write pattern occurrences to this file (ASCII)" << std::endl;
    std::cout << "   <index_file>               index file (with extension .move-r)" << std::endl;
    std::cout << "   <patterns_file>            file in pizza&chili format containing the patterns" << std::endl;
    exit(0);
}

void parse_args(char** argv, int argc, int& ptr)
{
    std::string s = argv[ptr];
    ptr++;

    if (s == "-c") {
        if (ptr >= argc - 1)
            help("error: missing parameter after -c option.");

        check_correctness = true;
        path_text_file = argv[ptr++];
    } else if (s == "-m") {
        if (ptr >= argc - 1)
            help("error: missing parameter after -o option.");

        std::string path_m_file = argv[ptr++];
        mf.open(path_m_file, std::filesystem::exists(path_m_file) ? std::ios::app : std::ios::out);

        if (!mf.good())
            help("error: cannot open measurement file");

        name_text_file = argv[ptr++];
    } else if (s == "-o") {
        if (ptr >= argc - 1)
            help("error: missing parameter after -o option.");

        output_occurrences = true;
        path_outputfile = argv[ptr++];
    } else {
        help("error: unrecognized '" + s + "' option");
    }
}

template <typename pos_t, move_r_support support>
void measure_locate()
{
    std::cout << std::setprecision(4);
    std::cout << "loading the index" << std::flush;
    auto t1 = now();
    move_r<support, char, pos_t> index;
    index.load(index_file);
    log_runtime(t1);
    index_file.close();
    std::cout << std::endl;
    index.log_data_structure_sizes();

    if (check_correctness) {
        std::cout << std::endl << "reading the original input file" << std::flush;
        input_file.seekg(0, std::ios::end);
        uint64_t n = input_file.tellg() + (std::streamoff)1;
        input_file.seekg(0, std::ios::beg);
        no_init_resize(input, n - 1);
        read_from_file(input_file, input.c_str(), n - 1);
        input_file.close();
    }

    std::cout << std::endl << "searching patterns ... " << std::endl;
    std::string header;
    std::getline(patterns_file, header);
    uint64_t num_patterns = number_of_patterns(header);
    uint64_t pattern_length = patterns_length(header);
    uint64_t perc;
    uint64_t last_perc = 0;
    uint64_t num_occurrences = 0;
    uint64_t time_locate = 0;
    std::chrono::steady_clock::time_point t2, t3;
    std::string pattern;
    no_init_resize(pattern, pattern_length);
    std::vector<pos_t> occurrences;
    bool is_sorted, equal;
    pos_t count;

    for (uint64_t i = 0; i < num_patterns; i++) {
        perc = (100 * i) / num_patterns;
        is_sorted = false;

        if (perc > last_perc) {
            std::cout << perc << "% done .." << std::endl;
            last_perc = perc;
        }

        patterns_file.read((char*)&pattern[0], pattern_length);
        t2 = now();
        index.locate(pattern, occurrences);
        t3 = now();
        time_locate += time_diff_ns(t2, t3);
        num_occurrences += occurrences.size();

        if (check_correctness) {
            ips4o::sort(occurrences.begin(), occurrences.end());
            is_sorted = true;

            if (occurrences.size() != (count = index.count(pattern))) {
                std::cout << "error: wrong number of located occurrences: " << occurrences.size() << "/" << count << std::endl;
            }

            for (pos_t occurrence : occurrences) {
                equal = true;

                for (pos_t pos = 0; pos < pattern_length; pos++) {
                    if (input[occurrence + pos] != pattern[pos]) {
                        equal = false;
                        break;
                    }
                }

                if (!equal) {
                    std::cout << "error: wrong occurrence: " << occurrence << " (" << num_occurrences << " occurrences) " << std::endl;

                    for (pos_t pos = 0; pos < pattern_length; pos++)
                        std::cout << input[occurrence + pos];

                    std::cout << std::endl << std::endl << "/" << std::endl << std::endl;

                    for (pos_t pos = 0; pos < pattern_length; pos++)
                        std::cout << pattern[pos];
                        
                    std::cout << std::endl;
                    break;
                }
            }
        }

        if (output_occurrences) {
            if (!is_sorted)
                ips4o::sort(occurrences.begin(), occurrences.end());
                
            output_file.write((char*)&occurrences[0], occurrences.size());
        }

        occurrences.clear();
    }

    if (num_occurrences == 0) {
        std::cout << "found no occurrences" << std::endl;
    } else {
        std::cout << "average occurrences per pattern: " << (num_occurrences / num_patterns) << std::endl;
        std::cout << "number of patterns: " << num_patterns << std::endl;
        std::cout << "pattern length: " << pattern_length << std::endl;
        std::cout << "total number of occurrences: " << num_occurrences << std::endl;
        std::cout << "locate time: " << format_time(time_locate) << std::endl;
        std::cout << "             " << format_time(time_locate / num_patterns) << "/pattern" << std::endl;
        std::cout << "             " << format_time(time_locate / num_occurrences) << "/occurrence" << std::endl;
    }

    if (mf.is_open()) {
        mf << "RESULT";
        mf << " type=locate";
        mf << " text=" << name_text_file;
        mf << " a=" << index.balancing_parameter();
        mf << " n=" << index.input_size();
        mf << " sigma=" << std::to_string(index.alphabet_size());
        mf << " r=" << index.num_bwt_runs();
        mf << " r_=" << index.M_LF().num_intervals();

        if constexpr (support == _locate_move) {
            mf << " r__=" << index.M_Phi_m1().num_intervals();
        } else if constexpr (support == _locate_rlzdsa) {
            mf << " z__=" << index.num_phrases_rlzdsa();
            mf << " z_l_=" << index.num_literal_phrases_rlzdsa();
            mf << " z_c_=" << index.num_copy_phrases_rlzdsa();
        }

        mf << " pattern_length=" << pattern_length;
        index.log_data_structure_sizes(mf);
        mf << " num_patterns=" << num_patterns;
        mf << " num_occurrences=" << num_occurrences;
        mf << " time_locate=" << time_locate;
        mf << std::endl;
        mf.close();
    }
}

int main(int argc, char** argv)
{
    if (argc < 3)
        help("");
    while (ptr < argc - 2)
        parse_args(argv, argc, ptr);

    path_index_file = argv[ptr];
    path_patterns_file = argv[ptr + 1];

    index_file.open(path_index_file);
    patterns_file.open(path_patterns_file);

    if (!index_file.good())
        help("error: could not read <index_file>");
    if (!patterns_file.good())
        help("error: could not read <patterns_file>");

    if (output_occurrences) {
        output_file.open(path_outputfile);
        if (!output_file.good())
            help("error: could not create <output_file>");
    }

    if (check_correctness) {
        input_file.open(path_text_file);
        if (!input_file.good())
            help("error: could not read <input_file>");
    }

    bool is_64_bit;
    index_file.read((char*)&is_64_bit, 1);
    move_r_support _support;
    index_file.read((char*)&_support, sizeof(move_r_support));
    index_file.seekg(0, std::ios::beg);

    if (_support == _count || _support == _locate_one) {
        std::cout << "error: this index does not support locate" << std::endl;
        exit(0);
    } else if (_support == _locate_move) {
        if (is_64_bit) {
            measure_locate<uint64_t, _locate_move>();
        } else {
            measure_locate<uint32_t, _locate_move>();
        }
    } else {
        if (is_64_bit) {
            measure_locate<uint64_t, _locate_rlzdsa>();
        } else {
            measure_locate<uint32_t, _locate_rlzdsa>();
        }
    }

    patterns_file.close();
    if (output_occurrences)
        output_file.close();
}