#include <iostream>
#include <filesystem>
#include <move_r/misc/utils.hpp>
#include <move_r/move_r.hpp>

std::string pathprefix_indexfile;
move_r_construction_mode construction_mode = move_r_construction_mode::runtime;
std::vector<move_r_support> support = full_support;
uint16_t p = omp_get_max_threads();
uint16_t a = 8;
uint64_t n;
std::ofstream measurement_file_index;
std::ofstream m_file_mds;
std::ifstream input_file;
std::ofstream index_file;
std::string path_inputfile;
std::string name_textfile;
std::string path_index_file;
int ptr = 1;

void help(std::string msg) {
    if (msg != "") std::cout << msg << std::endl;
    std::cout << "move-r-build: builds move-r." << std::endl << std::endl;
    std::cout << "usage: move-r-build [options] <input_file>" << std::endl;
    std::cout << "   -c <mode>          construction mode: runtime or space (default: runtime)" << std::endl;
    std::cout << "   -o <base_name>     names the index file base_name.move-r (default: input_file)" << std::endl;
    std::cout << "   -s <op1> <op2> ... supported operations: revert, count and locate" << std::endl;
    std::cout << "                      (default: all operations)" << std::endl;
    std::cout << "   -p <integer>       number of threads to use during the construction of the index" << std::endl;
    std::cout << "                      (default: all threads)" << std::endl;
    std::cout << "   -a <integer>       balancing parameter; a must be an integer number and a >= 2 (default: 8)" << std::endl;
    std::cout << "   -m_idx <m_file>    m_file is file to write measurement data of the index construction to" << std::endl;
    std::cout << "   -m_mds <m_file>    m_file is file to write measurement data of the construction of the move" << std::endl;
    std::cout << "                      data structures to" << std::endl;
    std::cout << "   <input_file>       input file" << std::endl;
    exit(0);
}

void parse_args(char** argv, int argc, int &ptr) {
    std::string s = argv[ptr];
    ptr++;

    if (s == "-o") {
        if (ptr >= argc-1) help("error: missing parameter after -o option");
        pathprefix_indexfile = argv[ptr++];
    } else if (s == "-p") {
        if (ptr >= argc-1) help("error: missing parameter after -p option");
        p = atoi(argv[ptr++]);
        if (p < 1) help("error: p < 1");
        if (p > omp_get_max_threads()) help("error: p > maximum number of threads");
    } else if (s == "-c") {
        if (ptr >= argc-1) help("error: missing parameter after -p option");
        std::string construction_mode_str = argv[ptr++];
        if (construction_mode_str == "runtime") construction_mode = move_r_construction_mode::runtime;
        else if (construction_mode_str == "space") construction_mode = move_r_construction_mode::space;
        else help("error: invalid option for -c");
    } else if (s == "-s") {
        if (ptr >= argc-1) help("error: missing parameter after -s option");
        support.clear();
        bool parsed_first_op = false;

        while (true) {
            std::string support_str = argv[ptr++];
            if (support_str == "revert") support.emplace_back(move_r_support::revert);
            else if (support_str == "count") support.emplace_back(move_r_support::count);
            else if (support_str == "locate") support.emplace_back(move_r_support::locate);
            else if (!parsed_first_op) help("error: unknown mode provided with -s option");
            else break;
            parsed_first_op = true;
        }
    } else if (s == "-a") {
        if (ptr >= argc-1) help("error: missing parameter after -a option");
        a = atoi(argv[ptr++]);
        if (a < 2) help("error: a < 2");
    } else if (s == "-m_idx") {
        if (ptr >= argc-1) help("error: missing parameter after -m_idx option");
        std::string path_measurement_file_index = argv[ptr++];
        measurement_file_index.open(path_measurement_file_index,std::filesystem::exists(path_measurement_file_index) ? std::ios::app : std::ios::out);
        if (!measurement_file_index.good()) help("error: cannot open or create measurement file");
    } else if (s == "-m_mds") {
        if (ptr >= argc-1) help("error: missing parameter after -m_mds option");
        std::string path_m_file_mds = argv[ptr++];
        m_file_mds.open(path_m_file_mds,std::filesystem::exists(path_m_file_mds) ? std::ios::app : std::ios::out);
        if (!m_file_mds.good()) help("error: cannot open or create at least one measurement file");
    } else {
        help("error: unrecognized '" + s + "' option");
    }
}

template <typename uint_t>
void build() {
    move_r<uint_t> index(
        input_file,support,construction_mode,p,a,true,
        measurement_file_index.is_open() ? &measurement_file_index : NULL,
        m_file_mds.is_open() ? &m_file_mds : NULL,
        path_inputfile
    );
    input_file.close();
    auto time = now();
    index.serialize(index_file,support);
    std::cout << "serializing the index" << std::flush;
    time = log_runtime(time);
}

int main(int argc, char** argv) {
    if (argc < 2) help("");
    while (ptr < argc - 1) parse_args(argv, argc, ptr);
    path_inputfile = argv[ptr];
    if (pathprefix_indexfile == "") pathprefix_indexfile = path_inputfile;

    std::cout << std::setprecision(4);
    name_textfile = path_inputfile.substr(path_inputfile.find_last_of("/\\") + 1);
    path_index_file = pathprefix_indexfile.append(".move-r");
    
    std::cout << "building move-r of " << path_inputfile;
    std::cout << " using " << format_threads(p) << " and a = " << a << std::endl;
    std::cout << "the index will be saved to " << path_index_file << std::endl << std::endl;

    input_file.open(path_inputfile);
    index_file.open(path_index_file);

    if (!input_file.good()) help("error: invalid input, could not read <input_file>");
    if (!index_file.good()) help("error: invalid input, could not create <index_file>");

    input_file.seekg(0,std::ios::end);
    n = input_file.tellg()+(std::streamsize)+1;
    input_file.seekg(0,std::ios::beg);

    if (p > 1 && 1000*p > n) {
        p = std::max((uint64_t)1,n/1000);
        std::cout << "n = " << n << ", warning: p > n/1000, setting p to n/1000 ~ "<< std::to_string(p) << std::endl;
    }

    if (measurement_file_index.is_open()) {
        measurement_file_index << "RESULT"
            << " type=build_index"
            << " text=" << name_textfile
            << " index_impl=move_r"
            << " p=" << p
            << " a=" << a;
    }

    if (n < UINT_MAX) {
        build<uint32_t>();
    } else {
        build<uint64_t>();
    }

    if (measurement_file_index.is_open()) measurement_file_index.close();
    if (m_file_mds.is_open()) m_file_mds.close();
    index_file.close();
}