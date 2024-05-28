#include <iostream>
#include <filesystem>
#include <move_r/move_r.hpp>

int ptr = 1;
uint16_t p = 1;
bool revert_in_memory = false;
std::string path_index_file;
std::string path_outputfile;
std::string name_text_file;
std::ifstream index_file;
std::ofstream output_file;
std::ofstream mf;

void help(std::string msg) {
    if (msg != "") std::cout << msg << std::endl;
    std::cout << "move-r-revert: reconstruct the original file from the index." << std::endl << std::endl;
    std::cout << "usage: move-r-revert [options] <index_file> <output_file>" << std::endl;
    std::cout << "   -im                        revert in memory; faster, but stores the whole" << std::endl;
    std::cout << "                              output in memory" << std::endl;
    std::cout << "   -p <integer>               number of threads to use while reverting" << std::endl;
    std::cout << "                              (default: greatest possible)" << std::endl;
    std::cout << "   -m <m_file> <text_name>    m_file is the file to write measurement data to," << std::endl;
    std::cout << "                              text_name should be the name of the original file" << std::endl;
    std::cout << "   <index_file>               index file (with extension .move-r)" << std::endl;
    std::cout << "   <output_file>              output file" << std::endl;
    exit(0);
}

void parse_args(char** argv, int argc, int &ptr) {
    std::string s = argv[ptr];
    ptr++;

    if (s == "-m") {
        if (ptr >= argc-1) help("error: missing parameter after -m option");
        std::string path_m_file = argv[ptr++];
        mf.open(path_m_file,std::filesystem::exists(path_m_file) ? std::ios::app : std::ios::out);
        if (!mf.good()) help("error: cannot open measurement file");
        name_text_file = argv[ptr++];
    } else if (s == "-im") {
        revert_in_memory = true;
    } else if (s == "-p") {
        if (ptr >= argc-1) help("error: missing parameter after -p option");
        p = atoi(argv[ptr++]);
        if (p < 1) help("error: p < 1");
        if (p > omp_get_max_threads()) help("error: p > number of available threads");
    } else  {
        help("error: unrecognized '" + s + "' option");
    }
}

template <typename pos_t, move_r_locate_supp locate_support>
void measure_revert() {
    std::cout << std::setprecision(4);
    std::cout << "loading the index" << std::flush;
    auto t1 = now();
    move_r<locate_support,char,pos_t> index;
    index.load(index_file,{_revert});
    log_runtime(t1);
    index_file.close();
    std::cout << std::endl;
    index.log_data_structure_sizes();
    std::cout << std::endl;
    std::chrono::steady_clock::time_point t2,t3,t4;
    std::string input;
    p = std::min({(uint16_t)omp_get_max_threads(),index.max_revert_threads(),p});

    if (revert_in_memory) {
        std::cout << "reverting the index in memory using " << format_threads(p) << std::flush;
        t2 = now();
        input = index.revert({.num_threads = p});
        t3 = now();
        log_runtime(t2,t3);
        std::cout << "writing the input to the file " << std::flush;
        write_to_file(output_file,input.c_str(),input.size());
        t4 = now();
        log_runtime(t3,t4);
    } else {
        std::cout << "reverting the index using " << format_threads(p) << std::flush;
        t2 = now();
        index.revert(output_file,{.num_threads = p});
        t4 = now();
        log_runtime(t2,t4);
    }
    
    uint64_t time_revert = time_diff_ns(t2,t4);

    if (mf.is_open()) {
        mf << "RESULT";
        mf << " type=revert";
        mf << " text=" << name_text_file;
        mf << " a=" << index.balancing_parameter();
        mf << " n=" << index.input_size();
        mf << " sigma=" << std::to_string(index.alphabet_size());
        mf << " r=" << index.num_bwt_runs();
        mf << " r_=" << index.M_LF().num_intervals();
        mf << " r__=" << index.M_Phi().num_intervals();
        index.log_data_structure_sizes(mf);
        mf << " time_revert=" << time_revert;
        mf << std::endl;
        mf.close();
    }
}

int main(int argc, char **argv) {
    if (argc < 3) help("");
    while (ptr < argc - 2) parse_args(argv, argc, ptr);

    path_index_file = argv[ptr];
    path_outputfile = argv[ptr+1];

    index_file.open(path_index_file);
    output_file.open(path_outputfile);

    if (!index_file.good()) help("error: could not read <index_file>");
    if (!output_file.good()) help("error: could not create <output_file>");

    bool is_64_bit;
    index_file.read((char*)&is_64_bit,1);
    move_r_locate_supp _locate_support;
    index_file.read((char*)&_locate_support,sizeof(move_r_locate_supp));
    index_file.seekg(0,std::ios::beg);

    if (_locate_support == _phi) {
        if (is_64_bit) {
            measure_revert<uint64_t,_phi>();
        } else {
            measure_revert<uint32_t,_phi>();
        }
    } else {
        if (is_64_bit) {
            measure_revert<uint64_t,_rlzdsa>();
        } else {
            measure_revert<uint32_t,_rlzdsa>();
        }
    }

    output_file.close();
}