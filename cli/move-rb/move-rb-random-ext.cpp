/**
 * part of LukasNalbach/Move-r
 * Adapted for random bidirectional extensions with COUNT/LOCATE switch.
 */

#include <filesystem>
#include <iostream>
#include <random>
#include <move_rb/move_rb.hpp>

static constexpr int min_args = 2;
int arg_idx = 1;
std::ofstream mf;
std::string path_index_file;
std::string path_patterns_file;
std::ifstream index_file;
std::ifstream patterns_file;
std::string name_text_file;
move_rb_query_support_t query_mode = COUNT;

void help(std::string msg)
{
    if (msg != "") std::cout << msg << std::endl;
    std::cout << "move-rb-random-ext: search all patterns using a random sequence of left/right extensions." << std::endl << std::endl;
    std::cout << "usage: move-rb-random-ext [...] [-q <count|locate>] <index_file> <patterns_file>" << std::endl;
    std::cout << "   -q <mode>                  query mode to use: 'count' or 'locate' (default: count)" << std::endl;
    std::cout << "   -m <m_file> <text_name>    m_file is the file to write measurement data to," << std::endl;
    std::cout << "                              text_name should be the name of the original file" << std::endl;
    exit(0);
}

bool parse_args(char** argv, int argc)
{
    if (arg_idx >= argc - min_args) return false;
    std::string s = argv[arg_idx];

    if (s == "-m") {
        arg_idx++;
        if (arg_idx >= argc - 1) help("error: missing parameter after -m option.");
        std::string path_m_file = argv[arg_idx++];
        mf.open(path_m_file, std::filesystem::exists(path_m_file) ? std::ios::app : std::ios::out);
        if (!mf.good()) help("error: cannot open nor create <m_file>");
        name_text_file = argv[arg_idx++];
    } else if (s == "-q") {
        arg_idx++;
        if (arg_idx >= argc - 1) help("error: missing parameter after -q option.");
        std::string q_str = argv[arg_idx++];
        if (q_str == "count") query_mode = COUNT;
        else if (q_str == "locate") query_mode = LOCATE;
        else help("error: invalid query mode, use 'count' or 'locate'.");
    } else return false;
    
    return true;
}

template <typename pos_t, move_r_support support, move_rb_query_support_t query_support>
void measure_random_ext()
{
    std::cout << std::setprecision(4);
    std::cout << "loading the index" << std::flush;
    auto time = now();
    using idx_t = move_rb<support, char, pos_t>;
    idx_t index;
    index.load(index_file);
    index_file.close();
    time = log_runtime(time);
    index.log_data_structure_sizes();

    std::cout << std::endl << "searching patterns via random extensions ... " << std::endl;
    std::string header;
    std::getline(patterns_file, header);
    uint64_t num_patterns = number_of_patterns(header);
    uint64_t pattern_length = patterns_length(header);

    uint64_t perc;
    uint64_t last_perc = 0;
    uint64_t num_occurrences = 0;
    uint64_t time_search = 0;
    std::chrono::steady_clock::time_point t2, t3;
    std::string pattern;
    no_init_resize(pattern, pattern_length);
    uint64_t checksum = 0;
    uint64_t baseline_alloc = malloc_count_current();
    malloc_count_reset_peak();

    uint64_t seed = num_patterns * 131ULL + pattern_length * 37ULL + 7ULL;
    std::mt19937 gen(seed);

    for (uint64_t i = 0; i < num_patterns; i++) {
        perc = (100 * i) / num_patterns;

        if (perc > last_perc) {
            std::cout << perc << "% done .." << std::endl;
            last_perc = perc;
        }

        patterns_file.read(pattern.data(), pattern_length);
        t2 = now();

        std::uniform_int_distribution<uint64_t> start_dist(0, pattern_length - 1);
        uint64_t start_pos = start_dist(gen);
        typename idx_t::template search_context_t<query_support> ctx(index);
        auto [next_ctx, success] = ctx.template extend<LEFT>(index, pattern[start_pos]);
        ctx = next_ctx;
        uint64_t left_bound = start_pos;
        uint64_t right_bound = start_pos;

        while (success && (right_bound - left_bound + 1 < pattern_length)) {
            direction_t dir;
            uint64_t next_char_idx;

            if (left_bound > 0 && right_bound < pattern_length - 1) {
                std::uniform_int_distribution<int> dir_dist(0, 1);

                if (dir_dist(gen) == 0) {
                    dir = LEFT;
                    left_bound--;
                    next_char_idx = left_bound;
                } else {
                    dir = RIGHT;
                    right_bound++;
                    next_char_idx = right_bound;
                }
            } else if (left_bound > 0) {
                dir = LEFT;
                left_bound--;
                next_char_idx = left_bound;
            } else {
                dir = RIGHT;
                right_bound++;
                next_char_idx = right_bound;
            }

            auto [extended_ctx, ext_success] = ctx.extend(index, pattern[next_char_idx], dir);
            success = ext_success;
            if (success) {
                ctx = extended_ctx;
            }
        }

        if (success) {
            num_occurrences += ctx.num_occ();

            if constexpr (query_support == LOCATE) {
                auto loc_ctx = ctx.locate_phase();
                loc_ctx.locate(index, ctx, [&](pos_t occ) {
                    checksum += occ;
                });
            } else {
                auto [b, e] = ctx.forward_sa_interval();
                checksum += b + e;
            }
        }

        t3 = now();
        time_search += time_diff_ns(t2, t3);
    }

    patterns_file.close();
    std::cout << "checksum: " << checksum << std::endl;
    std::cout << "average occurrences per pattern: " << (num_patterns ? (num_occurrences / num_patterns) : 0) << std::endl;
    std::cout << "number of patterns: " << num_patterns << std::endl;
    std::cout << "pattern length: " << pattern_length << std::endl;
    std::cout << "query mode: " << (query_support == LOCATE ? "locate" : "count") << std::endl;
    std::cout << "total number of occurrences: " << num_occurrences << std::endl;
    std::cout << "search time: " << format_time(time_search) << std::endl;
    std::cout << "             " << format_time(time_search / num_patterns) << "/pattern" << std::endl;
    std::cout << "             " << format_time(time_search / (num_patterns * pattern_length)) << "/character" << std::endl;
    std::cout << "             " << format_time(time_search / num_occurrences) << "/occurrence" << std::endl;

    if (mf.is_open()) {
        mf << "RESULT";
        mf << " algo=random_ext_" << (query_support == LOCATE ? "locate" : "count") << "_move_rb_" << move_r_support_suffix(support);
        mf << " text=" << name_text_file;
        mf << " n=" << index.forward_index().input_size();
        mf << " pattern_length=" << pattern_length;
        mf << " num_patterns=" << num_patterns;
        mf << " num_occurrences=" << num_occurrences;
        mf << " time_search=" << time_search;
        index.log_data_structure_sizes(mf);
        mf << std::endl;
        mf.close();
    }
}

int main(int argc, char** argv)
{
    if (argc - 1 < min_args) help("");
    while (parse_args(argv, argc));

    path_index_file = argv[arg_idx];
    path_patterns_file = argv[arg_idx + 1];

    index_file.open(path_index_file);
    patterns_file.open(path_patterns_file);

    if (!index_file.good()) help("error: could not read <index_file>");
    if (!patterns_file.good()) help("error: could not read <patterns_file>");

    bool is_64_bit;
    index_file.read((char*) &is_64_bit, 1);
    move_r_support _support;
    index_file.read((char*) &_support, sizeof(move_r_support));
    index_file.seekg(0, std::ios::beg);

    if (_support == _count) {
        if (query_mode == LOCATE) help("error: this index does not support locate.");
        if (is_64_bit) measure_random_ext<uint64_t, _count, COUNT>();
        else           measure_random_ext<uint32_t, _count, COUNT>();
    } else if (_support == _locate_move) {
        if (is_64_bit) {
            if (query_mode == LOCATE) measure_random_ext<uint64_t, _locate_move, LOCATE>();
            else                      measure_random_ext<uint64_t, _locate_move, COUNT>();
        } else {
            if (query_mode == LOCATE) measure_random_ext<uint32_t, _locate_move, LOCATE>();
            else                      measure_random_ext<uint32_t, _locate_move, COUNT>();
        }
    } else if (_support == _locate_rlzsa) {
        if (is_64_bit) {
            if (query_mode == LOCATE) measure_random_ext<uint64_t, _locate_rlzsa, LOCATE>();
            else                      measure_random_ext<uint64_t, _locate_rlzsa, COUNT>();
        } else {
            if (query_mode == LOCATE) measure_random_ext<uint32_t, _locate_rlzsa, LOCATE>();
            else                      measure_random_ext<uint32_t, _locate_rlzsa, COUNT>();
        }
    }

    return 0;
}