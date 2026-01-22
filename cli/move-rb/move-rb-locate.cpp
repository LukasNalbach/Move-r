/**
 * part of LukasNalbach/Move-r
 *
 * MIT License
 *
 * Copyright (c) Lukas Nalbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <filesystem>
#include <iostream>
#include <move_r/move_rb.hpp>

static constexpr int min_args = 6;
int arg_idx = 1;
int64_t k = -1;
std::string scheme_str;
distance_metric_t dist_metr = NO_METRIC;
search_scheme_t search_scheme;
bool output_occurrences = false;
bool filter_occurrences = false;
bool check_correctness = false;
std::string input;
std::ofstream mf;
std::string path_index_file;
std::string path_patterns_file;
std::string path_output_file;
std::ifstream index_file;
std::ifstream patterns_file;
std::ifstream input_file;
std::ofstream output_file;
std::string name_text_file;
std::string path_input_file;

void help(std::string msg)
{
    if (msg != "") std::cout << msg << std::endl;
    std::cout << "move-rb-locate: locate all approximate occurrences of the input patterns." << std::endl << std::endl;
    std::cout << "usage: move-rb-locate [...] -d <metric> -s <scheme> [-k <mismatches>] <index_file> <patterns_file>" << std::endl;
    std::cout << "   -m <m_file> <text_name>    m_file is the file to write measurement data to," << std::endl;
    std::cout << "                              text_name should be the name of the original file" << std::endl;
    std::cout << "   -i <input_file>            input_file must be the file the index was built for" << std::endl;
    std::cout << "                              (required for the -c option)" << std::endl;
    std::cout << "   -c                         checks correctness of each pattern occurrence on <input_file>" << std::endl;
    std::cout << "   -f                         filters out redundant occurrdences" << std::endl;
    std::cout << "   <metric>                   distance metric to use (hamming or edit)" << std::endl;
    std::cout << "   <scheme>                   search scheme to use (pigeon_hole, suffix_filter, 01 or path to a file)" << std::endl;
    std::cout << "   <mismatches>               maximum number of allowed mismatches; applies only to" << std::endl;
    std::cout << "                              pigeon_hole suffix_filter and 01 search schemes" << std::endl;
    std::cout << "   -o <output_file>           write pattern occurrences to this file (in ASCII format; one line per pattern)" << std::endl;
    std::cout << "   <index_file>               index file (with extension .move-r)" << std::endl;
    std::cout << "   <patterns_file>            file in move-rb-patterns format containing the patterns." << std::endl;
    exit(0);
}

bool parse_args(char** argv, int argc)
{
    if (arg_idx >= argc - min_args) return false;
    std::string s = argv[arg_idx];
    if (s == "-d") return false;
    arg_idx++;

    if (s == "-f") {
        filter_occurrences = true;
    } else if (s == "-c") {
        check_correctness = true;
    } else if (s == "-m") {
        if (arg_idx >= argc - 1) help("error: missing parameter after -o option.");
        std::string path_m_file = argv[arg_idx++];
        mf.open(path_m_file, std::filesystem::exists(path_m_file) ? std::ios::app : std::ios::out);
        if (!mf.good()) help("error: cannot open nor create <m_file>");
        name_text_file = argv[arg_idx++];
    } else if (s == "-i") {
        if (arg_idx >= argc - 1) help("error: missing parameter after -i option.");
        path_input_file = argv[arg_idx++];
        input_file.open(path_input_file);
        if (!input_file.good()) help("error: cannot open <input_file>");
    }  else if (s == "-o") {
        if (arg_idx >= argc - 1) help("error: missing parameter after -o option.");
        output_occurrences = true;
        path_output_file = argv[arg_idx++];
    } else {
        help("error: unrecognized '" + s + "' option");
    }

    return true;
}

template <typename pos_t, move_r_support support>
void measure_locate()
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

    if (check_correctness) {
        if (path_input_file == "") help("error: <input_file> not provided");
        std::cout << std::endl << "loading input file" << std::flush;
        no_init_resize(input, index.forward_index().input_size());
        read_from_file(input_file, input.data(), input.size());
        input_file.close();
        time = log_runtime(time);
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
    std::vector<aprx_occ_t<pos_t>> occurrences;
    std::vector<aprx_occ_t<pos_t>> redundant_occurrences;
    uint64_t checksum = 0;
    uint64_t baseline_alloc = malloc_count_current();
    malloc_count_reset_peak();

    for (uint64_t i = 0; i < num_patterns; i++) {
        perc = (100 * i) / num_patterns;

        if (perc > last_perc) {
            std::cout << perc << "% done .." << std::endl;
            last_perc = perc;
        }

        patterns_file.read(pattern.data(), pattern_length);
        t2 = now();

        if (dist_metr == HAMMING_DISTANCE) {
            index.locate_hamming_dist(pattern, search_scheme, [&](aprx_occ_t<pos_t> occ){redundant_occurrences.emplace_back(occ);});
        } else {
            index.locate_edit_dist(pattern, search_scheme, [&](aprx_occ_t<pos_t> occ){redundant_occurrences.emplace_back(occ);});
        }

        if (filter_occurrences) {
            ips4o::sort(redundant_occurrences.begin(), redundant_occurrences.end());
            static constexpr pos_t infty = std::numeric_limits<pos_t>::max();
            aprx_occ_t<pos_t> prev {.pos = infty, .len = infty, .err = k + 1};
            int64_t window = 2 * k;

            for (const auto& occ : redundant_occurrences) {
                int64_t dist = abs_diff<int64_t>(occ.pos, prev.pos);
                if (dist == 0) continue;

                if (dist <= window) {
                    if (occ.err > prev.err ||
                       (occ.err == prev.err && occ.len >= prev.len)
                    ) continue;

                    occurrences.pop_back();
                }

                occurrences.emplace_back(occ);
                prev = occ;
            }

            redundant_occurrences.clear();
        } else {
            std::swap(occurrences, redundant_occurrences);
        }
        
        num_occurrences += occurrences.size();

        t3 = now();
        time_locate += time_diff_ns(t2, t3);

        if (check_correctness) {
            for (auto occ : occurrences) {
                std::string_view occ_view(input.c_str() + occ.pos, occ.len);
                pos_t dist;

                if (dist_metr == HAMMING_DISTANCE) {
                    dist = hamming_dist_bounded<pos_t>(occ_view, pattern, k);
                } else {
                    dist = edit_dist_bounded<pos_t>(occ_view, pattern, k);
                }

                if (dist > k || occ.err != dist) {
                    std::cout << "error: wrong approximate occurrence " << occ.pos <<
                        " with length " << occ.len << " of pattern '" << pattern << "'" << std::endl;
                    exit(-1);
                }
            }
        }

        if (output_occurrences) {
            if (!filter_occurrences) {
                ips4o::sort(occurrences.begin(), occurrences.end());
            }

            for (auto occ : occurrences) {
                output_file <<
                    "(pos=" << occ.pos <<
                    " len=" << occ.len <<
                    " err=" << occ.err << ") ";
            }

            output_file << std::endl;
        }

        checksum += occurrences.size();

        for (auto occ : occurrences) {
            checksum += occ.pos + occ.err;
        }

        occurrences.clear();
    }

    patterns_file.close();
    std::cout << "checksum: " << checksum << std::endl;
    std::cout << "additional memory consumption during the search phase: " << format_size(malloc_count_peak() - baseline_alloc) << std::endl;
    std::cout << "average occurrences per pattern: " << (num_occurrences / num_patterns) << std::endl;
    std::cout << "number of patterns: " << num_patterns << std::endl;
    std::cout << "pattern length: " << pattern_length << std::endl;
    std::cout << "maximum number of mismatches: " << k << std::endl;
    std::cout << "search scheme: " << scheme_str << std::endl;
    std::cout << "total number of occurrences: " << num_occurrences << std::endl;
    std::cout << "locate time: " << format_time(time_locate) << std::endl;
    std::cout << "             " << format_time(time_locate / num_patterns) << "/pattern" << std::endl;
    std::cout << "             " << format_time(time_locate / (num_patterns * pattern_length)) << "/character" << std::endl;
    std::cout << "             " << format_time(time_locate / num_occurrences) << "/occurrence" << std::endl;

    if (mf.is_open()) {
        mf << "RESULT";
        mf << " algo=count_move_rb_" << move_r_support_suffix(support);
        mf << " text=" << name_text_file;
        mf << " a=" << index.forward_index().balancing_parameter();
        mf << " n=" << index.forward_index().input_size();
        mf << " sigma=" << std::to_string(index.forward_index().alphabet_size());
        mf << " r=" << index.forward_index().num_bwt_runs();
        mf << " r_rev=" << index.backward_index().num_bwt_runs();
        mf << " r_=" << index.forward_index().M_LF().num_intervals();
        mf << " r_rev_=" << index.backward_index().M_LF().num_intervals();

        if constexpr (idx_t::supports_multiple_locate) {
            if constexpr (idx_t::has_locate_move) {
                mf << " r__=" << index.forward_index().M_Phi_m1().num_intervals();
                mf << " r___=" << index.forward_index().M_Phi().num_intervals();
            } else if constexpr (idx_t::has_rlzsa) {
                mf << " z=" << index.forward_index().num_phrases_rlzsa();
                mf << " z_l=" << index.forward_index().num_literal_phrases_rlzsa();
                mf << " z_c=" << index.forward_index().num_copy_phrases_rlzsa();
            }
        }

        mf << " pattern_length=" << pattern_length;
        index.log_data_structure_sizes(mf);
        mf << " num_patterns=" << num_patterns;
        mf << " max_mismatches=" << k;
        mf << " num_occurrences=" << num_occurrences;
        mf << " time_locate=" << time_locate;
        mf << std::endl;
        mf.close();
    }
}

int main(int argc, char** argv)
{
    if (argc - 1 < min_args) help("");
    while (parse_args(argv, argc));

    std::string arg = argv[arg_idx++];
    if (arg != "-d") help("");
    std::string dist_str = argv[arg_idx++];
    if      (dist_str == "hamming") dist_metr = HAMMING_DISTANCE;
    else if (dist_str == "edit")    dist_metr = EDIT_DISTANCE;
    else help("error: invalid option after -d");

    arg = argv[arg_idx++];
    if (arg != "-s") help("");
    scheme_str = argv[arg_idx++];
    bool is_default_scheme = scheme_str == "pigeon_hole" || scheme_str == "suffix_filter" || scheme_str == "01";
    if (!is_default_scheme) {
        if (std::filesystem::exists(scheme_str)) {
            std::string file_content;
            uint64_t file_size = std::filesystem::file_size(scheme_str);
            no_init_resize(file_content, file_size);
            std::ifstream ifile(scheme_str);
            ifile.read(file_content.data(), file_size);
            search_scheme = parse_search_scheme(file_content);
        } else help("error: invalid option after -s");
    }
    
    if (is_default_scheme) {
        arg = argv[arg_idx++];
        if (arg != "-k") help("");
        k = atoi(argv[arg_idx++]);
        if (k < 0) help("error: invalid k value");
        if      (scheme_str == "pigeon_hole")   search_scheme = pigeon_hole_scheme(k);
        else if (scheme_str == "suffix_filter") search_scheme = suffix_filter_scheme(k);
        else if (scheme_str == "01")            search_scheme = zero_one_scheme(k);
    } else {
        k = search_scheme.k;
    }

    path_index_file = argv[arg_idx];
    path_patterns_file = argv[arg_idx + 1];

    index_file.open(path_index_file);
    patterns_file.open(path_patterns_file);

    if (!index_file.good()) help("error: could not read <index_file>");
    if (!patterns_file.good()) help("error: could not read <patterns_file>");

    if (output_occurrences) {
        output_file.open(path_output_file);
        if (!output_file.good()) help("error: could not create <output_file>");
    }

    bool is_64_bit;
    index_file.read((char*) &is_64_bit, 1);
    move_r_support _support;
    index_file.read((char*) &_support, sizeof(move_r_support));
    index_file.seekg(0, std::ios::beg);

    if (_support == _count) {
        help("error: this index does not support locate");
    } else if (_support == _locate_move) {
        if (is_64_bit) {
            measure_locate<uint64_t, _locate_move>();
        } else {
            measure_locate<uint32_t, _locate_move>();
        }
    } else if (_support == _locate_rlzsa) {
        if (is_64_bit) {
            measure_locate<uint64_t, _locate_rlzsa>();
        } else {
            measure_locate<uint32_t, _locate_rlzsa>();
        }
    }
}