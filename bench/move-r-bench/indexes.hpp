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

#ifndef MOVE_R_BENCH_MOVE_R_BENCH_INDEXES_HPP
#define MOVE_R_BENCH_MOVE_R_BENCH_INDEXES_HPP

#include "framework.hpp"
#include <misc/index_cli.hpp>
#include <misc/pattern_reader.hpp>
#include <BitsUtil.cpp>
#include <DynRleForRlbwt.hpp>
#include <DynSuccForRindex.hpp>
#include <OnlineRindex.hpp> // online-rlbwt
#include <r_index.hpp> // r-index
#include <r_index_f.hpp> // r-index-f
#include <rindex_types.hpp> // rcomp-glfig

// the number of random characters in the names of the temporary files of the compared tools
inline constexpr uint64_t temp_name_length = 10;

// the configurations of the compared tools below are part of the measurement setup

// r-index
using r_index = ri::r_index<>;

// rcomp-glfig
using rcomp_glfig = rcomp::rindex_types::glfig_serialized<16>::type;

// online-rlbwt
using BTreeNodeT = itmmti::BTreeNode<16>; // BTree arity = {16, 32, 64, 128}
using BtmNodeMT = itmmti::BtmNodeM_StepCode<BTreeNodeT, 32>; // BtmNode arity in {16, 32, 64, 128}.
using BtmMInfoT = itmmti::BtmMInfo_BlockVec<BtmNodeMT, 512>; // Each block has 512 btmNodeM.
using BtmNodeST = itmmti::BtmNodeS<BTreeNodeT, uint32_t, 8>; // CharT = uint32_t. BtmNode arity = {4, 8, 16, 32, 64, 128}.
using BtmSInfoT = itmmti::BtmSInfo_BlockVec<BtmNodeST, 1024>; // Each block has 1024 btmNodeS.
using DynRleT = itmmti::DynRleForRlbwt<itmmti::WBitsBlockVec<1024>, itmmti::Samples_WBitsBlockVec<1024>, BtmMInfoT, BtmSInfoT>;
using BtmNodeInSucc = itmmti::BtmNodeForPSumWithVal<16>; // BtmNode arity = {16, 32, 64, 128}.
using DynSuccT = itmmti::DynSuccForRindex<BTreeNodeT, BtmNodeInSucc>;
using online_rlbwt = itmmti::OnlineRlbwtIndex<DynRleT, DynSuccT>;

// ############################# move-r #############################

template <>
void build_index<move_r<_locate_move>>(move_r<_locate_move>& index, uint64_t num_threads)
{
    index = move_r<_locate_move>(input, { .mode = _pfp, .num_threads = num_threads });
}

template <>
uint64_t count_pattern<move_r<_locate_move>>(move_r<_locate_move>& index, std::string& pattern)
{
    return index.count(pattern);
}

template <>
void locate_pattern<move_r<_locate_move>>(move_r<_locate_move>& index, std::string& pattern, std::vector<uint64_t>& occurrences)
{
    std::vector<uint48_t> occ = index.locate(pattern);
    occurrences.assign(occ.begin(), occ.end());
}

// ############################# move-rlzsa #############################

template <>
void build_index<move_r<_locate_rlzsa>>(move_r<_locate_rlzsa>& index, uint64_t num_threads)
{
    index = move_r<_locate_rlzsa>(input, { .mode = _pfp, .num_threads = num_threads });
}

template <>
uint64_t count_pattern<move_r<_locate_rlzsa>>(move_r<_locate_rlzsa>& index, std::string& pattern)
{
    return index.count(pattern);
}

template <>
void locate_pattern<move_r<_locate_rlzsa>>(move_r<_locate_rlzsa>& index, std::string& pattern, std::vector<uint64_t>& occurrences)
{
    std::vector<uint48_t> occ = index.locate(pattern);
    occurrences.assign(occ.begin(), occ.end());
}

// ############################# rcomp-glfig #############################

template <>
void build_index<rcomp_glfig>(rcomp_glfig& index, uint64_t)
{
    for (uint64_t pos = 1; pos < input_size; pos++) {
        index.extend(input[input_size - 1 - pos]);
    }
}

template <>
uint64_t count_pattern<rcomp_glfig>(rcomp_glfig& index, std::string& pattern)
{
    std::string reversed_pattern(pattern.rbegin(), pattern.rend());
    return index.count(rcomp::make_range(reversed_pattern));
}

template <>
void locate_pattern<rcomp_glfig>(rcomp_glfig& index, std::string& pattern, std::vector<uint64_t>& occurrences)
{
    std::string reversed_pattern(pattern.rbegin(), pattern.rend());
    index.locate(rcomp::make_range(reversed_pattern), [&occurrences](uint64_t occurrence) { occurrences.emplace_back(input_size - 1 - occurrence); });
}

// ############################# r-index #############################

template <>
void build_index<r_index>(r_index& index, uint64_t num_threads)
{
    std::streambuf* cout_rdbuf = cout.rdbuf();
    std::cout.rdbuf(nullptr);
    std::string text_name = "r-index-" + random_alphanumeric_string(temp_name_length);
    direct_ofstream text_file(text_name);
    write_in_chunks(text_file, input.c_str(), input_size - 1);
    text_file.close();
    index = r_index(text_name, "bigbwt", num_threads);
    std::filesystem::remove(text_name);
    direct_ifstream log_file(text_name + ".log");
    update_peak_memory_usage(log_file);
    log_file.close();
    std::filesystem::remove(text_name + ".log");
    system("rm -f nul");
    std::cout.rdbuf(cout_rdbuf);
}

template <>
uint64_t count_pattern<r_index>(r_index& index, std::string& pattern)
{
    auto range = index.count(pattern);
    return range.second - range.first + 1;
}

template <>
void locate_pattern<r_index>(r_index& index, std::string& pattern, std::vector<uint64_t>& occurrences)
{
    index.locate_all<uint64_t>(pattern, occurrences);
}

// ############################# online-rlbwt #############################

template <>
void build_index<online_rlbwt>(online_rlbwt& index, uint64_t)
{
    for (uint64_t pos = 0; pos < input_size - 1; pos++) {
        index.extend(input[pos]);
    }
}

template <>
uint64_t count_pattern<online_rlbwt>(online_rlbwt& index, std::string& pattern)
{
    online_rlbwt::PatTracker tracker = index.getInitialPatTracker();

    for (uint64_t pos = 0; pos < pattern.size(); pos++) {
        index.lfMap(tracker, pattern[pos]);
    }

    return index.getNumOcc(tracker);
}

template <>
void locate_pattern<online_rlbwt>(online_rlbwt& index, std::string& pattern, std::vector<uint64_t>& occurrences)
{
    online_rlbwt::PatTracker tracker = index.getInitialPatTracker();

    for (uint64_t pos = 0; pos < pattern.size(); pos++) {
        index.lfMap(tracker, pattern[pos]);
    }

    uint64_t numOcc = index.getNumOcc(tracker);
    uint64_t curPos = index.calcFstOcc(tracker);
    occurrences.emplace_back(curPos - pattern.size());

    for (uint64_t pos = 1; pos < numOcc; pos++) {
        curPos = index.calcNextPos(curPos);
        occurrences.emplace_back(curPos - pattern.size());
    }
}

// ############################# r-index-f #############################

template <>
void build_index<r_index_f<>>(r_index_f<>& index, uint64_t)
{
    std::streambuf* cout_rdbuf = cout.rdbuf();
    std::cout.rdbuf(nullptr);
    std::string text_name = "r-index-f-" + random_alphanumeric_string(temp_name_length);
    direct_ofstream text_file(text_name);
    write_in_chunks(text_file, input.c_str(), input_size - 1);
    text_file.close();
    system(("newscanNT.x " + text_name + " >nul 2>nul").c_str());
    system(("build/external/pfp-thresholds/pfp-thresholds " + text_name +
        " -r > " + text_name + ".log 2>" + text_name + ".log").c_str());
    index = r_index_f<>(text_name);
    direct_ifstream log_file(text_name + ".log");
    update_peak_memory_usage(log_file);
    log_file.close();
    system(("rm -f " + text_name + "* >nul 2>nul").c_str());
    system("rm -f nul");
    std::cout.rdbuf(cout_rdbuf);
}

template <>
uint64_t count_pattern<r_index_f<>>(r_index_f<>& index, std::string& pattern) { return index.count(pattern); }

// ############################# block-rlbwt #############################

// the results of prepare_blockrlbwt
struct block_rlbwt_data {
    uint64_t time_build; // time in ns of grlBWT and of building block_RLBWT
    uint64_t peak_memory_usage; // peak memory usage of grlBWT in bytes
    std::string prefix_tmp_files; // prefix of the names of the temporary files
    uint64_t pattern_length; // length of the patterns of the first pattern file
};

// writes the patterns of a pattern file one per line (the input format of block_RLBWT) and returns their length
uint64_t write_blockrlbwt_patterns_file(std::istream& patterns_input_file, std::ostream& patterns_output_file)
{
    patterns_input_file.seekg(0);
    auto [num_queries, pattern_length] = read_patterns_header(patterns_input_file);
    std::string pattern;
    resize_no_init(pattern, pattern_length);

    for (uint64_t cur_query = 0; cur_query < num_queries; cur_query++) {
        patterns_input_file.read((char*) &pattern[0], pattern_length);
        if (chars_remapped) map_string(pattern);
        patterns_output_file.write((char*) &pattern[0], pattern_length);
        patterns_output_file << std::endl;
    }

    return pattern_length;
}

// computes the run-length BWT of the input with grlBWT, builds block_RLBWT for its alphabet and writes the patterns of
// the first pattern file for block_RLBWT
block_rlbwt_data prepare_blockrlbwt()
{
    std::cout << "############## running grlBWT " << " ##############" << std::flush;

    uint64_t time_build = 0;
    uint64_t peak_memory_usage = 0;
    std::string prefix_tmp_files = random_alphanumeric_string(temp_name_length);
    external_peak_memory_usage = 0;

    direct_ofstream grlbwt_input_file(prefix_tmp_files + ".grlbwt");
    grlbwt_input_file.write(&input[0], input_size);
    grlbwt_input_file.close();

    auto t1 = now();
    system(("build/external/grlBWT/grlbwt-cli -T . " + prefix_tmp_files +
        ".grlbwt -o " + prefix_tmp_files + " >log_1 2>log_2").c_str());
    time_build += time_diff_ns(t1, now());
    direct_ifstream log_file("log_2");
    update_peak_memory_usage(log_file);
    log_file.close();
    std::filesystem::remove(prefix_tmp_files + ".grlbwt");
    std::filesystem::remove("log_1");
    std::filesystem::remove("log_2");

    t1 = now();
    system(("build/external/grlBWT/grlbwt2rle " + prefix_tmp_files +
        ".rl_bwt " + prefix_tmp_files + " >log_1 2>log_2").c_str());
    time_build += time_diff_ns(t1, now());
    log_file.open("log_2");
    update_peak_memory_usage(log_file);
    log_file.close();
    std::filesystem::remove("log_1");
    std::filesystem::remove("log_2");
    std::filesystem::remove(prefix_tmp_files + ".rl_bwt");

    t1 = now();
    system(("build/external/block_RLBWT/make_alphabet_header -h " + prefix_tmp_files + ".syms -r " +
        prefix_tmp_files + ".len  > external/block_RLBWT/include/custom_alphabet.hpp 2>log_2").c_str());
    time_build += time_diff_ns(t1, now());
    std::filesystem::remove("log_2");
    t1 = now();
    system("cd external/block_RLBWT/ && make make_bwt >log_1 2>log_2 && make count_matches >log_1 2>log_2 && rm log_* && cd ../.. ");
    time_build += time_diff_ns(t1, now());
    system("mv external/block_RLBWT/make_bwt build/external/block_RLBWT/make_bwt");
    system("mv external/block_RLBWT/count_matches build/external/block_RLBWT/count_matches");

    direct_ofstream blockrlbwt_patterns_file((prefix_tmp_files + "-blockrlbwt-patterns-bal").c_str());
    uint64_t pattern_length = write_blockrlbwt_patterns_file(patterns_file_1, blockrlbwt_patterns_file);
    blockrlbwt_patterns_file.close();
    std::cout << std::endl << std::endl;

    peak_memory_usage = external_peak_memory_usage;
    external_peak_memory_usage = 0;
    return { time_build, peak_memory_usage, prefix_tmp_files, pattern_length };
}

// removes the temporary files of prepare_blockrlbwt
void cleanup_grlbwt(block_rlbwt_data& bd)
{
    std::filesystem::remove(bd.prefix_tmp_files + ".syms");
    std::filesystem::remove(bd.prefix_tmp_files + ".len");
    std::filesystem::remove(bd.prefix_tmp_files + "-blockrlbwt-patterns-bal");
}

// builds the block_RLBWT variant index_name (block_rlbwt_2, block_rlbwt_v or block_rlbwt_r), measures count with it
// and writes the RESULT lines
void measure_blockrlbwt(std::string index_name, block_rlbwt_data& bd)
{
    uint64_t time_build = bd.time_build;
    uint64_t size_index;
    uint64_t time_query;
    uint64_t num_queries;
    uint64_t num_occurrences;
    std::string blockrlbwt_param = "";

    if (index_name == "block_rlbwt_v") {
        blockrlbwt_param = "-s";
    } else if (index_name == "block_rlbwt_r") {
        blockrlbwt_param = "-c";
    }

    std::cout << "############## benchmarking " << index_name << " ##############" << std::endl << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(1));
    std::cout << "building " << index_name << std::flush;

    auto t1 = now();
    system(("build/external/block_RLBWT/make_bwt -h " + bd.prefix_tmp_files + ".syms -r " +
        bd.prefix_tmp_files + ".len " + blockrlbwt_param + " " + bd.prefix_tmp_files + ".rlbwt >log_1 2>log_2").c_str());
    time_build += time_diff_ns(t1, now());
    size_index = std::filesystem::file_size(bd.prefix_tmp_files + ".rlbwt") + std::filesystem::file_size(bd.prefix_tmp_files + "_data.rlbwt");
    direct_ifstream log_file("log_2");
    update_peak_memory_usage(log_file);
    log_file.close();
    std::filesystem::remove("log_1");
    std::filesystem::remove("log_2");

    std::cout << std::endl;
    std::cout << "build time: " << format_time(time_build) << std::endl;
    std::cout << "peak memory usage: " << format_size(std::max(bd.peak_memory_usage, external_peak_memory_usage)) << std::endl;
    std::cout << "index size: " << format_size(size_index) << std::endl << std::endl;

    std::this_thread::sleep_for(std::chrono::seconds(1));
    std::cout << "counting the first set of patterns" << std::flush;

    system(("build/external/block_RLBWT/count_matches " + blockrlbwt_param + " " +
        bd.prefix_tmp_files + ".rlbwt " + bd.prefix_tmp_files + "-blockrlbwt-patterns-bal >log_1 2>log_2").c_str());
    direct_ifstream results_file("blockrlbwt_count_result.txt");
    results_file.read((char*) &time_query, sizeof(uint64_t));
    results_file.read((char*) &num_queries, sizeof(uint64_t));
    results_file.read((char*) &num_occurrences, sizeof(uint64_t));
    results_file.close();

    std::filesystem::remove(bd.prefix_tmp_files + ".rlbwt");
    std::filesystem::remove(bd.prefix_tmp_files + "_data.rlbwt");
    std::filesystem::remove("blockrlbwt_count_result.txt");
    std::filesystem::remove("log_1");
    std::filesystem::remove("log_2");

    std::cout << ": " << format_query_throughput(num_queries, time_query);
    std::cout << std::endl << "total number of occurrences: " << num_occurrences << std::endl << std::endl;

    write_measurement_data(
        true, false, index_name, size_index,
        build_result { 1, time_build, std::max(bd.peak_memory_usage, external_peak_memory_usage), size_index },
        query_result { num_queries, bd.pattern_length, num_occurrences, time_query });
}

#endif // MOVE_R_BENCH_MOVE_R_BENCH_INDEXES_HPP
