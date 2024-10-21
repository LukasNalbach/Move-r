#pragma once

#include "../../external/r-index-rlzsa/internal/r_index.hpp" // r-index-rlzsa
#include <BitsUtil.cpp>
#include <DynRleForRlbwt.hpp>
#include <DynSuccForRindex.hpp>
#include <OnlineRindex.hpp> // online-rlbwt
#include <r_index.hpp> // r-index
#include <r_index_f.hpp> // r-index-f
#include <rindex_types.hpp> // rcomp-glfig

// r-index
using r_index = ri::r_index<>;

// r-index-rlzsa
using r_index_rlzsa = ri_rlzsa::r_index<>;

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
void build_index<move_r<_locate_move, char, uint32_t>>(move_r<_locate_move, char, uint32_t>& index, uint16_t num_threads)
{
    index = move_r<_locate_move, char, uint32_t>(input, { .mode = _bigbwt, .num_threads = num_threads });
}

template <>
void build_index<move_r<_locate_move, char, uint64_t>>(move_r<_locate_move, char, uint64_t>& index, uint16_t num_threads)
{
    index = move_r<_locate_move, char, uint64_t>(input, { .mode = _bigbwt, .num_threads = num_threads });
}

template <>
uint64_t build_index_sa_and_bwt<int32_t, move_r<_locate_move, char, uint32_t>>(move_r<_locate_move, char, uint32_t>& index, uint16_t num_threads)
{
    index = move_r<_locate_move, char, uint32_t>(get_sa<int32_t>(), BWT, { .num_threads = num_threads });
    return 0;
}

template <>
uint64_t build_index_sa_and_bwt<int64_t, move_r<_locate_move, char, uint32_t>>(move_r<_locate_move, char, uint32_t>& index, uint16_t num_threads)
{
    index = move_r<_locate_move, char, uint32_t>(get_sa<int64_t>(), BWT, { .num_threads = num_threads });
    return 0;
}

template <>
uint64_t build_index_sa_and_bwt<int64_t, move_r<_locate_move, char, uint64_t>>(move_r<_locate_move, char, uint64_t>& index, uint16_t num_threads)
{
    index = move_r<_locate_move, char, uint64_t>(get_sa<int64_t>(), BWT, { .num_threads = num_threads });
    return 0;
}

template <>
void destroy_index<move_r<_locate_move, char, uint32_t>>(move_r<_locate_move, char, uint32_t>& index)
{
    index = move_r<_locate_move, char, uint32_t>();
}

template <>
void destroy_index<move_r<_locate_move, char, uint64_t>>(move_r<_locate_move, char, uint64_t>& index)
{
    index = move_r<_locate_move, char, uint64_t>();
}

template <>
uint32_t count_pattern<uint32_t, move_r<_locate_move, char, uint32_t>>(move_r<_locate_move, char, uint32_t>& index, std::string& pattern)
{
    return index.count(pattern);
}

template <>
uint64_t count_pattern<uint64_t, move_r<_locate_move, char, uint64_t>>(move_r<_locate_move, char, uint64_t>& index, std::string& pattern)
{
    return index.count(pattern);
}

template <>
void locate_pattern<uint32_t, move_r<_locate_move, char, uint32_t>>(move_r<_locate_move, char, uint32_t>& index, std::string& pattern, std::vector<uint32_t>& occurrences)
{
    index.locate(pattern, occurrences);
}

template <>
void locate_pattern<uint64_t, move_r<_locate_move, char, uint64_t>>(move_r<_locate_move, char, uint64_t>& index, std::string& pattern, std::vector<uint64_t>& occurrences)
{
    index.locate(pattern, occurrences);
}

// ############################# move-rlzdsa #############################

template <>
void build_index<move_r<_locate_rlzdsa, char, uint32_t>>(move_r<_locate_rlzdsa, char, uint32_t>& index, uint16_t num_threads)
{
    index = move_r<_locate_rlzdsa, char, uint32_t>(input, { .mode = _bigbwt, .num_threads = num_threads });
}

template <>
void build_index<move_r<_locate_rlzdsa, char, uint64_t>>(move_r<_locate_rlzdsa, char, uint64_t>& index, uint16_t num_threads)
{
    index = move_r<_locate_rlzdsa, char, uint64_t>(input, { .mode = _bigbwt, .num_threads = num_threads });
}

template <>
uint64_t build_index_sa_and_bwt<int32_t, move_r<_locate_rlzdsa, char, uint32_t>>(move_r<_locate_rlzdsa, char, uint32_t>& index, uint16_t num_threads)
{
    index = move_r<_locate_rlzdsa, char, uint32_t>(get_sa<int32_t>(), BWT, { .num_threads = num_threads });
    return 0;
}

template <>
uint64_t build_index_sa_and_bwt<int64_t, move_r<_locate_rlzdsa, char, uint32_t>>(move_r<_locate_rlzdsa, char, uint32_t>& index, uint16_t num_threads)
{
    index = move_r<_locate_rlzdsa, char, uint32_t>(get_sa<int64_t>(), BWT, { .num_threads = num_threads });
    return 0;
}

template <>
uint64_t build_index_sa_and_bwt<int64_t, move_r<_locate_rlzdsa, char, uint64_t>>(move_r<_locate_rlzdsa, char, uint64_t>& index, uint16_t num_threads)
{
    index = move_r<_locate_rlzdsa, char, uint64_t>(get_sa<int64_t>(), BWT, { .num_threads = num_threads });
    return 0;
}

template <>
void destroy_index<move_r<_locate_rlzdsa, char, uint32_t>>(move_r<_locate_rlzdsa, char, uint32_t>& index)
{
    index = move_r<_locate_rlzdsa, char, uint32_t>();
}

template <>
void destroy_index<move_r<_locate_rlzdsa, char, uint64_t>>(move_r<_locate_rlzdsa, char, uint64_t>& index)
{
    index = move_r<_locate_rlzdsa, char, uint64_t>();
}

template <>
uint32_t count_pattern<uint32_t, move_r<_locate_rlzdsa, char, uint32_t>>(move_r<_locate_rlzdsa, char, uint32_t>& index, std::string& pattern)
{
    return index.count(pattern);
}

template <>
uint64_t count_pattern<uint64_t, move_r<_locate_rlzdsa, char, uint64_t>>(move_r<_locate_rlzdsa, char, uint64_t>& index, std::string& pattern)
{
    return index.count(pattern);
}

template <>
void locate_pattern<uint32_t, move_r<_locate_rlzdsa, char, uint32_t>>(move_r<_locate_rlzdsa, char, uint32_t>& index, std::string& pattern, std::vector<uint32_t>& occurrences)
{
    index.locate(pattern, occurrences);
}

template <>
void locate_pattern<uint64_t, move_r<_locate_rlzdsa, char, uint64_t>>(move_r<_locate_rlzdsa, char, uint64_t>& index, std::string& pattern, std::vector<uint64_t>& occurrences)
{
    index.locate(pattern, occurrences);
}

// ############################# rcomp-glfig #############################

template <>
void build_index<rcomp_glfig>(rcomp_glfig& index, uint16_t)
{
    for (uint64_t pos = 1; pos < input_size; pos++) {
        index.extend(input[input_size - 1 - pos]);
    }
}

template <>
void destroy_index<rcomp_glfig>(rcomp_glfig&) { }

template <>
uint32_t count_pattern<uint32_t, rcomp_glfig>(rcomp_glfig& index, std::string& pattern)
{
    std::string reversed_pattern(pattern.rbegin(), pattern.rend());
    return index.count(rcomp::make_range(reversed_pattern));
}

template <>
uint64_t count_pattern<uint64_t, rcomp_glfig>(rcomp_glfig& index, std::string& pattern)
{
    std::string reversed_pattern(pattern.rbegin(), pattern.rend());
    return index.count(rcomp::make_range(reversed_pattern));
}

template <>
void locate_pattern<uint32_t, rcomp_glfig>(rcomp_glfig& index, std::string& pattern, std::vector<uint32_t>& occurrences)
{
    std::string reversed_pattern(pattern.rbegin(), pattern.rend());
    index.locate(rcomp::make_range(reversed_pattern), [&occurrences](uint32_t occurrence) { occurrences.emplace_back(input_size - 1 - occurrence); });
}

template <>
void locate_pattern<uint64_t, rcomp_glfig>(rcomp_glfig& index, std::string& pattern, std::vector<uint64_t>& occurrences)
{
    std::string reversed_pattern(pattern.rbegin(), pattern.rend());
    index.locate(rcomp::make_range(reversed_pattern), [&occurrences](uint64_t occurrence) { occurrences.emplace_back(input_size - 1 - occurrence); });
}

// ############################# r-index #############################

template <>
void build_index<r_index>(r_index& index, uint16_t num_threads)
{
    std::streambuf* cout_rdbuf = cout.rdbuf();
    std::cout.rdbuf(NULL);
    std::string name_text_file = "r-index-" + random_alphanumeric_string(10);
    std::ofstream text_file(name_text_file);
    write_to_file(text_file, input.c_str(), input_size - 1);
    text_file.close();
    index = r_index(name_text_file, "bigbwt", num_threads);
    std::filesystem::remove(name_text_file);
    std::ifstream log_file(name_text_file + ".log");
    update_peak_memory_usage(log_file);
    log_file.close();
    std::filesystem::remove(name_text_file + ".log");
    system("rm -f nul");
    std::cout.rdbuf(cout_rdbuf);
}

template <>
void destroy_index<r_index>(r_index& index)
{
    index = r_index();
}

template <>
uint32_t count_pattern<uint32_t, r_index>(r_index& index, std::string& pattern)
{
    auto range = index.count(pattern);
    return range.second - range.first + 1;
}

template <>
uint64_t count_pattern<uint64_t, r_index>(r_index& index, std::string& pattern)
{
    auto range = index.count(pattern);
    return range.second - range.first + 1;
}

template <>
void locate_pattern<uint32_t, r_index>(r_index& index, std::string& pattern, std::vector<uint32_t>& occurrences)
{
    index.locate_all<uint32_t>(pattern, occurrences);
}

template <>
void locate_pattern<uint64_t, r_index>(r_index& index, std::string& pattern, std::vector<uint64_t>& occurrences)
{
    index.locate_all<uint64_t>(pattern, occurrences);
}

// ############################# r-index-rlzsa #############################

template <>
void build_index<r_index_rlzsa>(r_index_rlzsa& index, uint16_t num_threads)
{
    index = r_index_rlzsa(input);
}

template <>
void destroy_index<r_index_rlzsa>(r_index_rlzsa& index)
{
    index = r_index_rlzsa();
}

template <>
uint32_t count_pattern<uint32_t, r_index_rlzsa>(r_index_rlzsa& index, std::string& pattern)
{
    auto range = index.count(pattern);
    return range.second - range.first + 1;
}

template <>
uint64_t count_pattern<uint64_t, r_index_rlzsa>(r_index_rlzsa& index, std::string& pattern)
{
    auto range = index.count(pattern);
    return range.second - range.first + 1;
}

template <>
void locate_pattern<uint32_t, r_index_rlzsa>(r_index_rlzsa& index, std::string& pattern, std::vector<uint32_t>& occurrences)
{
    index.locate_all<uint32_t>(pattern, occurrences);
}

template <>
void locate_pattern<uint64_t, r_index_rlzsa>(r_index_rlzsa& index, std::string& pattern, std::vector<uint64_t>& occurrences)
{
    index.locate_all<uint64_t>(pattern, occurrences);
}

// ############################# online-rlbwt #############################

template <>
void build_index<online_rlbwt>(online_rlbwt& index, uint16_t)
{
    for (uint64_t pos = 0; pos < input_size - 1; pos++) {
        index.extend(input[pos]);
    }
}

template <>
void destroy_index<online_rlbwt>(online_rlbwt&) { }

template <>
uint32_t count_pattern<uint32_t, online_rlbwt>(online_rlbwt& index, std::string& pattern)
{
    online_rlbwt::PatTracker tracker = index.getInitialPatTracker();

    for (uint32_t pos = 0; pos < pattern.size(); pos++) {
        index.lfMap(tracker, pattern[pos]);
    }

    return index.getNumOcc(tracker);
}

template <>
uint64_t count_pattern<uint64_t, online_rlbwt>(online_rlbwt& index, std::string& pattern)
{
    online_rlbwt::PatTracker tracker = index.getInitialPatTracker();

    for (uint64_t pos = 0; pos < pattern.size(); pos++) {
        index.lfMap(tracker, pattern[pos]);
    }

    return index.getNumOcc(tracker);
}

template <>
void locate_pattern<uint32_t, online_rlbwt>(online_rlbwt& index, std::string& pattern, std::vector<uint32_t>& occurrences)
{
    online_rlbwt::PatTracker tracker = index.getInitialPatTracker();

    for (uint32_t pos = 0; pos < pattern.size(); pos++) {
        index.lfMap(tracker, pattern[pos]);
    }

    uint64_t numOcc = index.getNumOcc(tracker);
    uint64_t curPos = index.calcFstOcc(tracker);
    occurrences.emplace_back(curPos - pattern.size());

    for (uint32_t pos = 1; pos < numOcc; pos++) {
        curPos = index.calcNextPos(curPos);
        occurrences.emplace_back(curPos - pattern.size());
    }
}

template <>
void locate_pattern<uint64_t, online_rlbwt>(online_rlbwt& index, std::string& pattern, std::vector<uint64_t>& occurrences)
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
void build_index<r_index_f<>>(r_index_f<>& index, uint16_t)
{
    std::streambuf* cout_rdbuf = cout.rdbuf();
    std::cout.rdbuf(NULL);
    std::string name_text_file = "r-index-f-" + random_alphanumeric_string(10);
    std::ofstream text_file(name_text_file);
    write_to_file(text_file, input.c_str(), input_size - 1);
    text_file.close();
    system(("newscanNT.x " + name_text_file + " >nul 2>nul").c_str());
    system(("build/external/pfp-thresholds/pfp-thresholds " + name_text_file +
        " -r > " + name_text_file + ".log 2>" + name_text_file + ".log").c_str());
    index = r_index_f<>(name_text_file);
    std::ifstream log_file(name_text_file + ".log");
    update_peak_memory_usage(log_file);
    log_file.close();
    system(("rm -f " + name_text_file + "* >nul 2>nul").c_str());
    system("rm -f nul");
    std::cout.rdbuf(cout_rdbuf);
}

template <>
void destroy_index<r_index_f<>>(r_index_f<>& index)
{
    index = r_index_f<>();
}

template <>
uint32_t count_pattern<uint32_t, r_index_f<>>(r_index_f<>& index, std::string& pattern)
{
    return index.count(pattern);
}

template <>
uint64_t count_pattern<uint64_t, r_index_f<>>(r_index_f<>& index, std::string& pattern)
{
    return index.count(pattern);
}

// ############################# block-rlbwt #############################

struct block_rlbwt_data {
    uint64_t time_build;
    uint64_t peak_memory_usage;
    std::string prefix_tmp_files;
    uint64_t pattern_length;
};

uint64_t write_blockrlbwt_patterns_file(std::ifstream& patterns_input_file, std::ofstream& patterns_output_file)
{
    patterns_input_file.seekg(0);
    std::string header;
    std::getline(patterns_input_file, header);
    uint64_t num_queries = number_of_patterns(header);
    uint64_t pattern_length = patterns_length(header);
    std::string pattern;
    no_init_resize(pattern, pattern_length);

    for (uint64_t cur_query = 0; cur_query < num_queries; cur_query++) {
        patterns_input_file.read((char*)&pattern[0], pattern_length);
        if (chars_remapped)
            map_string(pattern);
        patterns_output_file.write((char*)&pattern[0], pattern_length);
        patterns_output_file << std::endl;
    }

    return pattern_length;
}

block_rlbwt_data prepare_blockrlbwt()
{
    std::cout << "############## running grlBWT " << " ##############" << std::flush;

    uint64_t time_build = 0;
    uint64_t peak_memory_usage = 0;
    std::string prefix_tmp_files = random_alphanumeric_string(10);
    external_peak_memory_usage = 0;

    std::ofstream grlbwt_input_file(prefix_tmp_files + ".grlbwt");
    grlbwt_input_file.write(&input[0], input_size);
    grlbwt_input_file.close();

    auto t1 = now();
    system(("build/external/grlBWT/grlbwt-cli -T . " + prefix_tmp_files +
        ".grlbwt -o " + prefix_tmp_files + " >log_1 2>log_2").c_str());
    time_build += time_diff_ns(t1, now());
    std::ifstream log_file("log_2");
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

    std::ofstream blockrlbwt_patterns_file((prefix_tmp_files + "-blockrlbwt-patterns-bal").c_str());
    uint64_t pattern_length = write_blockrlbwt_patterns_file(patterns_file_1, blockrlbwt_patterns_file);
    blockrlbwt_patterns_file.close();

    std::cout << std::endl
              << std::endl;

    peak_memory_usage = external_peak_memory_usage;
    external_peak_memory_usage = 0;
    return { time_build, peak_memory_usage, prefix_tmp_files, pattern_length };
}

void cleanup_grlbwt(block_rlbwt_data& bd)
{
    std::filesystem::remove(bd.prefix_tmp_files + ".syms");
    std::filesystem::remove(bd.prefix_tmp_files + ".len");
    std::filesystem::remove(bd.prefix_tmp_files + "-blockrlbwt-patterns-bal");
}

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

    std::cout << "############## benchmarking " << index_name << " ##############" << std::endl
              << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(1));
    std::cout << "building " << index_name << std::flush;

    auto t1 = now();
    system(("build/external/block_RLBWT/make_bwt -h " + bd.prefix_tmp_files + ".syms -r " +
        bd.prefix_tmp_files + ".len " + blockrlbwt_param + " " + bd.prefix_tmp_files + ".rlbwt >log_1 2>log_2").c_str());
    time_build += time_diff_ns(t1, now());
    size_index = std::filesystem::file_size(bd.prefix_tmp_files + ".rlbwt") + std::filesystem::file_size(bd.prefix_tmp_files + "_data.rlbwt");
    std::ifstream log_file("log_2");
    update_peak_memory_usage(log_file);
    log_file.close();
    std::filesystem::remove("log_1");
    std::filesystem::remove("log_2");

    std::cout << std::endl;
    std::cout << "build time: " << format_time(time_build) << std::endl;
    std::cout << "peak memory usage: " << format_size(std::max(bd.peak_memory_usage, external_peak_memory_usage)) << std::endl;
    std::cout << "index size: " << format_size(size_index) << std::endl
              << std::endl;

    std::this_thread::sleep_for(std::chrono::seconds(1));
    std::cout << "counting the first set of patterns" << std::flush;

    system(("build/external/block_RLBWT/count_matches " + blockrlbwt_param + " " +
        bd.prefix_tmp_files + ".rlbwt " + bd.prefix_tmp_files + "-blockrlbwt-patterns-bal >log_1 2>log_2").c_str());
    std::ifstream results_file("blockrlbwt_count_result.txt");
    results_file.read((char*)&time_query, sizeof(uint64_t));
    results_file.read((char*)&num_queries, sizeof(uint64_t));
    results_file.read((char*)&num_occurrences, sizeof(uint64_t));
    results_file.close();

    std::filesystem::remove(bd.prefix_tmp_files + ".rlbwt");
    std::filesystem::remove(bd.prefix_tmp_files + "_data.rlbwt");
    std::filesystem::remove("blockrlbwt_count_result.txt");
    std::filesystem::remove("log_1");
    std::filesystem::remove("log_2");

    std::cout << ": " << format_query_throughput(num_queries, time_query);
    std::cout << std::endl
              << "total number of occurrences: " << num_occurrences << std::endl
              << std::endl;

    write_measurement_data(
        true, false, index_name, size_index,
        build_result { 1, time_build, std::max(bd.peak_memory_usage, external_peak_memory_usage), size_index },
        query_result { num_queries, bd.pattern_length, num_occurrences, time_query });
}