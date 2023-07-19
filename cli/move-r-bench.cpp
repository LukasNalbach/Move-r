#include <iostream>
#include <filesystem>
#include <string_view>

#include <move_r/move_r.hpp> // move-r
#include <rindex_types.hpp> // rcomp
#include "../external/r-index-prezza/internal/r_index.hpp" // r-index-prezza
#include "../external/r-index-mun/internal/r_index.hpp" // r-index-mun
#include <OnlineRindex.hpp> // OnlineRlbwt
#include <DynRleForRlbwt.hpp>
#include <DynSuccForRindex.hpp>
#include <BitsUtil.cpp>
#include <dynamic/dynamic.hpp> // DYNAMIC
#include <r_index_f.hpp> // r-index-f

uint16_t max_num_threads;
std::string input;
std::string input_reverted;
bool check_correctness = false;
uint64_t external_peak_memory_usage;
std::ofstream measurement_file;
std::ifstream input_file;
std::ifstream patterns_file_1;
std::ifstream patterns_file_2;
std::string path_inputfile;
std::string path_patternsfile_1;
std::string path_patternsfile_2;
std::string name_textfile;
int ptr = 1;

void help(std::string msg) {
	if (msg != "") std::cout << msg << std::endl;
	std::cout << "move-r-bench: benchmarks construction-, revert- and query-performance of move-r, r-index-f, rcomp, r-index" << std::endl;
	std::cout << "              (Prezza), r-index (Mun), OnlineRLBWT and DYNAMIC; has to be executed from the base folder." << std::endl << std::endl;
	std::cout << "usage: move-r-bench [options] <input_file> <patterns_file_1> <patterns_file_2> <num_threads>" << std::endl;
	std::cout << "   -c                 check for correctnes if possible; disables the -m option; will not print" << std::endl;
	std::cout << "                      runtime data if the runtime could be affected by checking for correctness" << std::endl;
	std::cout << "   -m <m_file>        writes measurement data to m_file" << std::endl;
	std::cout << "   <input_file>       input file" << std::endl;
	std::cout << "   <patterns_file_1>  file containing patterns (pattern length ~ number of occurrences) from <input_file>" << std::endl;
	std::cout << "                      to count and locate" << std::endl;
	std::cout << "   <patterns_file_2>  file containing patterns (pattern length << number of occurrences) from <input_file>" << std::endl;
	std::cout << "                      to locate" << std::endl;
	std::cout << "   <num_threads>      maximum number of threads to use" << std::endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr) {
	std::string s = argv[ptr];
	ptr++;

	if (s == "-m") {
		if (ptr >= argc-1) help("error: missing parameter after -m option");
		std::string path_measurement_file = argv[ptr++];
		measurement_file.open(path_measurement_file,std::filesystem::exists(path_measurement_file) ? std::ios::app : std::ios::out);
		if (!measurement_file.good()) help("error: cannot open or create measurement file");
	} else if (s == "-c") {
		if (ptr >= argc-1) help("error: missing parameter after -c option");
		check_correctness = true;
	} else {
		help("error: unrecognized '" + s +  "' option");
	}
}

// ############################# benchmarking framework #############################

uint64_t peak_memory_usage(std::ifstream& log_file) {
	std::string log_file_content;
	log_file.seekg(0,std::ios::end);
	no_init_resize(log_file_content,log_file.tellg());
	log_file.seekg(0,std::ios::beg);
	log_file.read((char*)&log_file_content[0],log_file_content.size());
	int32_t pos = 0;
	uint64_t cur_peak = 0;
	std::string str_cur_peak;

	while ((pos = log_file_content.find(", peak: ",pos)) != -1) {
		pos += 8;

		while (log_file_content[pos] != ',') {
			str_cur_peak.push_back(log_file_content[pos]);
			pos++;
		}

		cur_peak = std::max(cur_peak,(uint64_t)stol(str_cur_peak));
		str_cur_peak.clear();
	}

	return cur_peak;
}

template <typename idx_t, bool alternative_build_mode>
void build_index(idx_t& index, uint16_t num_threads);

template <typename idx_t>
void destroy_index(idx_t& index);

template <typename idx_t>
void revert_index(idx_t& index, uint16_t num_threads);

template <typename uint_t, typename idx_t>
uint_t count_pattern(idx_t& index, std::string& pattern);

template <typename uint_t, typename idx_t>
void locate_pattern(idx_t& index, std::string& pattern, std::vector<uint_t>& occurrences);

struct build_result {
	uint16_t num_threads;
	uint64_t time_build;
	uint64_t peak_memory_usage;
	uint64_t index_size;
};

struct revert_result {
	uint16_t num_threads;
	uint64_t time_revert;
};

struct query_result {
	uint64_t num_queries;
	uint64_t pattern_length;
	uint64_t num_occurrences;
	uint64_t time_query;
};

template <typename uint_t, typename idx_t>
query_result count_patterns(idx_t& index) {
	patterns_file_1.seekg(0);
	std::string header;
	std::getline(patterns_file_1,header);
	uint64_t num_queries = number_of_patterns(header);
	uint64_t pattern_length = patterns_length(header);
	uint64_t num_occurrences = 0;
	uint64_t time_query = 0;
	std::string pattern;
	no_init_resize(pattern,pattern_length);
	std::chrono::steady_clock::time_point t1,t2;

	for (uint64_t cur_query=0; cur_query<num_queries; cur_query++) {
		patterns_file_1.read((char*)&pattern[0],pattern_length);
		t1 = now();
		num_occurrences += count_pattern<uint_t,idx_t>(index,pattern);
		t2 = now();
		time_query += time_diff_ns(t1,t2);
	}

	return query_result{num_queries,pattern_length,num_occurrences,time_query};
}

template <typename uint_t, typename idx_t>
query_result locate_patterns(idx_t& index, std::ifstream& patterns_file) {
	patterns_file.seekg(0);
	std::string header;
	std::getline(patterns_file,header);
	uint64_t num_queries = number_of_patterns(header);
	uint64_t pattern_length = patterns_length(header);
	uint64_t num_occurrences = 0;
	std::vector<uint_t> occurrences;
	uint64_t time_query = 0;
	std::string pattern;
	no_init_resize(pattern,pattern_length);
	std::chrono::steady_clock::time_point t1,t2;
	bool correct = true;

	for (uint64_t cur_query=0; cur_query<num_queries; cur_query++) {
		patterns_file.read((char*)&pattern[0],pattern_length);
		t1 = now();
		locate_pattern<uint_t,idx_t>(index,pattern,occurrences);
		t2 = now();
		time_query += time_diff_ns(t1,t2);
		num_occurrences += occurrences.size();

		if (check_correctness) {
			for (uint_t occurrence : occurrences) {
				for (uint_t pos=0; pos<pattern_length; pos++) {
					if (input[occurrence+pos] != pattern[pos]) {
						correct = false;
						break;
					}
				}

				if (!correct) break;
			}

			if (!correct || occurrences.size() != count_pattern<uint_t,idx_t>(index,pattern)) {
				correct = false;
				break;
			}
		}

		occurrences.clear();
	}

	if (check_correctness) {
		if (correct) std::cout << " (no wrong occurrences)";
		else std::cout << " (wrong occurrences)";
	}

	return query_result{num_queries,pattern_length,num_occurrences,time_query};
}

template <typename uint_t, typename idx_t, bool alternative_build_mode, bool measure_revert, bool measure_count, bool measure_locate>
void measure(std::string index_name, std::string index_log_name, uint16_t max_build_threads, uint16_t max_revert_threads) {
	idx_t index;
	std::chrono::steady_clock::time_point t1,t2;
	uint64_t m1,m2;
	std::vector<build_result> results_build;
	std::vector<revert_result> results_revert;
	query_result result_count,result_locate_1,result_locate_2;
	std::cout << "############## benchmarking " << index_name << " ##############" << std::endl << std::endl;

	for (uint16_t cur_num_threads=1; cur_num_threads<=max_build_threads; cur_num_threads*=2) {
		std::cout << "building " << index_name << " using " << format_threads(cur_num_threads) << std::flush;

		destroy_index<idx_t>(index);
		malloc_count_reset_peak();
		external_peak_memory_usage = 0;
		m1 = malloc_count_current();
		t1 = now();
		build_index<idx_t,alternative_build_mode>(index,cur_num_threads);
		t2 = now();
		m2 = malloc_count_current();
		
		build_result res{
			.num_threads = cur_num_threads,
			.time_build = time_diff_ns(t1,t2),
			.peak_memory_usage = std::max(malloc_count_peak()-m1,external_peak_memory_usage),
			.index_size = m2-m1
		};

		results_build.emplace_back(res);
		std::cout << std::endl;
		std::cout << "build time: " << format_time(res.time_build) << std::endl;
		std::cout << "peak memory usage: " << format_size(res.peak_memory_usage) << std::endl;
		std::cout << "index size: " << format_size(res.index_size) << std::endl << std::endl;
	}

	if constexpr (measure_revert) {
		no_init_resize(input_reverted,input.size());

		for (uint16_t cur_num_threads=1; cur_num_threads<=max_revert_threads; cur_num_threads*=2) {
			std::cout << "reverting the index using " << format_threads(cur_num_threads) << std::flush;

			std::memset(&input_reverted[0],0,input_reverted.size());
			t1 = now();
			revert_index<idx_t>(index,cur_num_threads);
			t2 = now();

			revert_result res{cur_num_threads,time_diff_ns(t1,t2)};
			results_revert.emplace_back(res);
			std::cout << std::endl << "revert time: " << format_time(res.time_revert) << std::endl;

			if (check_correctness) {
				std::cout << "checking correctness of the reverted text:" << std::flush;
				bool correct = true;

				for (uint_t pos=0; pos<input.size(); pos++) {
					if (input[pos] != input_reverted[pos]) {
						correct = false;
						break;
					}
				}

				if (correct) std::cout << " correct" << std::endl;
				else std::cout << " not correct" << std::endl;
			}
			
			std::cout << std::endl;
		}

		input_reverted.clear();
		input_reverted.shrink_to_fit();
	}
	
	if constexpr (measure_count) {
		std::cout << "counting the first set of patterns" << std::flush;
		result_count = count_patterns<uint_t>(index);
		if (!check_correctness) std::cout << ": " << format_query_throughput(result_count.num_queries,result_count.time_query);
		std::cout << std::endl << "total number of occurences: " << result_count.num_occurrences << std::endl;
	}

	if constexpr (measure_locate) {
		std::cout << "locating the first set of patterns" << std::flush;
		result_locate_1 = locate_patterns<uint_t,idx_t>(index,patterns_file_1);
		if (!check_correctness) std::cout << ": " << format_query_throughput(result_locate_1.num_queries,result_locate_1.time_query);
		std::cout << std::endl << "total number of occurences: " << result_locate_1.num_occurrences << std::endl;
		
		std::cout << "locating the second set of patterns" << std::flush;
		result_locate_2 = locate_patterns<uint_t,idx_t>(index,patterns_file_2);
		if (!check_correctness) std::cout << ": " << format_query_throughput(result_locate_2.num_queries,result_locate_2.time_query);
		std::cout << std::endl << "total number of occurences: " << result_locate_2.num_occurrences << std::endl;
	}

	if constexpr (measure_count || measure_locate) std::cout << std::endl;

	if (measurement_file.is_open()) {
		uint64_t index_size = results_build.front().index_size;

		for (build_result& res : results_build) {
			measurement_file << "RESULT"
				<< " type=comparison_build"
				<< " implementation=" << index_log_name
				<< " text=" << name_textfile
				<< " num_threads=" << res.num_threads
				<< " time_build=" << res.time_build
				<< " peak_memory_usage=" << res.peak_memory_usage
				<< " index_size=" << res.index_size
				<< std::endl;
		}

		for (revert_result& res : results_revert) {
			measurement_file << "RESULT"
				<< " type=comparison_revert"
				<< " implementation=" << index_log_name
				<< " text=" << name_textfile
				<< " num_threads=" << res.num_threads
				<< " time_revert=" << res.time_revert
				<< " index_size=" << index_size
				<< std::endl;
		}

		if constexpr (measure_count) {
			measurement_file << "RESULT"
				<< " type=comparison_count"
				<< " implementation=" << index_log_name
				<< " text=" << name_textfile
				<< " num_queries=" << result_count.num_queries
				<< " pattern_length=" << result_count.pattern_length
				<< " num_occurrences=" << result_count.num_occurrences
				<< " time_query=" << result_count.time_query
				<< " index_size=" << index_size
				<< std::endl;
		}

		if constexpr (measure_locate) {
			std::vector<query_result> locate_results = {result_locate_1,result_locate_2};

			for (query_result& res : locate_results) {
				measurement_file << "RESULT"
					<< " type=comparison_locate"
					<< " implementation=" << index_log_name
					<< " text=" << name_textfile
					<< " num_queries=" << res.num_queries
					<< " pattern_length=" << res.pattern_length
					<< " num_occurrences=" << result_count.num_occurrences
					<< " time_query=" << res.time_query
					<< " index_size=" << results_build.back().index_size
					<< std::endl;
			}
		}
	}
}

// ############################# move-r #############################

template <>
void build_index<move_r<uint32_t>,false>(move_r<uint32_t>& index, uint16_t num_threads) {
	index = std::move(move_r<uint32_t>(input,full_support,runtime,num_threads));
}

template <>
void build_index<move_r<uint64_t>,false>(move_r<uint64_t>& index, uint16_t num_threads) {
	index = std::move(move_r<uint64_t>(input,full_support,runtime,num_threads));
}

template <>
void build_index<move_r<uint32_t>,true>(move_r<uint32_t>& index, uint16_t num_threads) {
	index = std::move(move_r<uint32_t>(input,full_support,space,num_threads));
}

template <>
void build_index<move_r<uint64_t>,true>(move_r<uint64_t>& index, uint16_t num_threads) {
	index = std::move(move_r<uint64_t>(input,full_support,space,num_threads));
}

template <>
void destroy_index<move_r<uint32_t>>(move_r<uint32_t>& index) {
	index = std::move(move_r<uint32_t>());
}

template <>
void destroy_index<move_r<uint64_t>>(move_r<uint64_t>& index) {
	index = std::move(move_r<uint64_t>());
}

template <>
void revert_index<move_r<uint32_t>>(move_r<uint32_t>& index, uint16_t num_threads) {
	index.revert_range(input_reverted,0,index.input_size()-1,0,index.input_size()-1,num_threads);
}

template <>
void revert_index<move_r<uint64_t>>(move_r<uint64_t>& index, uint16_t num_threads) {
	index.revert_range(input_reverted,0,index.input_size()-1,0,index.input_size()-1,num_threads);
}

template <>
uint32_t count_pattern<uint32_t,move_r<uint32_t>>(move_r<uint32_t>& index, std::string& pattern) {
	return index.count(pattern);
}

template <>
uint64_t count_pattern<uint64_t,move_r<uint64_t>>(move_r<uint64_t>& index, std::string& pattern) {
	return index.count(pattern);
}

template <>
void locate_pattern<uint32_t,move_r<uint32_t>>(move_r<uint32_t>& index, std::string& pattern, std::vector<uint32_t>& occurrences) {
	index.locate(pattern,occurrences);
}

template <>
void locate_pattern<uint64_t,move_r<uint64_t>>(move_r<uint64_t>& index, std::string& pattern, std::vector<uint64_t>& occurrences) {
	index.locate(pattern,occurrences);
}

// ############################# rcomp #############################

using rcomp_lfig = rcomp::rindex_types::lfig_naive<7>::type;
using rcomp_glfig_16 = rcomp::rindex_types::glfig_serialized<16>::type;

template <>
void build_index<rcomp_lfig,false>(rcomp_lfig& index, uint16_t) {
	for (size_t pos=1; pos<=input.size(); pos++) {
        index.extend(input[input.size()-pos]);
    }
}

template <>
void build_index<rcomp_glfig_16,false>(rcomp_glfig_16& index, uint16_t) {
	for (size_t pos=1; pos<=input.size(); pos++) {
        index.extend(input[input.size()-pos]);
    }
}

template <>
void destroy_index<rcomp_lfig>(rcomp_lfig&) {}

template <>
void destroy_index<rcomp_glfig_16>(rcomp_glfig_16&) {}

template <>
void revert_index<rcomp_lfig>(rcomp_lfig& index, uint16_t) {
	uint64_t pos = input_reverted.size()-1;
    index.decode_text([&pos](char c){input_reverted[pos--] = c;});
}

template <>
void revert_index<rcomp_glfig_16>(rcomp_glfig_16& index, uint16_t) {
	uint64_t pos = input_reverted.size()-1;
    index.decode_text([&pos](char c){input_reverted[pos--] = c;});
}

template <>
uint32_t count_pattern<uint32_t,rcomp_lfig>(rcomp_lfig& index, std::string& pattern) {
	std::string reversed_pattern(pattern.rbegin(),pattern.rend());
	return index.count(rcomp::make_range(reversed_pattern));
}

template <>
uint64_t count_pattern<uint64_t,rcomp_lfig>(rcomp_lfig& index, std::string& pattern) {
	std::string reversed_pattern(pattern.rbegin(),pattern.rend());
	return index.count(rcomp::make_range(reversed_pattern));
}

template <>
uint32_t count_pattern<uint32_t,rcomp_glfig_16>(rcomp_glfig_16& index, std::string& pattern) {
	std::string reversed_pattern(pattern.rbegin(),pattern.rend());
	return index.count(rcomp::make_range(reversed_pattern));
}

template <>
uint64_t count_pattern<uint64_t,rcomp_glfig_16>(rcomp_glfig_16& index, std::string& pattern) {
	std::string reversed_pattern(pattern.rbegin(),pattern.rend());
	return index.count(rcomp::make_range(reversed_pattern));
}

template <>
void locate_pattern<uint32_t,rcomp_lfig>(rcomp_lfig& index, std::string& pattern, std::vector<uint32_t>& occurrences) {
	std::string reversed_pattern(pattern.rbegin(),pattern.rend());
	uint32_t input_size = index.get_num_chars();
    index.locate(rcomp::make_range(reversed_pattern),[&occurrences,&input_size](uint32_t occurrence){occurrences.emplace_back(input_size-occurrence);});
}

template <>
void locate_pattern<uint64_t,rcomp_lfig>(rcomp_lfig& index, std::string& pattern, std::vector<uint64_t>& occurrences) {
	std::string reversed_pattern(pattern.rbegin(),pattern.rend());
	uint64_t input_size = index.get_num_chars();
    index.locate(rcomp::make_range(reversed_pattern),[&occurrences,&input_size](uint64_t occurrence){occurrences.emplace_back(input_size-occurrence);});
}

template <>
void locate_pattern<uint32_t,rcomp_glfig_16>(rcomp_glfig_16& index, std::string& pattern, std::vector<uint32_t>& occurrences) {
	std::string reversed_pattern(pattern.rbegin(),pattern.rend());
	uint32_t input_size = index.get_num_chars();
    index.locate(rcomp::make_range(reversed_pattern),[&occurrences,&input_size](uint32_t occurrence){occurrences.emplace_back(input_size-occurrence);});
}

template <>
void locate_pattern<uint64_t,rcomp_glfig_16>(rcomp_glfig_16& index, std::string& pattern, std::vector<uint64_t>& occurrences) {
	std::string reversed_pattern(pattern.rbegin(),pattern.rend());
	uint64_t input_size = index.get_num_chars();
    index.locate(rcomp::make_range(reversed_pattern),[&occurrences,&input_size](uint64_t occurrence){occurrences.emplace_back(input_size-occurrence);});
}

// ############################# r-index-prezza #############################

template <>
void build_index<ri::r_index<>,false>(ri::r_index<>& index, uint16_t) {
    std::streambuf* cout_rfbuf = cout.rdbuf();
    std::cout.rdbuf(NULL);
	index = ri::r_index<>(input,true);
    std::cout.rdbuf(cout_rfbuf);
}

template <>
void destroy_index<ri::r_index<>>(ri::r_index<>& index) {
	index = std::move(ri::r_index<>());
}

template <>
void revert_index<ri::r_index<>>(ri::r_index<>& index, uint16_t) {
    index.revert(input_reverted);
}

template <>
uint32_t count_pattern<uint32_t,ri::r_index<>>(ri::r_index<>& index, std::string& pattern) {
	auto range = index.count(pattern);
	return range.second - range.first + 1;
}

template <>
uint64_t count_pattern<uint64_t,ri::r_index<>>(ri::r_index<>& index, std::string& pattern) {
	auto range = index.count(pattern);
	return range.second - range.first + 1;
}

template <>
void locate_pattern<uint32_t,ri::r_index<>>(ri::r_index<>& index, std::string& pattern, std::vector<uint32_t>& occurrences) {
    index.locate_all<uint32_t>(pattern,occurrences);
}

template <>
void locate_pattern<uint64_t,ri::r_index<>>(ri::r_index<>& index, std::string& pattern, std::vector<uint64_t>& occurrences) {
    index.locate_all<uint64_t>(pattern,occurrences);
}

// ############################# r-index-mun #############################

template <>
void build_index<ri_mun::r_index<>,false>(ri_mun::r_index<>& index, uint16_t num_threads) {
    std::streambuf* cout_rfbuf = cout.rdbuf();
    std::cout.rdbuf(NULL);
    std::string name_textfile = "r-index-" + random_alphanumeric_string(10);
	std::ofstream text_file(name_textfile);
	write_to_file(text_file,input.c_str(),input.size());
	text_file.close();
	index = ri_mun::r_index<>(name_textfile,"bigbwt",num_threads);
	std::filesystem::remove(name_textfile);
	std::ifstream log_file(name_textfile + ".log");
	external_peak_memory_usage = peak_memory_usage(log_file);
	log_file.close();
	std::filesystem::remove(name_textfile + ".log");
	system("rm -f nul");
    std::cout.rdbuf(cout_rfbuf);
}

template <>
void destroy_index<ri_mun::r_index<>>(ri_mun::r_index<>& index) {
	index = std::move(ri_mun::r_index<>());
}

template <>
void revert_index<ri_mun::r_index<>>(ri_mun::r_index<>& index, uint16_t) {
    index.revert(input_reverted);
}

template <>
uint32_t count_pattern<uint32_t,ri_mun::r_index<>>(ri_mun::r_index<>& index, std::string& pattern) {
	auto range = index.count(pattern);
	return range.second - range.first + 1;
}

template <>
uint64_t count_pattern<uint64_t,ri_mun::r_index<>>(ri_mun::r_index<>& index, std::string& pattern) {
	auto range = index.count(pattern);
	return range.second - range.first + 1;
}

template <>
void locate_pattern<uint32_t,ri_mun::r_index<>>(ri_mun::r_index<>& index, std::string& pattern, std::vector<uint32_t>& occurrences) {
    index.locate_all<uint32_t>(pattern,occurrences);
}

template <>
void locate_pattern<uint64_t,ri_mun::r_index<>>(ri_mun::r_index<>& index, std::string& pattern, std::vector<uint64_t>& occurrences) {
    index.locate_all<uint64_t>(pattern,occurrences);
}

// ############################# OnlineRlbwt #############################

using BTreeNodeT = itmmti::BTreeNode<16>; // BTree arity = {16, 32, 64, 128}
using BtmNodeMT = itmmti::BtmNodeM_StepCode<BTreeNodeT, 32>; // BtmNode arity in {16, 32, 64, 128}.
using BtmMInfoT = itmmti::BtmMInfo_BlockVec<BtmNodeMT, 512>; // Each block has 512 btmNodeM.
using BtmNodeST = itmmti::BtmNodeS<BTreeNodeT, uint32_t, 8>; // CharT = uint32_t. BtmNode arity = {4, 8, 16, 32, 64, 128}.
using BtmSInfoT = itmmti::BtmSInfo_BlockVec<BtmNodeST, 1024>; // Each block has 1024 btmNodeS.
using DynRleT = itmmti::DynRleForRlbwt<itmmti::WBitsBlockVec<1024>, itmmti::Samples_WBitsBlockVec<1024>, BtmMInfoT, BtmSInfoT>;
using BtmNodeInSucc = itmmti::BtmNodeForPSumWithVal<16>; // BtmNode arity = {16, 32, 64, 128}.
using DynSuccT = itmmti::DynSuccForRindex<BTreeNodeT, BtmNodeInSucc>;
using OnlineRlbwt = itmmti::OnlineRlbwtIndex<DynRleT, DynSuccT>;

template <>
void build_index<OnlineRlbwt,false>(OnlineRlbwt& index, uint16_t) {
	for (uint64_t pos=0; pos<input.size(); pos++) {
		index.extend(input[pos]);
	}
}

template <>
void destroy_index<OnlineRlbwt>(OnlineRlbwt&) {}

template <>
void revert_index(OnlineRlbwt& index, uint16_t) {
    index.revert(input_reverted);
}

template <>
uint32_t count_pattern<uint32_t,OnlineRlbwt>(OnlineRlbwt& index, std::string& pattern) {
    OnlineRlbwt::PatTracker tracker = index.getInitialPatTracker();

	for (uint32_t pos=0; pos<pattern.size(); pos++) {
		index.lfMap(tracker,pattern[pos]);
	}

	return index.getNumOcc(tracker);
}

template <>
uint64_t count_pattern<uint64_t,OnlineRlbwt>(OnlineRlbwt& index, std::string& pattern) {
    OnlineRlbwt::PatTracker tracker = index.getInitialPatTracker();

	for (uint64_t pos=0; pos<pattern.size(); pos++) {
		index.lfMap(tracker,pattern[pos]);
	}

	return index.getNumOcc(tracker);
}

template <>
void locate_pattern<uint32_t,OnlineRlbwt>(OnlineRlbwt& index, std::string& pattern, std::vector<uint32_t>& occurrences) {
   	OnlineRlbwt::PatTracker tracker = index.getInitialPatTracker();

	for (uint32_t pos=0; pos<pattern.size(); pos++) {
		index.lfMap(tracker,pattern[pos]);
	}

    uint64_t numOcc = index.getNumOcc(tracker);
    uint64_t curPos = index.calcFstOcc(tracker);
	occurrences.emplace_back(curPos-pattern.size());

	for (uint32_t pos=1; pos<numOcc; pos++) {
		curPos = index.calcNextPos(curPos);
      	occurrences.emplace_back(curPos-pattern.size());
    }
}

template <>
void locate_pattern<uint64_t,OnlineRlbwt>(OnlineRlbwt& index, std::string& pattern, std::vector<uint64_t>& occurrences) {
   	OnlineRlbwt::PatTracker tracker = index.getInitialPatTracker();

	for (uint64_t pos=0; pos<pattern.size(); pos++) {
		index.lfMap(tracker,pattern[pos]);
	}

    uint64_t numOcc = index.getNumOcc(tracker);
    uint64_t curPos = index.calcFstOcc(tracker);
	occurrences.emplace_back(curPos-pattern.size());

	for (uint64_t pos=1; pos<numOcc; pos++) {
		curPos = index.calcNextPos(curPos);
      	occurrences.emplace_back(curPos-pattern.size());
    }
}

// ############################# DYNAMIC #############################

using DYNAMIC = dyn::rle_bwt;

template <>
void build_index<DYNAMIC,false>(DYNAMIC& index, uint16_t) {
	for (uint64_t pos=0; pos<input.size(); pos++) {
		index.extend(input[pos]);
	}
}

template <>
void destroy_index<DYNAMIC>(DYNAMIC&) {}

template <>
void revert_index<DYNAMIC>(DYNAMIC& index, uint16_t) {
    index.revert(input_reverted);
}

template <>
uint32_t count_pattern<uint32_t,DYNAMIC>(DYNAMIC& index, std::string& pattern) {
	auto range = index.count(pattern);
	return range.second - range.first;
}

template <>
uint64_t count_pattern<uint64_t,DYNAMIC>(DYNAMIC& index, std::string& pattern) {
	auto range = index.count(pattern);
	return range.second - range.first;
}

// ############################# r-index-f #############################

template <>
void build_index<r_index_f<>,false>(r_index_f<>& index, uint16_t) {
    std::streambuf* cout_rfbuf = cout.rdbuf();
	std::cout.rdbuf(NULL);
    std::string name_textfile = "r-index-f-" + random_alphanumeric_string(10);
	std::ofstream text_file(name_textfile);
	write_to_file(text_file,input.c_str(),input.size());
	text_file.close();
	std::string cmd_newscan = "build/external/Big-BWT/newscanNT.x " + name_textfile + " >nul 2>nul";
	system(cmd_newscan.c_str());
	std::string cmd_pfp_thresholds = "build/external/pfp-thresholds/pfp-thresholds " +
	name_textfile + " -r > " + name_textfile + ".log 2>" + name_textfile + ".log";
	system(cmd_pfp_thresholds.c_str());
	index = std::move(r_index_f<>(name_textfile));
	std::ifstream log_file(name_textfile + ".log");
	external_peak_memory_usage = peak_memory_usage(log_file);
	log_file.close();
	std::string cmd_rm = "rm -f " + name_textfile + "* >nul 2>nul";
	system(cmd_rm.c_str());
	system("rm -f nul");
    std::cout.rdbuf(cout_rfbuf);
}

template <>
void destroy_index<r_index_f<>>(r_index_f<>& index) {
	index = std::move(r_index_f<>());
}

template <>
uint32_t count_pattern<uint32_t,r_index_f<>>(r_index_f<>& index, std::string& pattern) {
	return index.count(pattern);
}

template <>
uint64_t count_pattern<uint64_t,r_index_f<>>(r_index_f<>& index, std::string& pattern) {
	return index.count(pattern);
}

template <>
void revert_index<r_index_f<>>(r_index_f<>& index, uint16_t) {
    index.invert(input_reverted);
}

// ############################# benchmarking framework #############################

template <typename uint_t>
void measure_all() {
	measure<uint_t,move_r<uint_t>,false,true,true,true>("move-r (a=8)","move_r",max_num_threads,max_num_threads);
	measure<uint_t,move_r<uint_t>,true,false,false,false>("move-r (a=8, Big-BWT)","move_r_bigbwt",max_num_threads,0);
	measure<uint_t,r_index_f<>,false,true,true,false>("r-index-f","r_index_f",1,1);
	measure<uint_t,rcomp_lfig,false,true,true,true>("rcomp (LFIG)","rcomp_lfig",1,1);
	measure<uint_t,rcomp_glfig_16,false,true,true,true>("rcomp (GLFIG, g=16)","rcomp_glfig_g16",1,1);
	measure<uint_t,ri::r_index<>,false,true,true,true>("r-index (prezza)","r_index_prezza",1,1);
	measure<uint_t,ri_mun::r_index<>,false,false,false,false>("r-index (mun)","r_index_mun",1,1);
	measure<uint_t,OnlineRlbwt,false,true,true,true>("OnlineRlbwt","online_rlbwt",1,1);
	measure<uint_t,DYNAMIC,false,true,true,false>("DYNAMIC","dynamic",1,1);
}

int main(int argc, char** argv) {
	if (argc < 5) help("");
	while (ptr < argc - 4) parse_args(argv, argc, ptr);

	path_inputfile = argv[ptr];
	path_patternsfile_1 = argv[ptr+1];
	path_patternsfile_2 = argv[ptr+2];
	max_num_threads = atoi(argv[ptr+3]);

	input_file.open(path_inputfile);
	patterns_file_1.open(path_patternsfile_1);
	patterns_file_2.open(path_patternsfile_2);

	if (!input_file.good()) help("error: invalid input, could not read <input_file>");
	if (!patterns_file_1.good()) help("error: invalid input, could read <patterns_file_1>");
	if (!patterns_file_2.good()) help("error: invalid input, could read <patterns_file_2>");
	if (measurement_file.is_open() && check_correctness) help("error: cannot output measurement data when checking for correctness");
	if (max_num_threads == 0 || max_num_threads > omp_get_max_threads()) help("error: invalid number of threads");

	system("chmod +x build/external/Big-BWT/*");
	system("chmod +x build/external/pfp-thresholds/*");
	std::cout << std::setprecision(4);
	name_textfile = path_inputfile.substr(path_inputfile.find_last_of("/\\") + 1);
	std::cout << "benchmarking " << path_inputfile;

	input_file.seekg(0,std::ios::end);
	uint64_t input_size = input_file.tellg()+(std::streamsize)+1;
	input_file.seekg(0,std::ios::beg);
	input.reserve(input_size);
	no_init_resize(input,input_size-1);
	read_from_file(input_file,input.c_str(),input_size-1);
	input_file.close();
	
	std::cout << " (" << format_size(input_size-1) << ") using up to " << format_threads(max_num_threads) << std::endl;
	if (check_correctness) std::cout << "correctnes will be checked if possible" << std::endl;
	std::cout << std::endl;

	if (input_size < UINT_MAX) {
		measure_all<uint32_t>();
	} else {
		measure_all<uint64_t>();
	}

	patterns_file_1.close();
	patterns_file_2.close();

	if (measurement_file.is_open()) measurement_file.close();
}