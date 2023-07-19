#include <iostream>
#include <filesystem>
#include <move_r/misc/utils.hpp>
#include <move_r/move_r.hpp>

int ptr = 1;
bool output_occurences = false;
bool check_correctness = false;
std::string input;
std::ofstream measurement_file;
std::string path_index_file;
std::string path_patterns_file;
std::string path_textfile;
std::string path_outputfile;
std::ifstream index_file;
std::ifstream patterns_file;
std::ifstream input_file;
std::ofstream output_file;
std::string name_textfile;

void help(std::string msg) {
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

void parse_args(char **argv, int argc, int &ptr) {
	std::string s = argv[ptr];
	ptr++;

	if (s == "-c") {
		if (ptr >= argc-1) help("error: missing parameter after -c option.");
		check_correctness = true;
		path_textfile = argv[ptr++];
	} else if (s == "-m") {
		if (ptr >= argc - 1) help("error: missing parameter after -o option.");
		std::string path_m_file = argv[ptr++];
		measurement_file.open(path_m_file,std::filesystem::exists(path_m_file) ? std::ios::app : std::ios::out);
		if (!measurement_file.good()) help("error: cannot open measurement file");
		name_textfile = argv[ptr++];
	} else if (s == "-o") {
		if (ptr >= argc-1) help("error: missing parameter after -o option.");
		output_occurences = true;
		path_outputfile = argv[ptr++];
	} else  {
		help("error: unrecognized '" + s + "' option");
	}
}

template <typename uint_t>
void measure_locate() {
	std::cout << std::setprecision(4);
	std::cout << "loading the index" << std::flush;
	auto t1 = now();
	move_r<uint_t> index;
	index.load(index_file,{locate});
	log_runtime(t1);
	index_file.close();
	std::cout << std::endl;
	index.log_data_structure_sizes();

	if (check_correctness) {
		std::cout << std::endl << "reading the original input file" << std::flush;
		input_file.seekg(0,std::ios::end);
		uint64_t n = input_file.tellg()+(std::streamoff)1;
		input_file.seekg(0,std::ios::beg);
		no_init_resize(input,n-1);
		read_from_file(input_file,input.c_str(),n-1);
		input_file.close();
    }

	std::cout << std::endl << "searching patterns ... " << std::endl;
	std::string header;
	std::getline(patterns_file,header);
	uint64_t num_patterns = number_of_patterns(header);
	uint64_t pattern_length = patterns_length(header);
	uint64_t perc;
	uint64_t last_perc = 0;
	uint64_t num_occurrences = 0;
	uint64_t time_locate = 0;
	std::chrono::steady_clock::time_point t2,t3;
	std::string pattern;
	no_init_resize(pattern,pattern_length);
	std::vector<uint_t> occurrences;
	bool is_sorted, equal;
	uint_t count;

	for (uint64_t i=0; i<num_patterns; i++) {
		perc = (100*i) / num_patterns;
		is_sorted = false;

		if (perc > last_perc) {
			std::cout << perc << "% done .." << std::endl;
			last_perc = perc;
		}

		patterns_file.read((char*)&pattern[0],pattern_length);
		t2 = now();
		index.locate(pattern,[&occurrences](uint_t occurrence){occurrences.emplace_back(occurrence);});
		t3 = now();
		time_locate += time_diff_ns(t2,t3);
		num_occurrences += occurrences.size();

		if (check_correctness) {
			ips4o::sort(occurrences.begin(),occurrences.end());
			is_sorted = true;

			if (occurrences.size() != (count = index.count(pattern))) std::cout << "error: wrong number of located occurrences: " << occurrences.size() << "/" << count << std::endl;

			for (uint_t occurrence : occurrences) {
				equal = true;

				for (uint_t pos=0; pos<pattern_length; pos++) {
					if (input[occurrence+pos] != pattern[pos]) {
						equal = false;
						break;
					}
				}

				if (!equal) {
					std::cout << "error: wrong occurrence: " << occurrence << " ("  << num_occurrences << " occurrences) "<< std::endl;
					for (uint_t pos=0; pos<pattern_length; pos++) std::cout << input[occurrence+pos];
					std::cout << std::endl << std::endl << "/" << std::endl << std::endl;
					for (uint_t pos=0; pos<pattern_length; pos++) std::cout << pattern[pos];
					std::cout << std::endl;
					break;
				}
			}
		}

		if (output_occurences) {
			if (!is_sorted) ips4o::sort(occurrences.begin(),occurrences.end());
			output_file.write((char*)&occurrences[0],occurrences.size());
		}

		occurrences.clear();
	}

	if (num_occurrences == 0) {
		std::cout << "found no occurrences" << std::endl;
	} else {
		std::cout << "average occurrences per pattern: " << (num_occurrences/num_patterns) << std::endl;
		std::cout << "number of patterns: " << num_patterns << std::endl;
		std::cout << "pattern length: " << pattern_length << std::endl;
		std::cout << "total number of occurrences: " << num_occurrences << std::endl;
		std::cout << "locate time: " << format_time(time_locate) << std::endl;
		std::cout << "            " << format_time(time_locate/num_patterns) << "/pattern" << std::endl;
		std::cout << "            " << format_time(time_locate/num_occurrences) << "/occurrence" << std::endl;
	}

	if (measurement_file.is_open()) {
		measurement_file << "RESULT";
		measurement_file << " type=locate";
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
		measurement_file << " time_locate=" << time_locate;
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

	if (output_occurences) {
		output_file.open(path_outputfile);
		if (!output_file.good()) help("error: could not create <output_file>");
	}

	if (check_correctness) {
		input_file.open(path_textfile);
		if (!input_file.good()) help("error: could not read <input_file>");
	}

	bool is_64_bit;
	index_file.read((char*)&is_64_bit,1);
	index_file.seekg(0,std::ios::beg);

	if (is_64_bit) {
		measure_locate<uint64_t>();
	} else {
		measure_locate<uint32_t>();
	}

	patterns_file.close();
	if (output_occurences) output_file.close();
}