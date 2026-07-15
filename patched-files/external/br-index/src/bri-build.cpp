#include <chrono>
#include <iostream>
#include <string>

#include <malloc_count.h> // for malloc_count_peak / malloc_count_current

#include "br_index.hpp"
#include "br_index_nplcp.hpp"
#include "utils.hpp"

using namespace std;
using namespace bri;

string out_basename = string();
string input_file = string();
bool sais = true;
bool nplcp = false;

void help(){
	cout << "bri-build: builds the bidirectional r-index. Extension .bri/.brin is automatically added to output index file" << endl << endl;
	cout << "Usage: bri-build [options] <input_file_name>" << endl;
	cout << "   -o <basename>        use 'basename' as prefix for all index files. Default: basename is the specified input_file_name"<<endl;
	cout << "   -divsufsort          use divsufsort algorithm to build the BWT (fast, 7.5n Bytes of RAM). By default,"<<endl;
	cout << "                        SE-SAIS is used (about 4 time slower than divsufsort, 4n Bytes of RAM)."<<endl;
    cout << "   -nplcp               use the version without PLCP. When locating, calculate LF^d(p) first."<<endl;
    cout << "                        fast when occ is very high, but takes slightly larger space than the normal version."<<endl;
	cout << "   <input_file_name>    input text file." << endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	if (s.compare("-o") == 0) 
    {

		if(ptr >= argc-1){
			cout << "Error: missing parameter after -o option." << endl;
			help();
		}

		out_basename = string(argv[ptr]);
		ptr++;

	}
    else if (s.compare("-divsufsort") == 0 )
    {

		sais = false;

	}
    else if (s.compare("-nplcp") == 0)
    {

        nplcp = true;

    }
    else
    {
		cout << "Error: unrecognized '" << s << "' option." << endl;
		help();
	}

}

int main(int argc, char** argv) 
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;

    // heap usage of this process before construction; peak growth from here on is
    // reported via malloc_count (0 B when malloc_count is disabled at compile time),
    // mirroring how move-r and columba measure their construction peak.
    const size_t mem_baseline = malloc_count_current();

    auto t1 = high_resolution_clock::now();

    out_basename = string();
    input_file = string();

    int ptr = 1;

    if (argc < 2) help();

    while (ptr < argc-1) parse_args(argv, argc, ptr);

    input_file = string(argv[ptr]);

    if (out_basename.compare("") == 0)
        out_basename = string(input_file);
    
    string idx_file = out_basename;

    if (nplcp) idx_file.append(".brin");
    else idx_file.append(".bri");

    cout << "Building br-index of input file " << input_file << endl;
    cout << "Index will be saved to " << idx_file << endl;

    string input;

    {
        std::ifstream fs(input_file);
        std::stringstream buffer;
        buffer << fs.rdbuf();

        input = buffer.str();
    }

    std::ofstream out(idx_file);

    if (nplcp)
    {
        br_index_nplcp<> idx(input,sais);
        idx.serialize(out);
    } 
    else 
    {
        br_index<> idx(input,sais);
        idx.serialize(out);
    }

    auto t2 = high_resolution_clock::now();

    ulint total = duration_cast<duration<double, std::ratio<1>>>(t2-t1).count();
    cout << "Build time: " << get_time(total) << endl;

    const size_t raw_peak = malloc_count_peak();
    const uint64_t peak = raw_peak > mem_baseline ? raw_peak - mem_baseline : raw_peak;
    cout << "Peak memory usage during construction: " << peak << " bytes" << endl;

    out.close();
}