#include <move_r/move_r.hpp>

int main() {
    // build an index
    move_r<> index("This is a test string");

    // build a 64-bit index (intended for large input strings > UINT_MAX
    // bytes ~ 4GB) with only count support, use Big-BWT
    // construction algorithm, use at most 8 threads and set the 
    // balancing parameter a to 4
    move_r<uint64_t> index_2("a large string",{
        .support = {count},
        .mode = _bigbwt,
        .num_threads = 8,
        .a = 4
    });

    // print the number of bwt runs in the input string
    std::cout << index.num_bwt_runs() << std::endl;

    // print the index size
    std::cout << format_size(index.size_in_bytes()) << std::endl;

    // print the number of occurences of a pattern
    std::cout << index.count("test") << std::endl;

    // print all occurences of a pattern
    index.locate("is",[](auto o){std::cout << o << ", ";});

    // store all occurences of a pattern in a vector
    auto Occ = index.locate("test");

    std::cout << std::endl;
    for (auto o : Occ) {
        std::cout << o << ", ";
    }
}