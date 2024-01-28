#include <move_r/move_r.hpp>

int main() {
    // build an index
    move_r<> index("This is a test string");

    // retrieve the range [8,17] of the original text and store 
    // it in a string using at most 2 threads
    std::string reverted_range = index.revert({
        .l = 8, .r = 17, .num_threads = 2
    });
    for (auto c : reverted_range) std::cout << c;
    std::cout << std::endl;

    // print the original text from right to left without storing it
    // using 1 thread
    index.revert([](auto,auto c){std::cout << c;},{.num_threads = 1});
    std::cout << std::endl;

    // retrieve the suffix array values in the range [2,6] using at
    // most 4 threads and store them in a vector
    std::vector<uint32_t> SA_range = index.SA({
        .l = 2, .r = 6, .num_threads = 4
    });
    for (auto s : SA_range) std::cout << s << ", ";
    std::cout << std::endl;

    // print SA[1]
    std::cout << index.SA(1) << std::endl;

    // retrieve the BWT in the range [7,14] from left to right
    // using 1 thread
    index.BWT([](auto,auto s){std::cout << s << ", ";},{
        .l = 7, .r = 14, .num_threads = 1
    });
    std::cout << std::endl;

    // print BWT[16]
    std::cout << index.BWT(16) << std::endl;
}