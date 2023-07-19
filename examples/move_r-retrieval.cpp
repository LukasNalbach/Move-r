#include <move_r/move_r.hpp>

int main() {
    // build an index
    move_r<> index("This is a test string");

    // retrieve the range [8,17] of the original text and store 
    // it in a string using at most 2 threads
    std::string reverted_range = index.revert_range(8,17,2);
    
    std::cout << std::endl;
    for (auto c : reverted_range) {
        std::cout << c;
    }

    // print the original text from right to left without storing it
    // using 1 thread
    std::cout << std::endl;
    index.revert_range([](auto,auto c){std::cout << c;},0,20,1);

    // retrieve the suffix array values in the range [2,6] using at
    // most 4 threads and store them in a vector
    std::vector<uint32_t> sa_range = index.retrieve_sa_range(2,6,4);
    std::cout << std::endl;
    
    for (auto s : sa_range) {
        std::cout << s << ", ";
    }

    // print the suffix array values in the range [7,14] from right to left
    // using 1 thread
    std::cout << std::endl;
    index.retrieve_sa_range([](auto,auto s){std::cout << s << ", ";},7,14,1);
}