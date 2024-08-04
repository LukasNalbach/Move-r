#include <move_r/move_r.hpp>

int main() {
    // build an index
    move_r<> index("This is a test string");

    // build a 64-bit index (intended for large input strings > UINT_MAX
    // bytes ~ 4GB) with only count support, use the Big-BWT
    // construction algorithm, use at most 8 threads and set the 
    // balancing parameter a to 4
    move_r<_count,char,uint64_t> index_2("a large string",{
        .mode = _bigbwt, .num_threads = 8, .a = 4
    });

    // print the number of bwt runs in the input string
    std::cout << index.num_bwt_runs() << std::endl;

    // print the index size
    std::cout << format_size(index.size_in_bytes()) << std::endl;

    // print the number of occurences of a pattern
    std::cout << index.count("test") << std::endl;

    // store all occurences of a pattern in a vector
    auto Occ = index.locate("is");
    for (auto o : Occ) std::cout << o << ", ";
    std::cout << std::endl;

    // build an index for an integer vector using a relative
    // lempel-ziv encoded differential suffix array (rlzdsa)
    move_r<_locate_rlzdsa,int32_t> index_3({2,-1,5,-1,7,2,-1});

    // incrementally search the pattern [2,-1] in the input vector (from
    // right to left) and print the number of occurrences after each step
    auto query = index_3.query();
    query.prepend(-1);
    std::cout << query.num_occ() << std::endl;
    query.prepend(2);
    std::cout << query.num_occ() << std::endl;

    // print the suffix-array interval [b,e] of [2,-1]
    std::cout << "b = " << query.sa_interval().first
            << ", e = " << query.sa_interval().second << std::endl;

    // incrementally locate the occurrences of [2,-1] in the input vector
    while (query.num_occ_rem() > 0) {
        std::cout << query.next_occ() << ", " << std::flush;
    }

    // compute the longest suffix of [0,7,2] that occurs in the input vector
    std::vector<int32_t> pattern = {0,7,2};
    auto query_2 = index_3.query();
    uint32_t suffix = pattern.size();
    while (suffix > 0 && query_2.prepend(pattern[suffix-1])) suffix--;
    std::cout << std::endl << suffix << std::flush;
}