#include <move_r/misc/utils.hpp>
#include <move_r/data_structures/move_data_structure/move_data_structure_phi.hpp>
#include <move_r/data_structures/move_data_structure/move_data_structure_lf.hpp>

int main() {
    // Build a move data structure from the disjoint interval
    // sequence I = (0,1),(1,0); => n = 2
    move_data_structure_phi<> mds_phi({{0,1},{1,0}},2);

    // create a pair to perform move queries with
    std::pair<uint32_t,uint32_t> ix{0,0};

    // perform some move queries
    std::cout << to_string<>(ix = mds_phi.move(ix)) << std::endl;
    std::cout << to_string<>(ix = mds_phi.move(ix)) << std::endl;

    // build a move_data_structure_lf (intended for I_LF);
    // this move data structure additionally stores a string interleaved
    // with the arrays needed for performing move queries (intended for
    // storing the characters of the bwt (sub-)runs);

    // use at most 4 threads and set a := 2
    move_data_structure_lf<> mds_lf({{0,4},{1,5},{2,6},{3,7},{4,0}},8,4,2);

    // this disjoint interval sequence is not 2-balanced, because the output
    // interval [0,3] containes 4 >= 2a = 4 input intervals

    // print the pairs of the resulting disjoint interval sequence
    for (uint32_t i=0; i<mds_lf.num_intervals(); i++) {
        std::cout << to_string<>(std::make_pair(mds_lf.p(i),mds_lf.q(i)));
    }
    
    // the balancing algorithm has added the pair (6,2) to balance the sequence
}