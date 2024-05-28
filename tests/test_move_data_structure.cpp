#include <gtest/gtest.h>
#include <move_r/data_structures/move_data_structure/move_data_structure.hpp>
#include <move_r/data_structures/move_data_structure/move_data_structure_l_.hpp>

TEST(test_move_data_structure,fuzzy_test) {
    std::random_device rd;
    std::mt19937 gen(rd());
    uint16_t max_num_threads = omp_get_max_threads();

    std::uniform_int_distribution<uint32_t> input_size_distrib(1,200000);
    std::lognormal_distribution<double> avg_interval_length_distrib(4.0,2.0);
    std::uniform_int_distribution<uint16_t> num_threads_distrib(1,max_num_threads);
    std::lognormal_distribution<double> a_distrib(2.0,3.0);

    uint32_t input_size;
    uint32_t num_intervals;
    std::vector<std::pair<uint32_t,uint32_t>> interval_sequence;
    std::vector<uint32_t> interval_permutation;

    auto start_time = now();

    // generate random disjoint interval sequences and test the move data structure until one hour has passed
    while (time_diff_min(start_time,now()) < 60) {
        // choose a random input size
        input_size = input_size_distrib(gen);
        num_intervals = 0;

        // choose input intervals of random lengths
        std::uniform_int_distribution<uint32_t> interval_length_distrib(1,2*avg_interval_length_distrib(gen));
        uint32_t interval_length_prefix_sum = 0;
        interval_sequence.emplace_back(std::make_pair(0,0));
        while (interval_length_prefix_sum < input_size) {
            interval_length_prefix_sum = std::min<uint32_t>(
                interval_length_prefix_sum+interval_length_distrib(gen),
                input_size
            );
            interval_sequence.emplace_back(std::make_pair(interval_length_prefix_sum,0));
            num_intervals++;
        }

        // permute the input intervals randomly into the output intervals
        interval_permutation.resize(num_intervals);
        #pragma omp parallel for num_threads(max_num_threads)
        for (uint32_t i=0; i<num_intervals; i++) {
            interval_permutation[i] = i;
        }
        std::shuffle(interval_permutation.begin(),interval_permutation.end(),gen);
        interval_length_prefix_sum = 0;
        for (uint32_t i=0; i<num_intervals; i++) {
            interval_sequence[interval_permutation[i]].second = interval_length_prefix_sum;
            interval_length_prefix_sum += interval_sequence[interval_permutation[i]+1].first-interval_sequence[interval_permutation[i]].first;
        }
        interval_sequence.pop_back();

        // choose a random value for the balancing parameter a
        uint16_t a = std::min<uint16_t>(2+a_distrib(gen),32767);

        // build a move data structure from the disjoint interval sequence
        move_data_structure<uint32_t> mds(interval_sequence,input_size,{
            .num_threads = num_threads_distrib(gen), .a = a
        });

        // check if the number of input/output intervals has increased too much
        EXPECT_TRUE(mds.num_intervals()/(double)num_intervals <= (a/(double)(a-1))*1.125);

        // check if there is an a-heavy output interval
        #pragma omp parallel for num_threads(max_num_threads)
        for (uint32_t i=0; i<mds.num_intervals(); i++) {
            uint32_t j = bin_search_min_geq<uint32_t>(mds.q(i),0,mds.num_intervals()-1,[&mds](auto x){return mds.p(x);});
            uint16_t num_intervals_in_output_interval = 0;
            while (mds.p(j) < mds.q(i)+(mds.p(i+1)-mds.p(i))) {
                num_intervals_in_output_interval++;
                j++;
            }
            EXPECT_LE(num_intervals_in_output_interval,2*a);
        }

        // perform at most ~10000 move queries starting with roughly evenly distributed input values in the range [0,input_size)
        // using the move data structure and the original interval sequence and compare the output values
        uint32_t avg_step_size = std::max<uint32_t>(2,input_size/10000);
        std::uniform_int_distribution<uint32_t> step_size_distrib(avg_step_size/1.5,1.5*avg_step_size);
        #pragma omp parallel for num_threads(max_num_threads)
        for (uint32_t i=0; i<input_size; i+=step_size_distrib(gen)) {
            std::pair<uint32_t,uint32_t> ix_mds{i,
                bin_search_max_leq<uint32_t>(i,0,mds.num_intervals()-1,[&mds](uint32_t x){return mds.p(x);})
            };
            std::pair<uint32_t,uint32_t> ix_is{i,
                bin_search_max_leq<uint32_t>(i,0,num_intervals-1,[&interval_sequence](uint32_t x){return interval_sequence[x].first;})
            };
            ix_mds = mds.move(ix_mds);
            ix_is.first = interval_sequence[ix_is.second].second + (ix_is.first - interval_sequence[ix_is.second].first);
            EXPECT_EQ(ix_mds.first,ix_is.first);
            // also check if the index of the input interval containing the output value has been calculated correctly
            EXPECT_TRUE(mds.p(ix_mds.second) <= ix_mds.first && ix_mds.first < mds.p(ix_mds.second+1));
        }

        interval_sequence.clear();
    }
}