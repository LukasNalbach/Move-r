#pragma once

#include <move_r/data_structures/move_data_structure/move_data_structure.hpp>

template <typename pos_t>
void move_data_structure<pos_t>::construction::verify_correctness() {
    std::cout << "verifying correctness of the interval sequence:" << std::endl;
    bool correct = true;

    // check if the input interval starting positions ascend
    #pragma omp parallel for num_threads(p)
    for (uint64_t i=0; i<k_; i++) {
        if (!(mds.p(i) < mds.p(i+1))) {
            #pragma omp critical
            {
                std::cout << "input interval starting positions do not ascend:" << std::endl;
                std::cout << "i = " << i << std::endl;
                std::cout << "p_i = " << mds.p(i) << std::endl;
                std::cout << "p_{i+1} = " << mds.p(i+1) << std::endl << std::endl;
                correct = false;
            }
        }
    }

    // check if an input interval is too long
    #pragma omp parallel for num_threads(p)
    for (uint64_t i=0; i<k_; i++) {
        if (mds.p(i+1) - mds.p(i) > l_max) {
            #pragma omp critical
            {
                std::cout << "too long input interval (> l_max = " << l_max << "):" << std::endl;
                std::cout << "i = " << i << std::endl;
                std::cout << "p_i = " << mds.p(i) << std::endl;
                std::cout << "p_{i+1} = " << mds.p(i+1) << std::endl << std::endl;
                correct = false;
            }
        }
    }
    
    // check if the output interval lengths do not match the input interval lengths
    #pragma omp parallel for num_threads(p)
    for (uint64_t i=0; i<k_; i++) {
        if (D_q[pi[i+1]] - D_q[pi[i]] != mds.p(pi[i]+1) - mds.p(pi[i])) {
            #pragma omp critical
            {
                std::cout << "input- and output interval lengths do not match:" << std::endl;
                std::cout << "i = " << i << std::endl;
                std::cout << "p_{pi[i]} = " << mds.p(pi[i]) << std::endl;
                std::cout << "p_{pi[i]+1} = " << mds.p(pi[i]+1) << std::endl;
                std::cout << "q_{pi[i]} = " << D_q[pi[i]] << std::endl;
                std::cout << "q_{pi[i+1]} = " << D_q[pi[i+1]] << std::endl << std::endl;
                correct = false;
            }
        }
    }

    pos_t i = 0;
    pos_t j = 0;
    pos_t e;

    // check if there is an a-heavy output interval by iterating over the input- and output intervals with pi
    while (i < k && j < k) {
        while (i < k && mds.p(i) < D_q[pi[j]]) {
            i++;
        }

        if (mds.p(i) < D_q[pi[j+1]]) {
            e = 1;
            
            while (i+1 < k && mds.p(i+1) < D_q[pi[j+1]]) {
                i++;
                e++;

                if (e >= two_a) {
                    std::cout << "a-heavy output interval:" << std::endl;
                    std::cout << "i = " << i-two_a+1 << std::endl;
                    std::cout << "j = " << j << std::endl;
                    std::cout << "q_{pi[j]} = " << D_q[pi[j]] << std::endl;
                    std::cout << "p_{i} = " << mds.p(i-two_a+1) << std::endl;
                    std::cout << "p_{i+2a} = " << mds.p(i) << std::endl;
                    std::cout << "q_{pi[j+1]} = " << D_q[pi[j+1]] << std::endl << std::endl;
                    correct = false;
                    j++;
                    break;
                }
            }

            i++;
        } else {
            j++;
        }
    }
    
    // check whether D_idx has been calculated correctly
    #pragma omp parallel for num_threads(p)
    for (uint64_t i=0; i<k_; i++) {
        if (!(mds.p(mds.idx(i)) <= D_q[i] && D_q[i] < mds.p(mds.idx(i)+1))) {
            #pragma omp critical
            {
                std::cout << "wrong value:" << std::endl;
                std::cout << "i = " << i << std::endl;
                std::cout << "D_p[D_idx[i]] = " << mds.p(mds.idx(i)) << std::endl;
                std::cout << "D_q(i) = " << D_q[i] << std::endl;
                std::cout << "D_p[D_idx[i]+1] = " << mds.p(mds.idx(i)+1) << std::endl << std::endl;
                correct = false;
            }
        }
    }

    // check if D_offs has been calculated correctly (by checking if we can recalculate D_q together with D_idx and D_p)
    #pragma omp parallel for num_threads(p)
    for (uint64_t i=0; i<k_; i++) {
        if (mds.q(i) != D_q[i]) {
            #pragma omp critical
            {
                std::cout << "wrong value:" << std::endl;
                std::cout << "i = " << i << std::endl;
                std::cout << "D_idx[i] = " << mds.idx(i) << std::endl;
                std::cout << "D_p[D_idx[i]] = " << mds.p(mds.idx(i)) << std::endl;
                std::cout << "D_offs[i] = " << mds.offs(i) << std::endl;
                std::cout << "mds.q(i) = D_p[D_idx[i]] + D_offs[i] = " << mds.q(i);
                std::cout << " != D_q[i] = " << D_q[i] << std::endl << std::endl;
                correct = false;
            }
        }
    }

    std::cout << "the interval sequence has " << (correct ? "" : "not ") << "been built correctly" << std::endl;
}