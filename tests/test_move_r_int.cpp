#include <gtest/gtest.h>
#include <move_r/move_r.hpp>

std::random_device rd;
std::mt19937 gen(rd());
uint16_t max_num_threads = omp_get_max_threads();

std::lognormal_distribution<double> avg_input_rep_length_distrib(4.0, 2.0);
std::uniform_real_distribution<double> prob_distrib(0.0, 1.0);
std::uniform_int_distribution<uint32_t> input_size_distrib(1, 200000);
std::uniform_int_distribution<uint16_t> num_threads_distrib(1, max_num_threads);
std::lognormal_distribution<double> a_distrib(2.0, 3.0);

uint32_t input_size;
uint32_t alphabet_size;
std::vector<int32_t> input;
std::vector<int32_t> input_libsais;
std::vector<int32_t> input_reverted;
std::vector<int32_t> bwt;
std::vector<int32_t> bwt_retrieved;
std::vector<int32_t> suffix_array;
std::vector<uint32_t> suffix_array_retrieved;
uint32_t max_pattern_length;
uint32_t num_queries;

template <move_r_support support>
void test_move_r_int()
{
    // choose a random input length
    input_size = input_size_distrib(gen);

    // choose a random alphabet
    alphabet_size = std::uniform_int_distribution<int32_t>(1, input_size)(gen);
    std::uniform_int_distribution<int32_t> symbol_distrib(INT_MIN, INT_MAX);
    std::vector<int32_t> alphabet;
    alphabet.reserve(alphabet_size);

    for (uint32_t i = 0; i < alphabet_size; i++)
        alphabet.emplace_back(symbol_distrib(gen));

    ips4o::parallel::sort(alphabet.begin(), alphabet.end());
    alphabet.erase(std::unique(alphabet.begin(), alphabet.end()), alphabet.end());
    alphabet_size = alphabet.size();

    // choose a random input based on the alphabet
    std::uniform_int_distribution<uint32_t> symbol_idx_distrib(0, alphabet_size - 1);
    double avg_input_rep_length = 1.0 + avg_input_rep_length_distrib(gen);
    int32_t cur_symbol = alphabet[symbol_idx_distrib(gen)];
    input.reserve(input_size);

    for (uint32_t i = 0; i < input_size; i++) {
        if (prob_distrib(gen) < 1 / avg_input_rep_length)
            cur_symbol = alphabet[symbol_idx_distrib(gen)];

        input.emplace_back(cur_symbol);
    }

    // build move-r and choose a random number of threads and balancing parameter, but always use libsais,
    // because there are bugs in Big-BWT that come through during fuzzing but not really in practice
    move_r<support, int32_t, uint32_t> index(input, {
        .mode = _suffix_array,
        .num_threads = num_threads_distrib(gen),
        .a = std::min<uint16_t>(2 + a_distrib(gen), 32767) 
    });

    // revert the index and compare the output with the input string
    input_reverted = index.revert({ .num_threads = num_threads_distrib(gen) });

    #pragma omp parallel for num_threads(max_num_threads)
    for (uint32_t i = 0; i < input_size; i++)
        EXPECT_EQ(input[i], input_reverted[i]);

    // to build the suffix array, remap the symbols in the input to [1,2,...,alphabet_size]
    gtl::flat_hash_map<int32_t, int32_t> map_int;
    int32_t sym_cur = 1;

    for (uint32_t i = 0; i < alphabet_size; i++)
        map_int[alphabet[i]] = sym_cur++;

    no_init_resize(input_libsais, input_size + 1);

    #pragma omp parallel for num_threads(max_num_threads)
    for (uint32_t i = 0; i < input_size; i++)
        input_libsais[i] = map_int[input[i]];

    input_libsais[input_size] = 0;

    // retrieve the suffix array and compare it with the correct suffix array
    no_init_resize(suffix_array, input_size + 1 + 6 * (alphabet_size + 1));
    suffix_array[input_size] = 0;
    libsais_int_omp(&input_libsais[0], &suffix_array[0], input_size + 1, alphabet_size + 1, 6 * (alphabet_size + 1), max_num_threads);
    no_init_resize(suffix_array, input_size + 1);
    suffix_array_retrieved = index.SA({ .num_threads = num_threads_distrib(gen) });

    #pragma omp parallel for num_threads(max_num_threads)
    for (uint32_t i = 0; i <= input_size; i++)
        EXPECT_EQ(suffix_array[i], suffix_array_retrieved[i]);

    // compute each suffix array value separately and check if it is correct
    #pragma omp parallel for num_threads(max_num_threads)
    for (uint32_t i = 0; i <= input_size; i++)
        EXPECT_EQ(index.SA(i), suffix_array[i]);

    // retrieve the bwt and compare it with the correct bwt
    no_init_resize(bwt, input_size + 1);

    #pragma omp parallel for num_threads(max_num_threads)
    for (uint32_t i = 0; i <= input_size; i++) {
        bwt[i] = suffix_array[i] == 0 ? 0 : input[suffix_array[i] - 1];
    }

    bwt_retrieved = index.BWT({ .num_threads = num_threads_distrib(gen) });

    #pragma omp parallel for num_threads(max_num_threads)
    for (uint32_t i = 0; i <= input_size; i++)
        EXPECT_EQ(bwt[i], bwt_retrieved[i]);

    // compute each bwt character separately and check if it is correct
    #pragma omp parallel for num_threads(max_num_threads)
    for (uint32_t i = 0; i <= input_size; i++)
        EXPECT_EQ(index.BWT(i), bwt[i]);

    // generate patterns from the input and test count- and locate queries
    std::uniform_int_distribution<uint32_t> pattern_pos_distrib(0, input_size - 1);
    max_pattern_length = std::min<uint32_t>(10000, std::max<uint32_t>(100, input_size / 1000));
    std::uniform_int_distribution<uint32_t> pattern_length_distrib(1, max_pattern_length);
    num_queries = std::min<uint32_t>(10000, std::max<uint32_t>(1000, input_size / 100));
    num_queries = std::max<uint32_t>(1, num_queries / max_num_threads);

    #pragma omp parallel num_threads(max_num_threads)
    {
        std::random_device rd_thr;
        std::mt19937 gen_thr(rd());
        uint32_t pattern_pos;
        uint32_t pattern_length;
        std::vector<int32_t> pattern;
        std::vector<uint32_t> correct_occurrences;
        std::vector<uint32_t> occurrences;
        bool match;

        for (uint32_t cur_query = 0; cur_query < num_queries; cur_query++) {
            pattern_pos = pattern_pos_distrib(gen_thr);
            pattern_length = std::min<uint32_t>(input_size - pattern_pos, pattern_length_distrib(gen));
            no_init_resize(pattern, pattern_length);

            for (uint32_t i = 0; i < pattern_length; i++)
                pattern[i] = input[pattern_pos + i];

            for (uint32_t i = 0; i <= input_size - pattern_length; i++) {
                match = true;

                for (uint32_t j = 0; j < pattern_length; j++) {
                    if (input[i + j] != pattern[j]) {
                        match = false;
                        break;
                    }
                }

                if (match)
                    correct_occurrences.emplace_back(i);
            }

            EXPECT_EQ(index.count(pattern), correct_occurrences.size());

            occurrences = index.locate(pattern);
            ips4o::sort(occurrences.begin(), occurrences.end());
            EXPECT_EQ(occurrences, correct_occurrences);
            occurrences.clear();

            auto query = index.query();
            for (int32_t i = pattern.size() - 1; i >= 0; i--)
                query.prepend(pattern[i]);
            EXPECT_EQ(query.num_occ(), correct_occurrences.size());

            query.locate(occurrences);
            ips4o::sort(occurrences.begin(), occurrences.end());
            EXPECT_EQ(occurrences, correct_occurrences);
            correct_occurrences.clear();
        }
    }

    input.clear();
}

TEST(test_move_r, fuzzy_test)
{
    auto start_time = now();

    while (time_diff_min(start_time, now()) < 60) {
        if (prob_distrib(gen) < 0.5) {
            test_move_r_int<_locate_move>();
        } else {
            test_move_r_int<_locate_rlzdsa>();
        }
    }
}