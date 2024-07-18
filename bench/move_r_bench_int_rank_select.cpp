#include <ctime>
#include <experimental/random>

#include <move_r/data_structures/rank_select_support.hpp>
#include <move_r/misc/utils.hpp>

static constexpr uint32_t num_queries = 1 << 24;
static constexpr uint32_t target_input_size = 1 << 24;
std::random_device rd;
std::mt19937 gen(rd());
std::vector<uint32_t> input;
std::vector<uint32_t> num_occ;
struct rank_query {uint32_t sym;uint32_t pos;};
struct select_query {uint32_t sym;uint32_t rank;};
std::vector<rank_query> rank_queries;
std::vector<select_query> select_queries;
enum rank_mode {gmr,naive,hybrid};

template <rank_mode rm>
void bench_rank(
    rank_select_support<uint32_t>& rank_select,
    uint32_t avg_occ, uint32_t alphabet_size
) {
    std::cout << "benchmarking rank queries ("
        << "avg_occ: " << avg_occ
        << ", alphabet size: " << alphabet_size
        << ", mode: " << (rm == gmr ? "gmr" : (
            rm == naive ? "naive" : "hybrid"))
        << "): " << std::flush;

    uint32_t dummy_var;
    auto time_start = now();

    for (uint32_t query=0; query<num_queries; query++) {
        if constexpr (rm == gmr) {
            dummy_var += rank_select.rank<true>(
                rank_queries[query].sym,
                rank_queries[query].pos
            );
        } else if constexpr (rm == naive) {
            dummy_var += rank_select.rank<false>(
                rank_queries[query].sym,
                rank_queries[query].pos
            );
        } else if constexpr (rm == hybrid) {
            dummy_var += rank_select.rank(
                rank_queries[query].sym,
                rank_queries[query].pos
            );
        }
    }

    auto time_end = now();
    std::to_string(dummy_var);

    std::cout << format_query_throughput(
        num_queries,time_diff_ns(time_start,time_end))
        << std::endl;
}

template <bool use_gmr>
void bench_select(
    rank_select_support<uint32_t>& rank_select,
    uint32_t avg_occ, uint32_t alphabet_size
) {
    std::cout << "benchmarking select queries ("
        << "avg_occ: " << avg_occ
        << ", alphabet size: " << alphabet_size
        << ", mode: " << (use_gmr ? "gmr" : "naive")
        << "): " << std::flush;

    uint32_t dummy_var;
    auto time_start = now();

    for (uint32_t query=0; query<num_queries; query++) {
        dummy_var += rank_select.select<use_gmr>(
            select_queries[query].sym,
            select_queries[query].rank
        );
    }

    auto time_end = now();
    std::to_string(dummy_var);

    std::cout << format_query_throughput(
        num_queries,time_diff_ns(time_start,time_end))
        << std::endl;
}

int main() {
    std::srand(std::time(0));
    omp_set_num_threads(1);
    rank_select_support<uint32_t> rank_select;

    for (uint32_t avg_occ=1; avg_occ<=target_input_size; avg_occ*=2) {
        uint32_t alphabet_size = target_input_size/avg_occ;
        std::uniform_int_distribution<uint32_t> occ_distrib(1,2*avg_occ);
        num_occ.resize(alphabet_size);
        
        for (uint32_t sym=0; sym<alphabet_size; sym++) {
            num_occ[sym] = occ_distrib(gen);

            for (uint32_t occ=0; occ<num_occ[sym]; occ++) {
                input.emplace_back(sym);
            }
        }

        std::shuffle(input.begin(),input.end(),gen);
        rank_select = rank_select_support<uint32_t>(input,alphabet_size);
        std::uniform_int_distribution<uint32_t> sym_distrib(0,alphabet_size-1);
        std::uniform_int_distribution<uint32_t> pos_distrib(1,input.size());

        for (uint32_t query=0; query<num_queries; query++) {
            uint32_t sym = sym_distrib(gen);
            rank_queries.emplace_back(sym,pos_distrib(gen));
            select_queries.emplace_back(sym,
                std::experimental::randint(uint32_t{1},num_occ[sym]));
        }

        bench_rank<gmr>(rank_select,avg_occ,alphabet_size);
        bench_rank<naive>(rank_select,avg_occ,alphabet_size);
        bench_rank<hybrid>(rank_select,avg_occ,alphabet_size);

        bench_select<true>(rank_select,avg_occ,alphabet_size);
        bench_select<false>(rank_select,avg_occ,alphabet_size);

        input.clear();
        num_occ.clear();
        rank_queries.clear();
        select_queries.clear();

        std::cout << std::endl;
    }
}