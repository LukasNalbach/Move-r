/**
 * part of LukasNalbach/Move-r
 *
 * MIT License
 *
 * Copyright (c) Lukas Nalbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef MOVE_R_BENCH_MOVE_R_BENCH_OPS_MOVE_R_BENCH_OPS_TPP
#define MOVE_R_BENCH_MOVE_R_BENCH_OPS_MOVE_R_BENCH_OPS_TPP

#include <algorithm>
#include <cstdint>
#include <functional>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <misc/index_cli.hpp>
#include <move_r/move_r.hpp>

#include "../bench_ops.hpp"

namespace move_r_bench_ops {

using bench_ops::locate_cell_t;

// draw_symbol draws at most this many positions to find a symbol other than the terminator (measurement setup)
inline constexpr uint64_t max_symbol_draws = 8;
// benchmark_access: the number of unmeasured warm-up queries of BWT(i) and of SA(i), and the number of SA(i) queries
// between two checks of the time budget
inline constexpr uint64_t bwt_warmup_queries = 4096;
inline constexpr uint64_t sa_warmup_queries = 64;
inline constexpr uint64_t sa_budget_check_period = 64;

struct config_t {
    // number of patterns per pattern length; each operation is run once per pattern
    uint64_t num_samples = bench_ops::default_num_patterns;
    uint64_t seed = bench_ops::default_seed; // seed for sampling the patterns and for the order of the queries
};

// the operations measured per pattern length
enum class query_kind_t : uint8_t {
    prepend,  // one backward-search step
    one_occ,  // computing one occurrence of the matched pattern
    next_occ  // computing the second occurrence of the matched pattern (the second call of next_occ())
};

// returns the name of an operation in the RESULT lines
inline std::string_view query_kind_name(query_kind_t kind)
{
    switch (kind) {
        case query_kind_t::prepend:  return "prepend";
        case query_kind_t::one_occ:  return "one-occ";
        case query_kind_t::next_occ: return "next-occ";
    }
    return "unknown";
}

// number of queries and total time of one operation for one pattern length
struct cell_t {
    uint64_t queries = 0;
    uint64_t time_ns = 0;
};

// results of all operations for all pattern lengths
struct grid_t {
    std::vector<uint64_t> lengths; // the measured pattern lengths m
    std::vector<query_kind_t> kinds; // the measured operations
    // the results of the operations: [0..kinds.size() * lengths.size() - 1] cells[kind * lengths.size() + length]
    // is the result of the operation kinds[kind] for the pattern length lengths[length]
    std::vector<cell_t> cells;
    std::vector<locate_cell_t> locate; // the results of locate: [0..lengths.size()-1] for the pattern length lengths[i]
    uint64_t samples = 0; // number of sampled patterns per pattern length
    uint64_t checksum = 0; // hash of the results of the queries
    double clock_overhead_ns = 0; // time in ns of the two clock reads around a query

    // returns the result of the operation kinds[kind] for the pattern length lengths[length]
    const cell_t& at(uint64_t kind, uint64_t length) const
    { return cells[kind * lengths.size() + length]; }

    // returns the average time per query in ns without the clock overhead (at least 0)
    double ns_per_query(uint64_t kind, uint64_t length) const
    {
        const cell_t& cell = at(kind, length);
        if (cell.queries == 0) return 0.0;
        return std::max(0.0, double(cell.time_ns) - double(cell.queries) * clock_overhead_ns)
            / double(cell.queries);
    }
};

// result of the SA(i) or BWT(i) queries at uniformly random positions i
struct access_cell_t {
    uint64_t queries = 0;
    uint64_t time_ns = 0;
    uint64_t checksum = 0;

    double ns_per_query() const { return queries == 0 ? 0.0 : double(time_ns) / double(queries); }
};

// results of the SA(i) and BWT(i) queries
struct access_grid_t {
    bool has_sa = false;
    access_cell_t sa;
    access_cell_t bwt;
};

namespace detail {

using bench_ops::next_random;

template <typename index_t>
using context_t = typename index_t::query_context_t;

// a sampled pattern P and the symbol prepend is measured with
template <typename index_t>
struct sample_t {
    context_t<index_t> context; // query context of P
    uint64_t length = 0; // the requested length of P (context.length() is smaller if P could not be sampled)
    uint64_t occurrences = 0; // number of occurrences of P
    uint64_t symbol = 0; // a symbol c, s.t. cP occurs in the input (used to measure prepend)
    bool has_symbol = false; // true <=> symbol could be drawn
};

/**
 * @brief draws the symbol L[i] at a uniformly random position i in the suffix array interval of the matched
 *        pattern; the terminator symbol is not returned
 * @param index the index
 * @param context query context of the matched pattern; it is modified by sa_interval(), pass a copy to keep
 *        the state after the backward search
 * @param random_state state of the random number generator
 * @param symbol the drawn symbol (valid only if true is returned)
 * @return whether a symbol was drawn
 */
template <typename index_t>
bool draw_symbol(const index_t& index, context_t<index_t>& context, uint64_t& random_state,
    uint64_t& symbol)
{
    auto [begin, end] = context.sa_interval();
    if (end < begin) return false;

    // up to max_symbol_draws attempts to draw a symbol other than the terminator symbol
    for (uint64_t attempt = 0; attempt < max_symbol_draws; attempt++) {
        const uint64_t i = begin + next_random(random_state) % (end - begin + 1);
        symbol = index.BWT(i);
        if (symbol != 0) return true;
    }

    return false;
}

/**
 * @brief samples a pattern of the given length by prepending symbols drawn with draw_symbol
 * @param index the index
 * @param length length of the pattern
 * @param random_state state of the random number generator
 * @param ctx query context of the sampled pattern
 * @return length of the sampled pattern (smaller than length if no symbol could be drawn)
 */
template <typename index_t>
uint64_t sample_pattern(const index_t& index, uint64_t length, uint64_t& random_state,
    context_t<index_t>& ctx)
{
    ctx = index.query();

    for (uint64_t i = 0; i < length; i++) {
        uint64_t symbol = 0;
        context_t<index_t> copy = ctx; // draw_symbol modifies the context
        if (!draw_symbol(index, copy, random_state, symbol)) return i;
        if (!ctx.prepend(symbol)) return i;
    }

    return length;
}

// samples a pattern of the given length and a symbol to prepend to it
template <typename index_t>
sample_t<index_t> make_sample(const index_t& index, uint64_t length, uint64_t& random_state)
{
    sample_t<index_t> sample { .context = index.query(), .length = length };
    if (sample_pattern(index, length, random_state, sample.context) != length) return sample;

    context_t<index_t> copy = sample.context; // draw_symbol modifies the context
    sample.occurrences = copy.num_occ();
    sample.has_symbol = draw_symbol(index, copy, random_state, sample.symbol);
    return sample;
}

// samples requested patterns of one pattern length (see make_sample)
template <typename index_t>
std::vector<sample_t<index_t>> make_samples(
    const index_t& index, uint64_t length, uint64_t requested, uint64_t seed)
{
    uint64_t state = seed == 0 ? 1 : seed;
    std::vector<sample_t<index_t>> samples;
    samples.reserve(requested);

    for (uint64_t i = 0; i < requested; i++)
        samples.emplace_back(make_sample(index, length, state));

    return samples;
}

// returns the operations the index supports
template <typename index_t>
std::vector<query_kind_t> supported_kinds()
{
    std::vector<query_kind_t> kinds { query_kind_t::prepend };

    if constexpr (index_t::supports_locate) kinds.emplace_back(query_kind_t::one_occ);
    if constexpr (index_t::supports_multiple_locate) kinds.emplace_back(query_kind_t::next_occ);

    return kinds;
}

// executes an operation on the context of a sample and returns its result
template <typename index_t>
uint64_t execute(query_kind_t kind, context_t<index_t>& context, const sample_t<index_t>& sample)
{
    switch (kind) {
        case query_kind_t::prepend:
            return uint64_t(context.prepend(sample.symbol));

        case query_kind_t::one_occ:
            if constexpr (index_t::supports_locate) return context.one_occ();
            break;

        case query_kind_t::next_occ:
            if constexpr (index_t::supports_multiple_locate) return context.next_occ();
            break;
    }

    throw std::logic_error("operation kind is not supported by this index");
}

// returns whether next_occ() is called once on the context before the operation is measured
inline bool requires_first_occ(query_kind_t kind)
{
    return kind == query_kind_t::next_occ;
}

// returns whether an operation can be measured with a sample: the pattern must have the requested length,
// prepend requires a drawn symbol and next_occ at least two occurrences
template <typename index_t>
bool is_usable(query_kind_t kind, const sample_t<index_t>& sample)
{
    if (sample.context.length() != sample.length) return false;
    if (kind == query_kind_t::prepend) return sample.has_symbol;
    if (kind == query_kind_t::next_occ) return sample.occurrences >= 2;
    return true;
}

// measures locate on the patterns of one pattern length with at most bench_ops::max_occurrences_per_call
// occurrences, starting with a random pattern, until target_time_ns has passed
template <typename index_t>
locate_cell_t measure_locate(const std::vector<sample_t<index_t>>& samples,
    uint64_t target_time_ns, uint64_t seed)
{
    locate_cell_t result;

    if constexpr (!index_t::supports_multiple_locate) {
        return result;
    } else {
        if (samples.empty()) return result;
        uint64_t state = seed == 0 ? 1 : seed;
        uint64_t cursor = next_random(state) % samples.size();
        uint64_t skipped = 0;

        do {
            const sample_t<index_t>& sample = samples[cursor++ % samples.size()];

            if (sample.context.length() != sample.length ||
                sample.occurrences > bench_ops::max_occurrences_per_call) {
                // stop if all patterns were skipped
                if (++skipped >= samples.size()) break;
                continue;
            }

            skipped = 0;
            auto context = sample.context;
            result.calls++;
            auto begin = now();
            context.locate([&](uint64_t, uint64_t occurrence) {
                result.checksum ^= mix_64(occurrence + result.occurrences);
                result.occurrences++;
            });
            result.time_ns += time_diff_ns(begin, now());
        } while (result.time_ns < std::max<uint64_t>(1, target_time_ns));

        return result;
    }
}

} // namespace detail

/**
 * @brief measures the operations supported by the index: each operation is run once on each pattern of each
 *        pattern length, in random order, and each query is timed separately; locate is measured separately
 *        per pattern length (see measure_locate)
 * @param index the index
 * @param lengths the pattern lengths
 * @param config number of patterns per pattern length and seed
 * @param report_progress function called with each pattern length before its patterns are sampled
 * @return the results
 */
template <typename index_t>
grid_t benchmark(const index_t& index, const std::vector<uint64_t>& lengths, const config_t& config,
    const std::function<void(uint64_t)>& report_progress)
{
    using sample_type = detail::sample_t<index_t>;

    grid_t grid;
    grid.lengths = lengths;
    grid.kinds = detail::supported_kinds<index_t>();
    grid.samples = config.num_samples;

    if (grid.kinds.empty()) throw std::logic_error("the index supports no measurable operation");
    grid.cells.resize(grid.kinds.size() * lengths.size());
    grid.locate.resize(lengths.size());

    std::vector<std::vector<sample_type>> samples_by_length;
    samples_by_length.reserve(lengths.size());

    for (uint64_t length : lengths) {
        report_progress(length);
        samples_by_length.emplace_back(detail::make_samples(
            index, length, config.num_samples, config.seed ^ mix_64(length)));
    }

    // the queries to measure, in random order
    std::mt19937_64 generator(config.seed);
    std::vector<bench_ops::schedule_entry_t> schedule;
    schedule.reserve(grid.kinds.size() * lengths.size() * config.num_samples);

    for (uint32_t kind = 0; kind < grid.kinds.size(); kind++) {
        for (uint32_t length = 0; length < lengths.size(); length++) {
            for (uint32_t sample = 0; sample < samples_by_length[length].size(); sample++) {
                if (!detail::is_usable(grid.kinds[kind], samples_by_length[length][sample])) continue;
                schedule.push_back({ kind, length, sample });
            }
        }
    }

    std::shuffle(schedule.begin(), schedule.end(), generator);
    grid.clock_overhead_ns = bench_ops::measure_clock_overhead();

    for (const bench_ops::schedule_entry_t& entry : schedule) {
        const sample_type& sample = samples_by_length[entry.length][entry.sample];

        // prepare the query (not measured)
        auto context = sample.context;
        uint64_t prepared = 0;

        if constexpr (index_t::supports_multiple_locate) {
            if (detail::requires_first_occ(grid.kinds[entry.kind]))
                prepared = context.next_occ();
        }

        auto begin = now();
        uint64_t result = detail::execute(grid.kinds[entry.kind], context, sample);
        uint64_t elapsed = time_diff_ns(begin, now());

        cell_t& cell = grid.cells[entry.kind * lengths.size() + entry.length];
        cell.queries++;
        cell.time_ns += elapsed;
        grid.checksum ^= result + prepared + cell.queries;
    }

    // measure locate per pattern length for the average time of the other operations
    if constexpr (index_t::supports_multiple_locate) {
        for (uint64_t length = 0; length < lengths.size(); length++) {
            grid.locate[length] = detail::measure_locate(samples_by_length[length],
                bench_ops::mean_operation_time(grid, length), config.seed ^ mix_64(lengths[length] + 1));
        }
    }

    return grid;
}

/**
 * @brief measures BWT(i) and (if supported) SA(i) at uniformly random positions i
 * @param index the index
 * @param config seed
 * @param queries number of positions
 * @param sa_budget_ns time in ns after which the SA(i) queries are stopped
 * @return the results
 */
template <typename index_t>
access_grid_t benchmark_access(const index_t& index, const config_t& config, uint64_t queries,
    uint64_t sa_budget_ns)
{
    access_grid_t grid;
    grid.has_sa = index_t::has_SA_access;

    std::mt19937_64 generator(config.seed ^ bench_ops::access_seed_salt);
    std::uniform_int_distribution<uint64_t> pos_distrib(0, index.input_size());
    std::vector<uint64_t> positions(queries);
    for (uint64_t& pos : positions) pos = pos_distrib(generator);

    // warm-up queries (not measured)
    for (uint64_t i = 0; i < std::min<uint64_t>(bwt_warmup_queries, queries); i++)
        grid.bwt.checksum += index.BWT(positions[i]);

    auto begin = now();
    for (uint64_t pos : positions) grid.bwt.checksum += index.BWT(pos);
    grid.bwt.time_ns = time_diff_ns(begin, now());
    grid.bwt.queries = queries;

    if constexpr (index_t::has_SA_access) {
        for (uint64_t i = 0; i < std::min<uint64_t>(sa_warmup_queries, queries); i++)
            grid.sa.checksum += index.SA(positions[i]);

        begin = now();

        for (uint64_t pos : positions) {
            grid.sa.checksum += index.SA(pos);
            grid.sa.queries++;
            // check the time budget every sa_budget_check_period queries
            if (grid.sa.queries % sa_budget_check_period == 0 &&
                time_diff_ns(begin, now()) >= sa_budget_ns) break;
        }

        grid.sa.time_ns = time_diff_ns(begin, now());
    }

    return grid;
}

} // namespace move_r_bench_ops

#endif // MOVE_R_BENCH_MOVE_R_BENCH_OPS_MOVE_R_BENCH_OPS_TPP
