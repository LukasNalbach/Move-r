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

#ifndef MOVE_R_BENCH_MOVE_RB_BENCH_OPS_MOVE_RB_BENCH_OPS_TPP
#define MOVE_R_BENCH_MOVE_RB_BENCH_OPS_MOVE_RB_BENCH_OPS_TPP

#include <algorithm>
#include <functional>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <misc/index_cli.hpp>
#include <misc/move_rb_cli.hpp>
#include <move_rb/move_rb.hpp>

#include "../bench_ops.hpp"

namespace move_rb_bench_ops {

using bench_ops::locate_cell_t;

struct config_t {
    // number of contexts sampled per pattern length; each operation is measured once on each of them
    uint64_t num_samples = bench_ops::default_num_patterns;
    uint64_t seed = bench_ops::default_seed; // seed for sampling the contexts and for the order of the queries
};

/* the measured operations (see sample_t::context_for for the context an operation is measured on):
   - extend_<dir>: one extension in the direction dir of a context whose last symbol has been added in the
     direction dir
   - _switch: the last symbol of the context has been added in the other direction, i.e., the extension first
     updates the input intervals of its direction (see move_rb's valid_interval_dirs)
   - _all: enumeration of all extensions in the direction dir
   - _chain: chain_length consecutive extensions in the direction dir, each of the result of the previous one */
enum class query_kind_t : uint8_t {
    extend_left,
    extend_left_switch,
    extend_right,
    extend_right_switch,
    extend_left_all,
    extend_right_all,
    extend_left_chain,
    extend_right_chain
};

// number of operations in a chain (see query_kind_t; measurement setup)
static constexpr uint64_t chain_length = 16;

// largest time of a query in ns that is stored in cell_t::per_query (the largest value of its entries; larger times
// are stored as this value)
static constexpr double max_stored_query_time_ns = std::numeric_limits<uint32_t>::max();

// returns the name of an operation in the RESULT lines
inline std::string_view query_kind_name(query_kind_t kind)
{
    switch (kind) {
        case query_kind_t::extend_left:           return "extend-left";
        case query_kind_t::extend_left_switch:    return "extend-left-switch";
        case query_kind_t::extend_right:          return "extend-right";
        case query_kind_t::extend_right_switch:   return "extend-right-switch";
        case query_kind_t::extend_left_all:       return "extend-left-all";
        case query_kind_t::extend_right_all:      return "extend-right-all";
        case query_kind_t::extend_left_chain:     return "extend-left-chain";
        case query_kind_t::extend_right_chain:    return "extend-right-chain";
    }
    return "unknown";
}

// returns whether the time of a query is reported per enumerated extension (_all) resp. per operation (_chain)
inline bool reports_per_extension(query_kind_t kind)
{
    switch (kind) {
        case query_kind_t::extend_left_all:
        case query_kind_t::extend_right_all:
        case query_kind_t::extend_left_chain:
        case query_kind_t::extend_right_chain:
            return true;
        default:
            return false;
    }
}

// measurements of an operation for a pattern length
struct cell_t {
    uint64_t queries = 0; // number of measured queries
    uint64_t time_ns = 0; // total time in ns
    uint64_t extensions = 0; // number of enumerated extensions resp. operations (see reports_per_extension)
    // the times of the queries: [0..queries-1] time of the i-th query in ns, without the clock overhead (per
    // extension resp. operation, see reports_per_extension), at most 4 * 10^9
    std::vector<uint32_t> per_query;
};

// measurements of all operations for all pattern lengths
struct grid_t {
    bool locate_contexts = false; // true <=> the operations have been measured on LOCATE contexts
    // true <=> locate has been measured (LOCATE contexts and move_rb::supports_multiple_locate)
    bool locate_measured = false;
    std::vector<uint64_t> lengths; // the measured pattern lengths m
    std::vector<query_kind_t> kinds; // the measured operations
    // the measurements of the operations: [0..kinds.size() * lengths.size() - 1] cells[kind * lengths.size() +
    // length] belongs to the operation kinds[kind] and the pattern length lengths[length]
    std::vector<cell_t> cells;
    // the measurements of locate: [0..lengths.size()-1] for the pattern length lengths[i] (only with LOCATE contexts)
    std::vector<locate_cell_t> locate;
    uint64_t samples = 0; // number of sampled contexts per pattern length
    uint64_t checksum = 0; // hash of the results of the queries
    double clock_overhead_ns = 0; // time in ns of the two clock reads around a query

    // returns the measurements of the operation kinds[kind] for the pattern length lengths[length]
    const cell_t& at(uint64_t kind, uint64_t length) const
    { return cells[kind * lengths.size() + length]; }

    // returns the mean time of a query in ns, without the clock overhead (per extension resp. operation, see
    // reports_per_extension)
    double ns_per_query(uint64_t kind, uint64_t length) const
    {
        const cell_t& cell = at(kind, length);
        if (cell.queries == 0) return 0.0;
        double net = std::max(0.0, double(cell.time_ns) - double(cell.queries) * clock_overhead_ns);

        if (reports_per_extension(kinds[kind])) {
            return cell.extensions == 0 ? 0.0 : net / double(cell.extensions);
        }

        return net / double(cell.queries);
    }

    // returns the q-quantile (0 <= q <= 1) of the query times of a cell in ns (see cell_t::per_query)
    double quantile_ns(uint64_t kind, uint64_t length, double q) const
    {
        const cell_t& cell = at(kind, length);
        if (cell.per_query.empty()) return 0.0;
        std::vector<uint32_t> sorted = cell.per_query;
        uint64_t pos = std::min<uint64_t>(sorted.size() - 1, uint64_t(q * double(sorted.size())));
        std::nth_element(sorted.begin(), sorted.begin() + pos, sorted.end());
        return double(sorted[pos]);
    }

    double median_ns(uint64_t kind, uint64_t length) const { return quantile_ns(kind, length, 0.5); }
    double p90_ns(uint64_t kind, uint64_t length) const { return quantile_ns(kind, length, 0.9); }
    double max_ns(uint64_t kind, uint64_t length) const { return quantile_ns(kind, length, 1.0); }
};

namespace detail {

using bench_ops::next_random;

template <typename index_t, query_support_t query_support>
using context_t = typename index_t::template search_context_t<query_support>;

// the contexts and symbols the operations are measured with for one sampled pattern
template <typename index_t, query_support_t query_support>
struct sample_t {
    // two contexts of the same length; the last symbol of context has been added on the left, the last symbol
    // of context_r on the right
    context_t<index_t, query_support> context;
    context_t<index_t, query_support> context_r;
    uint64_t length = 0;
    uint8_t left_neighbor = 0, right_neighbor = 0;     // symbols context can be extended by
    uint8_t left_neighbor_r = 0, right_neighbor_r = 0; // symbols context_r can be extended by
    std::vector<uint8_t> chain_left;  // symbols of a chain of extensions of context to the left
    std::vector<uint8_t> chain_right; // symbols of a chain of extensions of context_r to the right

    // returns the context an operation is measured on
    const context_t<index_t, query_support>& context_for(query_kind_t kind) const
    {
        switch (kind) {
            case query_kind_t::extend_right:
            case query_kind_t::extend_right_all:
            case query_kind_t::extend_right_chain:
            case query_kind_t::extend_left_switch:
                return context_r;
            default:
                return context;
        }
    }
};

// the extensions of a context, in the order for_each_extension reports them
template <typename context_type>
struct children_t {
    std::vector<std::pair<uint64_t, context_type>> entries; // (symbol, extended context)
    uint64_t weight = 0;                                    // sum of the numbers of occurrences of the extended contexts

    void clear() { entries.clear(); weight = 0; }
};

// stores the extensions of node in direction dir in out
template <direction_t dir, typename context_type>
void list_children(context_type& node, children_t<context_type>& out)
{
    out.clear();

    node.template for_each_extension<dir, NO_TERM>([&](auto symbol, auto child) {
        uint64_t occurrences = child.num_occ();
        out.weight += occurrences;
        out.entries.emplace_back(uint64_t(symbol), std::move(child));
        return true;
    });
}

// returns a random extension, chosen with a probability proportional to its number of occurrences (nullptr if
// there is no extension)
template <typename context_type>
const std::pair<uint64_t, context_type>* draw_child(
    const children_t<context_type>& children, uint64_t& random_state)
{
    if (children.entries.empty() || children.weight == 0) return nullptr;
    uint64_t threshold = next_random(random_state) % children.weight;
    uint64_t running = 0;

    for (const auto& entry : children.entries) {
        running += entry.second.num_occ();
        if (threshold < running) return &entry;
    }

    return &children.entries.back();
}

// returns a random extension of from in direction dir (see draw_child); the extensions are enumerated on a
// copy of from
template <direction_t dir, typename context_type>
const std::pair<uint64_t, context_type>* draw_extension(
    const context_type& from, children_t<context_type>& children, uint64_t& random_state)
{
    context_type copy = from;
    list_children<dir>(copy, children);
    return draw_child(children, random_state);
}

// returns the symbols of up to steps consecutive random extensions of ctx in direction dir (see draw_child)
template <direction_t dir, typename context_type>
std::vector<uint8_t> record_chain(context_type ctx, uint64_t steps, children_t<context_type>& children,
    uint64_t& random_state)
{
    std::vector<uint8_t> symbols;

    for (uint64_t t = 0; t < steps; t++) {
        list_children<dir>(ctx, children);
        const auto* chosen = draw_child(children, random_state);
        if (chosen == nullptr) break;
        symbols.push_back(uint8_t(chosen->first));
        ctx = chosen->second;
    }

    return symbols;
}

/**
 * @brief samples a context of a given length from the index by extending the empty context with random
 *        symbols (see draw_child)
 * @param index the index
 * @param length length of the context
 * @param random_state state of the random number generator
 * @param ctx the sampled context
 * @return whether the context has the given length
 */
template <query_support_t query_support, typename index_t>
bool sample_context(const index_t& index, uint64_t length, uint64_t& random_state,
    context_t<index_t, query_support>& ctx)
{
    using ctx_t = context_t<index_t, query_support>;

    children_t<ctx_t> children;
    ctx = index.template empty_context<query_support>();

    for (uint64_t depth = 0; depth < length; depth++) {
        // the context is extended alternately to the right and to the left
        if (depth % 2 == 0) list_children<RIGHT>(ctx, children);
        else list_children<LEFT>(ctx, children);

        const auto* chosen = draw_child(children, random_state);
        if (chosen == nullptr) return false;

        ctx = chosen->second;
    }

    return ctx.pattern_length() == length;
}

/**
 * @brief samples the contexts and the symbols that the operations are measured with for a pattern length
 * @param index the index
 * @param length pattern length
 * @param random_state state of the random number generator
 * @return the sample (see sample_t)
 */
template <query_support_t query_support, typename index_t>
sample_t<index_t, query_support> make_sample(
    const index_t& index, uint64_t length, uint64_t& random_state)
{
    using ctx_t = context_t<index_t, query_support>;
    using sample_type = sample_t<index_t, query_support>;

    sample_type sample;
    sample.length = length;
    children_t<ctx_t> children;

    // sample a context of length length - 1 and extend it by one symbol to the left (context) and to the
    // right (context_r)
    ctx_t prefix = index.template empty_context<query_support>();
    sample_context<query_support>(index, length - 1, random_state, prefix);
    sample.context = prefix;
    sample.context_r = prefix;

    if (const auto* c = draw_extension<LEFT>(prefix, children, random_state)) sample.context = c->second;
    if (const auto* c = draw_extension<RIGHT>(prefix, children, random_state)) sample.context_r = c->second;

    // choose the symbols of the measured extensions among the symbols the contexts can be extended by
    if (const auto* c = draw_extension<LEFT>(sample.context, children, random_state))
        sample.left_neighbor = uint8_t(c->first);
    if (const auto* c = draw_extension<LEFT>(sample.context_r, children, random_state))
        sample.left_neighbor_r = uint8_t(c->first);
    sample.chain_left = record_chain<LEFT>(sample.context, chain_length, children, random_state);

    if (const auto* c = draw_extension<RIGHT>(sample.context, children, random_state))
        sample.right_neighbor = uint8_t(c->first);
    if (const auto* c = draw_extension<RIGHT>(sample.context_r, children, random_state))
        sample.right_neighbor_r = uint8_t(c->first);
    sample.chain_right = record_chain<RIGHT>(sample.context_r, chain_length, children, random_state);

    return sample;
}

/**
 * @brief creates the samples of a pattern length (see make_sample)
 * @param index the index
 * @param length pattern length
 * @param requested number of samples
 * @param seed seed of the random number generator
 */
template <query_support_t query_support, typename index_t>
std::vector<sample_t<index_t, query_support>> make_samples(
    const index_t& index, uint64_t length, uint64_t requested, uint64_t seed)
{
    uint64_t state = seed == 0 ? 1 : seed;
    std::vector<sample_t<index_t, query_support>> samples;
    samples.reserve(requested);

    for (uint64_t i = 0; i < requested; i++) {
        samples.emplace_back(make_sample<query_support>(index, length, state));
    }

    return samples;
}

// returns the operations that are measured (see query_kind_t)
inline std::vector<query_kind_t> supported_kinds()
{
    return { query_kind_t::extend_left, query_kind_t::extend_left_switch, query_kind_t::extend_left_all,
        query_kind_t::extend_left_chain, query_kind_t::extend_right, query_kind_t::extend_right_switch,
        query_kind_t::extend_right_all, query_kind_t::extend_right_chain };
}

// returns a hash of the forward suffix array interval and the pattern length of a context
template <typename context_type>
uint64_t context_checksum(const context_type& context)
{
    auto [begin, end] = context.forward_sa_interval();
    return mix_64(begin) ^ (mix_64(end) << 1) ^ (mix_64(context.pattern_length()) << 2);
}

/**
 * @brief executes an operation
 * @param kind the operation
 * @param context copy of sample.context_for(kind)
 * @param sample the sample
 * @param extensions the number of enumerated extensions (_all) resp. executed operations (_chain) is stored in it
 * @return checksum of the result
 */
template <query_support_t query_support, typename index_t>
uint64_t execute(
    query_kind_t kind,
    context_t<index_t, query_support>& context,
    const sample_t<index_t, query_support>& sample,
    uint64_t& extensions)
{
    switch (kind) {
        case query_kind_t::extend_left:
        case query_kind_t::extend_left_switch: {
            const uint8_t symbol = kind == query_kind_t::extend_left ? sample.left_neighbor : sample.left_neighbor_r;
            auto [extended, occurs] = context.template extend<LEFT>(symbol);
            return context_checksum(extended) ^ uint64_t(occurs);
        }

        case query_kind_t::extend_right:
        case query_kind_t::extend_right_switch: {
            const uint8_t symbol = kind == query_kind_t::extend_right_switch
                ? sample.right_neighbor : sample.right_neighbor_r;
            auto [extended, occurs] = context.template extend<RIGHT>(symbol);
            return context_checksum(extended) ^ (uint64_t(occurs) << 1);
        }

        case query_kind_t::extend_left_chain:
            for (uint8_t symbol : sample.chain_left) {
                auto [extended, occurs] = context.template extend<LEFT>(symbol);
                context = std::move(extended);
            }

            extensions = sample.chain_left.size();
            return context_checksum(context);

        case query_kind_t::extend_right_chain:
            for (uint8_t symbol : sample.chain_right) {
                auto [extended, occurs] = context.template extend<RIGHT>(symbol);
                context = std::move(extended);
            }

            extensions = sample.chain_right.size();
            return context_checksum(context);

        case query_kind_t::extend_left_all:
        case query_kind_t::extend_right_all: {
            uint64_t checksum = 0;

            auto report = [&](uint64_t symbol, const auto& child) {
                checksum ^= mix_64(uint64_t(symbol) + context_checksum(child));
                extensions++;
                return true;
            };

            if (kind == query_kind_t::extend_left_all) context.template for_each_extension<LEFT, NO_TERM>(report);
            else context.template for_each_extension<RIGHT, NO_TERM>(report);
            return checksum;
        }
    }

    throw std::logic_error("operation kind is not supported by this index");
}

/**
 * @brief measures locate on the contexts of a pattern length until the time budget is used up. Contexts with
 *        more than bench_ops::max_occurrences_per_call occurrences are skipped; if all contexts have more
 *        occurrences, max_occurrences_per_call occurrences per context are located with next_occ() instead.
 * @param samples the samples of the pattern length
 * @param target_time_ns time budget in ns
 * @param seed seed of the random number generator
 * @return the measurements
 */
template <query_support_t query_support, typename index_t>
locate_cell_t measure_locate(
    const std::vector<sample_t<index_t, query_support>>& samples,
    uint64_t target_time_ns, uint64_t seed)
{
    locate_cell_t result;
    if constexpr (query_support != LOCATE) {
        return result;
    } else {
        using bench_ops::max_occurrences_per_call;
        uint64_t state = seed == 0 ? 1 : seed;
        uint64_t cursor = next_random(state) % samples.size();
        const uint64_t budget = std::max<uint64_t>(1, target_time_ns);

        // only the calls of locate() are measured; the budget is checked after each call
        uint64_t skipped = 0;
        do {
            const auto& sample = samples[cursor++ % samples.size()];
            auto locate = sample.context.locate_phase();

            if (locate.num_occ_rem() > max_occurrences_per_call) {
                // stop if all contexts have more than max_occurrences_per_call occurrences
                if (++skipped >= samples.size()) break;
                continue;
            }
            skipped = 0;

            result.calls++;
            auto begin = now();
            locate.locate([&](uint64_t occurrence) {
                result.checksum ^= mix_64(occurrence + result.occurrences);
                result.occurrences++;
            });
            result.time_ns += time_diff_ns(begin, now());
        } while (result.time_ns < budget);

        // no context has been located: locate at most max_occurrences_per_call occurrences per context with next_occ()
        if (result.calls == 0) {
            do {
                const auto& sample = samples[cursor++ % samples.size()];
                auto locate = sample.context.locate_phase();
                uint64_t remaining = std::min(locate.num_occ_rem(), max_occurrences_per_call);
                if (remaining == 0) continue;

                result.calls++;
                auto begin = now();

                for (uint64_t k = 0; k < remaining; k++) {
                    result.checksum ^= mix_64(locate.next_occ() + result.occurrences);
                    result.occurrences++;
                }

                result.time_ns += time_diff_ns(begin, now());
            } while (result.time_ns < budget && result.calls < samples.size());
        }

        return result;
    }
}

} // namespace detail

/**
 * @brief measures all operations supported by an index for all given pattern lengths. Each operation is
 *        executed once on each sample of each pattern length; all queries are executed in random order and
 *        measured individually. Locate is measured separately for each pattern length (see measure_locate).
 * @tparam query_support COUNT or LOCATE contexts
 * @param index the index
 * @param lengths pattern lengths
 * @param config number of samples per pattern length and seed
 * @param report_progress function called with each pattern length before its samples are created
 * @return the measurements
 */
template <query_support_t query_support, typename index_t>
grid_t benchmark(
    const index_t& index,
    const std::vector<uint64_t>& lengths,
    const config_t& config,
    const std::function<void(uint64_t)>& report_progress)
{
    using sample_type = detail::sample_t<index_t, query_support>;

    grid_t grid;
    grid.locate_contexts = query_support == LOCATE;
    grid.locate_measured = query_support == LOCATE && index_t::supports_multiple_locate;
    grid.lengths = lengths;
    grid.kinds = detail::supported_kinds();
    grid.samples = config.num_samples;

    if (grid.kinds.empty()) throw std::logic_error("the index supports no measurable operation");
    grid.cells.resize(grid.kinds.size() * lengths.size());
    grid.locate.resize(lengths.size());

    std::vector<std::vector<sample_type>> samples_by_length;
    samples_by_length.reserve(lengths.size());

    for (uint64_t length : lengths) {
        report_progress(length);
        samples_by_length.emplace_back(detail::make_samples<query_support>(
            index, length, config.num_samples,
            config.seed ^ mix_64(length) ^
            (query_support == LOCATE ? golden_ratio_64 : 0)));
    }

    std::mt19937_64 generator(config.seed ^ (query_support == LOCATE ? bench_ops::schedule_seed_salt : 0));
    std::vector<bench_ops::schedule_entry_t> schedule;
    schedule.reserve(grid.kinds.size() * lengths.size() * config.num_samples);

    // one query per operation, pattern length and sample
    for (uint32_t kind = 0; kind < grid.kinds.size(); kind++) {
        for (uint32_t length = 0; length < lengths.size(); length++) {
            for (uint32_t sample = 0; sample < samples_by_length[length].size(); sample++) {
                schedule.push_back({ kind, length, sample });
            }
        }
    }

    std::shuffle(schedule.begin(), schedule.end(), generator);
    grid.clock_overhead_ns = bench_ops::measure_clock_overhead();

    for (const bench_ops::schedule_entry_t& entry : schedule) {
        const sample_type& sample = samples_by_length[entry.length][entry.sample];

        // the context is copied before the time measurement starts
        auto context = sample.context_for(grid.kinds[entry.kind]);
        uint64_t extensions = 0;
        auto begin = now();
        uint64_t result = detail::execute<query_support>(
            grid.kinds[entry.kind], context,
            sample, extensions);
        uint64_t elapsed = time_diff_ns(begin, now());

        cell_t& cell = grid.cells[entry.kind * lengths.size() + entry.length];
        cell.queries++;
        cell.time_ns += elapsed;
        cell.extensions += extensions;
        grid.checksum ^= result + cell.queries;

        double net = std::max(0.0, double(elapsed) - grid.clock_overhead_ns);
        if (reports_per_extension(grid.kinds[entry.kind])) net = extensions == 0 ? 0.0 : net / double(extensions);
        cell.per_query.push_back(uint32_t(std::min(net, max_stored_query_time_ns)));
    }

    // the time budget of locate is the average time the operations took in total at the pattern length
    if constexpr (query_support == LOCATE && index_t::supports_multiple_locate) {
        for (uint64_t length = 0; length < lengths.size(); length++) {
            grid.locate[length] = detail::measure_locate<query_support>(samples_by_length[length],
                bench_ops::mean_operation_time(grid, length), config.seed ^ mix_64(lengths[length] + 1));
        }
    }

    return grid;
}

} // namespace move_rb_bench_ops

#endif // MOVE_R_BENCH_MOVE_RB_BENCH_OPS_MOVE_RB_BENCH_OPS_TPP
