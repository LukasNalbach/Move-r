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

#ifndef MOVE_R_ALGORITHMS_MOVE_RB_APM_MODEL_TPP
#define MOVE_R_ALGORITHMS_MOVE_RB_APM_MODEL_TPP

#include <algorithms/move_rb/apm/model.hpp>

inline std::pair<uint64_t, uint64_t> apm_range_t::absolute(uint64_t m) const
{
    auto abs = [m](int64_t x) -> uint64_t {
        if (x < 0) return x + int64_t(m) < 0 ? 0 : uint64_t(x + int64_t(m));
        return std::min<uint64_t>(uint64_t(x), m);
    };

    uint64_t b = abs(beg);
    uint64_t e = abs(end);
    return {b, std::max(b, e)};
}

inline uint64_t parse_apm_value(std::string_view str, const std::string& what)
{
    if (str == "inf") return apm_forbidden;
    if (str.empty() || str.find_first_not_of("0123456789") != std::string_view::npos)
        throw std::invalid_argument("invalid " + what + " '" + std::string(str) + "'");
    return std::stoull(std::string(str));
}

// parses a subset of the letters s, i, d, o into a mask of apm_op_t
inline uint8_t parse_apm_ops(std::string_view str)
{
    uint8_t ops = 0;

    for (char c : str) {
        if (c == 's') ops |= APM_SUB;
        else if (c == 'i') ops |= APM_INS;
        else if (c == 'd') ops |= APM_DEL;
        else if (c == 'o') ops |= APM_OPEN;
        else throw std::invalid_argument("invalid operation '" + std::string(1, c) + "' (allowed: s, i, d, o)");
    }

    if (ops == 0) throw std::invalid_argument("no operation given");
    return ops;
}

inline apm_range_t parse_apm_range(std::string_view str)
{
    auto parse_int = [&](std::string_view s) -> int64_t {
        bool neg = !s.empty() && s[0] == '-';
        std::string_view digits = neg ? s.substr(1) : s;
        if (digits.empty() || digits.find_first_not_of("0123456789") != std::string_view::npos)
            throw std::invalid_argument("invalid range '" + std::string(str) + "'");
        int64_t v = std::stoll(std::string(digits));
        return neg ? -v : v;
    };

    apm_range_t range;
    uint64_t colon = str.find(':');

    if (colon == std::string_view::npos) {
        range.beg = parse_int(str);
        range.end = range.beg == -1 ? std::numeric_limits<int64_t>::max() : range.beg + 1;
        return range;
    }

    if (colon > 0) range.beg = parse_int(str.substr(0, colon));
    if (colon + 1 < str.size()) range.end = parse_int(str.substr(colon + 1));
    return range;
}

// splits "<head>[@<range>]" into head and range
inline std::string_view split_apm_range(std::string_view str, apm_range_t& range)
{
    uint64_t at = str.find('@');
    if (at == std::string_view::npos) return str;
    range = parse_apm_range(str.substr(at + 1));
    return str.substr(0, at);
}

inline apm_cost_rule_t parse_apm_cost_rule(std::string_view str)
{
    apm_cost_rule_t rule;
    uint64_t slash = str.find('/');

    if (slash != std::string_view::npos) {
        std::string_view syms = str.substr(slash + 1);
        uint64_t comma = syms.find(',');
        if (comma == std::string_view::npos) throw std::invalid_argument("missing ',' in '" + std::string(syms) + "'");
        rule.from = std::string(syms.substr(0, comma));
        rule.to = std::string(syms.substr(comma + 1));
        str = str.substr(0, slash);
    }

    str = split_apm_range(str, rule.range);
    uint64_t eq = str.find('=');
    if (eq == std::string_view::npos) throw std::invalid_argument("missing '=' in '" + std::string(str) + "'");
    rule.ops = parse_apm_ops(str.substr(0, eq));
    uint64_t cost = parse_apm_value(str.substr(eq + 1), "cost");
    if (cost > apm_forbidden) throw std::invalid_argument("cost too large");
    rule.cost = uint16_t(cost);
    return rule;
}

inline apm_counter_t parse_apm_counter(std::string_view str)
{
    apm_counter_t counter;
    str = split_apm_range(str, counter.range);
    uint64_t eq = str.find('=');
    if (eq == std::string_view::npos) throw std::invalid_argument("missing '=' in '" + std::string(str) + "'");
    counter.ops = parse_apm_ops(str.substr(0, eq));
    uint64_t max = parse_apm_value(str.substr(eq + 1), "counter bound");
    if (max > apm_max_counter_bound)
        throw std::invalid_argument("counter bound must be at most " + std::to_string(apm_max_counter_bound));
    counter.max = uint8_t(max);
    return counter;
}

inline apm_class_t parse_apm_class(std::string_view str)
{
    apm_class_t cls;
    str = split_apm_range(str, cls.range);
    uint64_t comma = str.find(',');
    if (comma == std::string_view::npos) throw std::invalid_argument("missing ',' in '" + std::string(str) + "'");
    cls.from = std::string(str.substr(0, comma));
    cls.to = std::string(str.substr(comma + 1));
    return cls;
}

template <std::ranges::contiguous_range pat_t>
void apm_costs_t::compile(const apm_model_t& model, const pat_t& P, std::string_view quality, bool reverse_strand,
    uint16_t bound, distance_metric_t metric, const std::vector<uint8_t>& alph_map)
{
    static_assert(sizeof(std::ranges::range_value_t<pat_t>) == 1);
    _m = std::ranges::size(P);
    _bound = bound;
    _tau = model.tau;
    _metric = metric;
    _num_counters = 0;
    _counter_bias = 0;
    _affine = false;

    if (bound >= apm_forbidden) throw std::invalid_argument("bound too large");
    if (model.quality && quality.size() != _m) throw std::invalid_argument("the model needs one quality per symbol");
    if (model.counters.size() > apm_max_counters)
        throw std::invalid_argument("at most " + std::to_string(apm_max_counters) + " counters");

    // the effective alphabet of the index
    if (&alph_map != _alph_map) {
        _alph_map = &alph_map;
        _num_ranks = 0;
        _rank.assign(256, 0);
        for (uint64_t c = 0; c < 256; c++) _num_ranks = std::max<uint64_t>(_num_ranks, alph_map[c]);
        _unrank.assign(_num_ranks, 0);

        for (uint64_t c = 0; c < 256; c++) {
            _rank[c] = alph_map[c] == 0 ? uint8_t(_num_ranks) : uint8_t(alph_map[c] - 1);
            if (alph_map[c] != 0) _unrank[alph_map[c] - 1] = uint8_t(c);
        }
    }

    // a range, a symbol and a quality of the read, as they are seen at the positions of P
    auto range_of = [&](const apm_range_t& range) -> std::pair<uint64_t, uint64_t> {
        auto [b, e] = range.absolute(_m);
        if (reverse_strand) return {_m - e, _m - b};
        return {b, e};
    };

    auto symbols_of = [&](const std::string& syms) { return reverse_strand ? reverse_complement(syms) : syms; };

    auto contains = [](const std::string& set, uint8_t c) { return set.empty() || set.find(char(c)) != std::string::npos; };
    _pattern.resize(_m);
    for (uint64_t i = 0; i < _m; i++) _pattern[i] = std::bit_cast<uint8_t>(std::ranges::data(P)[i]);
    auto sym = [&](uint64_t i) { return _pattern[i]; };

    _matches.assign(_m * _num_ranks, 0);
    _sub.assign(_m * _num_ranks, 1);
    _ins.assign(_m, 1);
    _del.assign(_m + 1, 1);
    _ins_open.assign(_m, 0);
    _del_open.assign(_m + 1, 0);
    _sub_counts.assign(_m, 0);
    _ins_counts.assign(_m, 0);
    _del_counts.assign(_m + 1, 0);
    _ins_open_counts.assign(_m, 0);
    _del_open_counts.assign(_m + 1, 0);

    for (uint64_t i = 0; i < _m; i++) {
        if (_rank[sym(i)] < _num_ranks) _matches[i * _num_ranks + _rank[sym(i)]] = 1;

        if (model.iupac) {
            for (char base : iupac_bases(sym(i))) {
                if (_rank[uint8_t(base)] < _num_ranks) _matches[i * _num_ranks + _rank[uint8_t(base)]] = 1;
            }
        }

        if (model.quality) {
            uint8_t q = uint8_t(quality[reverse_strand ? _m - 1 - i : i]);
            uint16_t phred = q > '!' ? uint16_t(q - '!') : 0;
            std::fill_n(_sub.begin() + i * _num_ranks, _num_ranks, phred);
        }
    }

    for (const apm_class_t& cls : model.classes) {
        auto [b, e] = range_of(cls.range);
        std::string from = symbols_of(cls.from), to = symbols_of(cls.to);

        for (uint64_t i = b; i < e; i++) {
            if (!contains(from, sym(i))) continue;
            for (char c : to) if (_rank[uint8_t(c)] < _num_ranks) _matches[i * _num_ranks + _rank[uint8_t(c)]] = 1;
        }
    }

    for (const apm_cost_rule_t& rule : model.costs) {
        auto [b, e] = range_of(rule.range);
        std::string from = symbols_of(rule.from), to = symbols_of(rule.to);
        if (!rule.to.empty() && (rule.ops & ~APM_SUB))
            throw std::invalid_argument("text symbols are only allowed for substitutions");
        if (!rule.from.empty() && (rule.ops & (APM_DEL | APM_OPEN)))
            throw std::invalid_argument("pattern symbols are only allowed for substitutions and insertions");

        for (uint64_t i = b; i < e; i++) {
            if (!contains(from, sym(i))) continue;

            if (rule.ops & APM_SUB) {
                for (uint64_t r = 0; r < _num_ranks; r++) if (contains(to, _unrank[r])) _sub[i * _num_ranks + r] = rule.cost;
            }

            if (rule.ops & APM_INS) _ins[i] = rule.cost;
            if (rule.ops & APM_OPEN) _ins_open[i] = rule.cost;
        }

        for (uint64_t x = 1; x < _m; x++) {
            if (!in_range(x, b, e)) continue;
            if (rule.ops & APM_DEL) _del[x] = rule.cost;
            if (rule.ops & APM_OPEN) _del_open[x] = rule.cost;
        }
    }

    // counters: a counter with bound 0 forbids its operations, the others are packed into byte lanes
    for (const apm_counter_t& counter : model.counters) {
        if (counter.max > apm_max_counter_bound) throw std::invalid_argument("counter bound too large");
        auto [b, e] = range_of(counter.range);

        if (counter.max == 0) {
            for (uint64_t i = b; i < e; i++) {
                if (counter.ops & APM_SUB) std::fill_n(_sub.begin() + i * _num_ranks, _num_ranks, apm_forbidden);
                if (counter.ops & APM_INS) _ins[i] = apm_forbidden;
                if (counter.ops & APM_OPEN) _ins_open[i] = apm_forbidden;
            }

            for (uint64_t x = 1; x < _m; x++) {
                if (!in_range(x, b, e)) continue;
                if (counter.ops & APM_DEL) _del[x] = apm_forbidden;
                if (counter.ops & APM_OPEN) _del_open[x] = apm_forbidden;
            }

            continue;
        }

        uint64_t j = _num_counters++;
        uint64_t one = uint64_t(1) << (lane_bits * j);
        _counter_ops[j] = counter.ops;
        _counter_beg[j] = b;
        _counter_end[j] = e;
        _counter_max[j] = counter.max;
        _counter_bias |= (apm_max_counter_bound - counter.max) << (lane_bits * j);

        for (uint64_t i = b; i < e; i++) {
            if (counter.ops & APM_SUB) _sub_counts[i] += one;
            if (counter.ops & APM_INS) _ins_counts[i] += one;
            if (counter.ops & APM_OPEN) _ins_open_counts[i] += one;
        }

        for (uint64_t x = 1; x < _m; x++) {
            if (!in_range(x, b, e)) continue;
            if (counter.ops & APM_DEL) _del_counts[x] += one;
            if (counter.ops & APM_OPEN) _del_open_counts[x] += one;
        }
    }

    // a matched symbol is no substitution; an alignment does not start or end with a deletion; Hamming distance
    // has neither insertions nor deletions
    for (uint64_t x = 0; x < _m * _num_ranks; x++) if (_matches[x]) _sub[x] = 0;
    _del[0] = _del[_m] = apm_forbidden;

    if (metric == HAMMING_DISTANCE) {
        std::fill(_ins.begin(), _ins.end(), apm_forbidden);
        std::fill(_del.begin(), _del.end(), apm_forbidden);
    }

    for (uint64_t i = 0; i < _m; i++) if (_ins[i] == 0) throw std::invalid_argument("an insertion must cost at least 1");
    for (uint64_t x = 1; x < _m; x++) if (_del[x] == 0) throw std::invalid_argument("a deletion must cost at least 1");

    for (uint64_t i = 0; i < _m; i++) _affine |= _ins_open[i] != 0 || _ins_open_counts[i] != 0;
    for (uint64_t x = 0; x <= _m; x++) _affine |= _del_open[x] != 0 || _del_open_counts[x] != 0;

    // a gap of insertions can touch both ends of a part, a gap of deletions one of them
    uint32_t max_ins_open = 0, max_del_open = 0;
    for (uint64_t i = 0; i < _m; i++) max_ins_open = std::max<uint32_t>(max_ins_open, _ins_open[i]);
    for (uint64_t x = 1; x < _m; x++) max_del_open = std::max<uint32_t>(max_del_open, _del_open[x]);
    _prune_slack = uint16_t(std::min<uint32_t>(apm_forbidden, 2 * max_ins_open + max_del_open));

    // with errors of equal cost, every substitution is counted as an error in the search
    const uint16_t error_cost = compute_uniform_cost();

    if (model.tau != 0 || error_cost != 0) {
        set_tau(model.tau != 0 ? model.tau : 1);
        _uniform_cost = error_cost >= _tau ? error_cost : 0;
        choose_cost_dp(model.cost_dp);
        return;
    }

    // tau maximizes the number of symbols without free alternatives that the part without errors of a search with
    // max_errors() + 1 parts has on average; the candidates are 1 and each cost of a substitution plus 1
    static thread_local std::vector<uint8_t> is_candidate;
    is_candidate.assign(uint64_t(_bound) + 2, 0);
    is_candidate[1] = 1;
    if (_bound >= 1) is_candidate[2] = 1;

    for (const apm_cost_rule_t& rule : model.costs) {
        if ((rule.ops & APM_SUB) && rule.cost <= _bound) is_candidate[rule.cost + 1] = 1;
    }

    if (model.quality) {
        for (char q : quality) if (uint8_t(q) > '!' && uint8_t(q) - '!' <= _bound) is_candidate[uint8_t(q) - '!' + 1] = 1;
    }
    std::vector<uint16_t> candidates;
    for (uint64_t c = 1; c < is_candidate.size(); c++) if (is_candidate[c]) candidates.emplace_back(c);
    uint16_t best_tau = 1, applied = 0;
    uint64_t best_exact = 0, best_parts = 1;

    // every tau counts at least the insertions and deletions as errors
    const uint64_t min_parts = max_num_errors(0, _m, true) + 1;

    for (uint16_t tau : candidates) {
        // the positions without free alternatives; a candidate is skipped if it cannot win with the fewest parts
        uint64_t exact = 0;

        for (uint64_t i = 0; i < _m; i++) {
            bool has_free = false;

            for (uint64_t r = 0; r < _num_ranks && !has_free; r++) {
                has_free = _unrank[r] != _pattern[i] && (_matches[i * _num_ranks + r] || _sub[i * _num_ranks + r] < tau);
            }

            exact += !has_free;
        }

        if (exact * best_parts <= best_exact * min_parts) continue;
        set_tau(tau);
        applied = tau;
        uint64_t parts = _max_errors + 1;

        if (exact * best_parts > best_exact * parts) {
            best_tau = tau;
            best_exact = exact;
            best_parts = parts;
        }
    }

    if (applied != best_tau) set_tau(best_tau);
    _uniform_cost = 0;
    choose_cost_dp(model.cost_dp);
}

inline uint16_t apm_costs_t::compute_uniform_cost() const
{
    if (_num_counters > 0 || _affine || _m == 0) return 0;
    const uint16_t cost = _ins[0];
    if (cost == apm_forbidden) return 0;
    for (uint64_t i = 0; i < _m; i++) if (_ins[i] != cost) return 0;
    for (uint64_t x = 1; x < _m; x++) if (_del[x] != cost) return 0;

    for (uint64_t i = 0; i < _m; i++) {
        for (uint64_t r = 0; r < _num_ranks; r++) if (!_matches[i * _num_ranks + r] && _sub[i * _num_ranks + r] != cost) return 0;
    }

    return cost;
}

inline void apm_costs_t::choose_cost_dp(apm_cost_dp_t cost_dp)
{
    // with AUTO, the search prunes without gap openings and counters, or bounds the cost if the costs alone allow
    // at most one insertion or deletion and the search more errors
    const uint64_t indels = max_num_errors(0, _m, true, false);
    const bool auto_on = indels <= 1 && _max_errors > indels;
    const bool auto_prune = !_affine && _num_counters == 0;
    _prunes = can_prune() && (cost_dp == APM_COST_DP_PRUNE || (cost_dp == APM_COST_DP_AUTO && auto_prune));
    _bounds_costs = _prunes || cost_dp == APM_COST_DP_ON || cost_dp == APM_COST_DP_PRUNE
        || (cost_dp == APM_COST_DP_AUTO && auto_on);
}

inline void apm_costs_t::set_tau(uint16_t tau)
{
    _tau = tau;
    _free_syms.clear();
    _free_beg.assign(_m + 1, 0);
    _min_error_cost.assign(_m, apm_forbidden);
    _all_errors_counted.assign(_m, 1);

    for (uint64_t i = 0; i < _m; i++) {
        _free_beg[i] = _free_syms.size();
        _free_syms.emplace_back(_pattern[i]);

        for (uint64_t r = 0; r < _num_ranks; r++) {
            uint16_t c = _sub[i * _num_ranks + r];
            if (_unrank[r] == _pattern[i]) continue;
            if (is_free(i, r)) { _free_syms.emplace_back(_unrank[r]); continue; }
            if (c > _bound) continue;
            _min_error_cost[i] = std::min(_min_error_cost[i], c);
            if (_sub_counts[i] == 0) _all_errors_counted[i] = 0;
        }

        if (_ins[i] <= _bound) {
            _min_error_cost[i] = std::min(_min_error_cost[i], _ins[i]);
            if (_ins_counts[i] == 0) _all_errors_counted[i] = 0;
        }
    }

    _free_beg[_m] = _free_syms.size();
    _max_errors = max_num_errors(0, _m, false);
    _max_indels = max_num_errors(0, _m, true);
}

inline uint64_t apm_costs_t::max_num_errors(uint64_t beg, uint64_t end, bool indels_only, bool counters) const
{
    // the cheapest error at each position
    auto item = [&](uint64_t i) -> uint16_t { return indels_only ? _ins[i] : _min_error_cost[i]; };

    // with the costs alone: the cheapest errors first; a deletion at a boundary can be repeated
    static thread_local std::vector<uint64_t> buckets;
    buckets.assign(_bound + 1, 0);
    for (uint64_t i = beg; i < end; i++) if (item(i) <= _bound) buckets[item(i)]++;
    uint16_t del_min = apm_forbidden;
    for (uint64_t x = beg + 1; x < end; x++) del_min = std::min(del_min, _del[x]);
    uint64_t by_cost = 0, budget = _bound;

    // an alignment of P with an insertion or a deletion opens a gap
    if (indels_only && beg == 0 && end == _m) {
        uint64_t open_min = apm_forbidden;
        for (uint64_t i = 0; i < _m; i++) if (_ins[i] <= _bound) open_min = std::min<uint64_t>(open_min, _ins_open[i]);
        for (uint64_t x = 1; x < _m; x++) if (_del[x] <= _bound) open_min = std::min<uint64_t>(open_min, _del_open[x]);
        if (open_min > budget) return 0;
        budget -= open_min;
    }

    for (uint64_t c = 1; c <= _bound; c++) {
        if (c > budget) break;

        if (del_min <= c) {
            by_cost += budget / del_min;
            budget %= del_min;
            break;
        }

        uint64_t take = std::min<uint64_t>(buckets[c], budget / c);
        by_cost += take;
        budget -= take * c;
    }

    if (!counters) return by_cost;

    // with the counters: every error that increments a counter consumes one unit of a counter that overlaps the range
    uint64_t by_counters = 0;
    bool unbounded = false;

    for (uint64_t i = beg; i < end; i++) {
        if (item(i) > _bound) continue;
        if (indels_only ? _ins_counts[i] == 0 : !_all_errors_counted[i]) by_counters++;
    }

    for (uint64_t x = beg + 1; x < end; x++) unbounded |= _del[x] <= _bound && _del_counts[x] == 0;

    for (uint64_t j = 0; j < _num_counters; j++) {
        uint8_t ops = indels_only ? APM_INS | APM_DEL : APM_SUB | APM_INS | APM_DEL;
        if (_counter_beg[j] < end && beg < _counter_end[j] && (_counter_ops[j] & ops)) by_counters += _counter_max[j];
    }

    return unbounded ? by_cost : std::min(by_cost, by_counters);
}

template <std::ranges::contiguous_range str_t>
std::pair<uint64_t, uint64_t> apm_costs_t::align_without_counters(const str_t& M, cigar_t<uint8_t>* cigar,
    uint64_t* counters, uint64_t reuse) const
{
    const uint64_t kb = _max_indels;
    const uint64_t n = std::min<uint64_t>(std::ranges::size(M), _m + kb);
    const auto* T = std::ranges::data(M);
    const uint64_t states = _affine ? 3 : 1;
    const uint64_t band = 2 * kb + 1;
    const uint64_t q_ins = _affine ? 1 : 0, q_del = _affine ? 2 : 0;
    const uint64_t no_column = std::numeric_limits<uint64_t>::max();
    reuse = std::min(reuse, n);

    // D[(i * band + j + kb - i) * states + q]: the cheapest cost of P[0, i) against M[0, j) ending in state q;
    // rows_done: the rows computed for the last string; first_live[i]: the first column of row i with a cell of at most
    // the bound, or no_column
    static thread_local std::vector<uint16_t> D;
    static thread_local std::vector<uint8_t> ranks;
    static thread_local std::vector<uint64_t> first_live;
    static thread_local uint64_t rows_done = 0;

    if (reuse == 0) {
        D.assign((_m + 1) * band * states, apm_forbidden);
        rows_done = 0;
    }

    first_live.resize(_m + 1);
    ranks.resize(n);
    for (uint64_t j = reuse; j < n; j++) ranks[j] = _rank[std::bit_cast<uint8_t>(T[j])];
    auto at = [&](uint64_t i, uint64_t j, uint64_t q) -> uint16_t& { return D[(i * band + (j + kb - i)) * states + q]; };
    auto in_band = [&](uint64_t i, uint64_t j) { return j + kb >= i && j <= i + kb && j <= n; };

    // the cost of a predecessor cell plus an operation, or apm_forbidden
    auto plus = [&](uint64_t v, uint64_t c) -> uint64_t { return v == apm_forbidden || c == apm_forbidden ? apm_forbidden : v + c; };

    // the sum of a cell and the cost of an operation, or apm_forbidden if it exceeds the bound
    auto sum = [&](uint32_t v, uint32_t c) -> uint16_t { return apm_bounded_sum(v, c, _bound); };

    D[kb * states] = 0;

    for (uint64_t i = 0; i <= _m; i++) {
        // the cell (i, j) has the index x = j + kb - i in row i; (i - 1, j - 1) has x, (i - 1, j) has x + 1 in row i - 1
        uint16_t* cur = D.data() + i * band * states;
        const uint16_t* prev = i > 0 ? cur - band * states : cur;
        const uint16_t* sub_row = i > 0 ? _sub.data() + (i - 1) * _num_ranks : nullptr;
        const uint32_t ins = i > 0 ? _ins[i - 1] : apm_forbidden, ins_open = i > 0 ? _ins_open[i - 1] : 0;
        const uint32_t del = _del[i], del_open = _del_open[i];
        uint64_t j_beg = i > kb ? i - kb : 0;
        uint64_t j_end = std::min(n, i + kb);
        uint64_t j_from = std::max<uint64_t>(j_beg, i == 0 ? 1 : 0);
        uint64_t live = i == 0 ? 0 : no_column;

        // the cells in the columns [j_from, reuse] are those of the last string
        if (i < rows_done && j_from <= reuse) {
            if (first_live[i] <= reuse) live = std::min(live, first_live[i]);
            j_from = reuse + 1;
        }

        for (uint64_t j = j_from; j <= j_end; j++) {
            const uint64_t x = j + kb - i;
            const uint32_t sub = i > 0 && j > 0 && ranks[j - 1] < _num_ranks ? sub_row[ranks[j - 1]] : apm_forbidden;
            const bool from_above = i > 0 && x + 1 < band, from_left = j > 0 && x > 0;

            if (states == 1) {
                uint16_t v = sum(prev[x], sub);
                if (from_above) v = std::min(v, sum(prev[x + 1], ins));
                if (from_left) v = std::min(v, sum(cur[x - 1], del));
                cur[x] = v;
                if (v != apm_forbidden) live = std::min(live, j);
                continue;
            }

            uint16_t* c3 = cur + 3 * x;
            const uint16_t* d3 = prev + 3 * x;
            c3[0] = std::min({sum(d3[0], sub), sum(d3[1], sub), sum(d3[2], sub)});

            if (from_above) {
                const uint16_t* a3 = prev + 3 * (x + 1);
                c3[1] = std::min({sum(a3[0], ins + ins_open), sum(a3[1], ins), sum(a3[2], ins + ins_open)});
            }

            if (from_left) {
                const uint16_t* l3 = cur + 3 * (x - 1);
                c3[2] = std::min({sum(l3[0], del + del_open), sum(l3[1], del + del_open), sum(l3[2], del)});
            }

            if (std::min({c3[0], c3[1], c3[2]}) != apm_forbidden) live = std::min(live, j);
        }

        first_live[i] = live;

        if (live == no_column) {
            rows_done = i + 1;
            return {uint64_t(_bound) + 1, 0};
        }
    }

    rows_done = _m + 1;

    // the cheapest cell in row m, the shortest window among them
    uint64_t best_cost = uint64_t(_bound) + 1, best_len = 0, best_q = 0;

    for (uint64_t j = std::max<uint64_t>(1, _m > kb ? _m - kb : 0); j <= std::min(n, _m + kb); j++) {
        for (uint64_t q = 0; q < states; q++) {
            if (at(_m, j, q) < best_cost) { best_cost = at(_m, j, q); best_len = j; best_q = q; }
        }
    }

    if ((cigar == nullptr && counters == nullptr) || best_cost > _bound) return {best_cost, best_len};

    // trace the alignment back from the cheapest cell, collecting its operations and its counters
    static thread_local std::vector<cigar_run_t<uint8_t>> runs;
    runs.clear();
    uint64_t path_counters = 0;
    uint64_t i = _m, j = best_len, q = best_q;

    while (i > 0 || j > 0) {
        uint64_t v = at(i, j, q);
        bool found = false;

        if (q == 0 && i > 0 && j > 0 && in_band(i - 1, j - 1)) {
            uint8_t t = std::bit_cast<uint8_t>(T[j - 1]);
            uint64_t r = _rank[t];
            uint64_t c = matches(i - 1, r) ? 0 : sub(i - 1, r);

            for (uint64_t q0 = 0; q0 < states && !found; q0++) {
                if (plus(at(i - 1, j - 1, q0), c) != v) continue;
                runs.emplace_back(t == _pattern[i - 1] ? MATCH<uint8_t>() : MISMATCH(t));
                if (!matches(i - 1, r)) path_counters += _sub_counts[i - 1];
                i--; j--; q = q0; found = true;
            }
        }

        if (!found && q == q_ins && i > 0 && in_band(i - 1, j)) {
            for (uint64_t q0 = 0; q0 < states && !found; q0++) {
                if (plus(at(i - 1, j, q0), plus(_ins[i - 1], q0 != q_ins ? _ins_open[i - 1] : 0)) != v) continue;
                runs.emplace_back(INS(_pattern[i - 1]));
                path_counters += _ins_counts[i - 1] + (q0 != q_ins ? _ins_open_counts[i - 1] : 0);
                i--; q = q0; found = true;
            }
        }

        if (!found && q == q_del && j > 0 && in_band(i, j - 1)) {
            for (uint64_t q0 = 0; q0 < states && !found; q0++) {
                if (plus(at(i, j - 1, q0), plus(_del[i], q0 != q_del ? _del_open[i] : 0)) != v) continue;
                runs.emplace_back(DEL(std::bit_cast<uint8_t>(T[j - 1])));
                path_counters += _del_counts[i] + (q0 != q_del ? _del_open_counts[i] : 0);
                j--; q = q0; found = true;
            }
        }

        assert(found);
        if (!found) break;
    }

    if (counters != nullptr) *counters = path_counters;

    if (cigar != nullptr) {
        cigar->clear();
        for (auto it = runs.rbegin(); it != runs.rend(); ++it) cigar_push(*cigar, *it);
    }

    return {best_cost, best_len};
}

template <std::ranges::contiguous_range str_t>
void apm_costs_t::align_prefixes(const str_t& strings, const std::vector<uint64_t>& beg,
    std::vector<std::pair<uint64_t, uint64_t>>& out) const
{
    const auto* S = std::ranges::data(strings);
    const uint64_t num = beg.size() - 1;
    auto string = [&](uint64_t x) { return std::span(S + beg[x], beg[x + 1] - beg[x]); };
    std::vector<uint32_t> order(num);
    std::iota(order.begin(), order.end(), 0);

    std::sort(order.begin(), order.end(), [&](uint32_t x, uint32_t y) {
        auto a = string(x), b = string(y);
        return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
    });

    out.resize(num);

    for (uint64_t y = 0; y < num; y++) {
        auto M = string(order[y]);
        uint64_t lcp = 0;

        if (y > 0) {
            auto L = string(order[y - 1]);
            lcp = std::mismatch(M.begin(), M.end(), L.begin(), L.end()).first - M.begin();
        }

        out[order[y]] = align_reusing(M, nullptr, lcp);
    }
}

template <std::ranges::contiguous_range str_t>
std::pair<uint64_t, uint64_t> apm_costs_t::align_reusing(const str_t& M, cigar_t<uint8_t>* cigar, uint64_t reuse) const
{
    if (_num_counters == 0) return align_without_counters(M, cigar, nullptr, reuse);

    // without the counters, the cost is at most the cost with them; if the cheapest alignment without them keeps them,
    // it is the cheapest alignment with them
    uint64_t counters;
    auto [cost, len] = align_without_counters(M, cigar, &counters, reuse);
    if (cost > _bound || !exceeds(counters)) return {cost, len};

    enum : uint8_t { step_start, step_match, step_sub, step_ins, step_del };
    const uint64_t n = std::ranges::size(M);
    const auto* T = std::ranges::data(M);
    const uint64_t kb = _max_indels;
    const uint64_t states = _affine ? 3 : 1;
    const uint64_t band = 2 * kb + 1;

    static thread_local std::vector<label_t> labels;
    static thread_local std::vector<std::pair<uint32_t, uint32_t>> cells;
    static thread_local std::vector<label_t> tmp;
    labels.clear();
    cells.assign((_m + 1) * band * states, {0, 0});

    // cell (i, j, q) holds the labels of P[0, i) against M[0, j) ending in state q; j lies in [i - kb, i + kb]
    auto cell = [&](uint64_t i, uint64_t j, uint64_t q) -> std::pair<uint32_t, uint32_t>& {
        return cells[(i * band + (j + kb - i)) * states + q];
    };

    auto in_band = [&](uint64_t i, uint64_t j) { return j + kb >= i && j <= i + kb && j <= n; };

    auto dominates = [](const label_t& a, const label_t& b) {
        return a.cost <= b.cost && (((b.counters | lane_high_bits) - a.counters) & lane_high_bits) == lane_high_bits;
    };

    auto offer = [&](const label_t& l) {
        for (const label_t& t : tmp) if (dominates(t, l)) return;
        std::erase_if(tmp, [&](const label_t& t) { return dominates(l, t); });
        tmp.emplace_back(l);
    };

    // extends the labels of cell (i0, j0) (all states) by an operation into state q
    auto extend = [&](uint64_t i0, uint64_t j0, uint8_t step, uint8_t q, uint16_t cost, uint64_t step_counters,
        uint16_t open_cost, uint64_t open_counters
    ) {
        if (cost == apm_forbidden) return;

        for (uint64_t q0 = 0; q0 < states; q0++) {
            auto [beg, end] = cell(i0, j0, q0);

            for (uint32_t x = beg; x < end; x++) {
                const label_t& l = labels[x];
                bool opens = (step == step_ins || step == step_del) && l.state != q;
                uint64_t c = uint64_t(l.cost) + cost + (opens ? open_cost : 0);
                if (c > _bound) continue;
                uint64_t k = l.counters + step_counters + (opens ? open_counters : 0);
                if (exceeds(k)) continue;
                offer(label_t{uint16_t(c), step, q, x, k});
            }
        }
    };

    labels.emplace_back(label_t{0, step_start, 0, no_label, 0});
    cell(0, 0, 0) = {0, 1};

    for (uint64_t i = 0; i <= _m; i++) {
        uint64_t j_beg = i > kb ? i - kb : 0;
        uint64_t j_end = std::min(n, i + kb);
        uint64_t row_beg = labels.size();

        for (uint64_t j = j_beg; j <= j_end; j++) {
            if (i == 0 && j == 0) continue;

            for (uint64_t q = 0; q < states; q++) {
                tmp.clear();
                uint8_t q_ins = _affine ? 1 : 0, q_del = _affine ? 2 : 0;

                if (q == 0 && i > 0 && j > 0 && in_band(i - 1, j - 1)) {
                    uint64_t r = _rank[std::bit_cast<uint8_t>(T[j - 1])];
                    bool match = matches(i - 1, r);
                    extend(i - 1, j - 1, match ? step_match : step_sub, 0, match ? 0 : sub(i - 1, r),
                        match ? 0 : _sub_counts[i - 1], 0, 0);
                }

                if (q == q_ins && i > 0 && in_band(i - 1, j)) {
                    extend(i - 1, j, step_ins, q_ins, _ins[i - 1], _ins_counts[i - 1], _ins_open[i - 1], _ins_open_counts[i - 1]);
                }

                if (q == q_del && j > 0 && in_band(i, j - 1)) {
                    extend(i, j - 1, step_del, q_del, _del[i], _del_counts[i], _del_open[i], _del_open_counts[i]);
                }

                uint32_t beg = labels.size();
                labels.insert(labels.end(), tmp.begin(), tmp.end());
                cell(i, j, q) = {beg, uint32_t(labels.size())};
            }
        }

        if (i > 0 && labels.size() == row_beg) return {uint64_t(_bound) + 1, 0};
    }

    // the cheapest label in row m, the shortest window among them
    uint64_t best_cost = uint64_t(_bound) + 1, best_len = 0;
    uint32_t best = no_label;

    for (uint64_t j = std::max<uint64_t>(1, _m > kb ? _m - kb : 0); j <= std::min(n, _m + kb); j++) {
        for (uint64_t q = 0; q < states; q++) {
            auto [beg, end] = cell(_m, j, q);

            for (uint32_t x = beg; x < end; x++) {
                if (labels[x].cost < best_cost) { best_cost = labels[x].cost; best_len = j; best = x; }
            }
        }
    }

    if (cigar != nullptr && best != no_label) {
        cigar->clear();
        std::vector<cigar_run_t<uint8_t>> runs;
        uint64_t i = _m, j = best_len;

        for (uint32_t x = best; labels[x].step != step_start; x = labels[x].prev) {
            switch (labels[x].step) {
                case step_match:
                    if (std::bit_cast<uint8_t>(T[j - 1]) == _pattern[i - 1]) runs.emplace_back(MATCH<uint8_t>());
                    else runs.emplace_back(MISMATCH(std::bit_cast<uint8_t>(T[j - 1])));
                    i--; j--;
                    break;
                case step_sub: runs.emplace_back(MISMATCH(std::bit_cast<uint8_t>(T[j - 1]))); i--; j--; break;
                case step_ins: runs.emplace_back(INS(_pattern[i - 1])); i--; break;
                case step_del: runs.emplace_back(DEL(std::bit_cast<uint8_t>(T[j - 1]))); j--; break;
            }
        }

        for (auto it = runs.rbegin(); it != runs.rend(); ++it) cigar_push(*cigar, *it);
    }

    return {best_cost, best_len};
}

inline void apm_cost_band::init(const apm_costs_t& costs, uint64_t beg, uint64_t end, direction_t dir, uint16_t cost_in,
    uint64_t counters_in, uint64_t num_rows, uint64_t restart_rows)
{
    _costs = &costs;
    _bound = costs.bound();
    _beg = beg;
    _len = end - beg;
    _reads_right = dir == RIGHT;
    _cols_after = costs.max_indels();
    _cols_before = _cols_after + restart_rows;
    _band = _cols_before + _cols_after + 1;
    _states = costs.is_affine() ? 3 : 1;
    _rows.assign(num_rows * _band * _states, apm_forbidden);
    _rows[_cols_before * _states] = cost_in;

    // the costs of the operations on the part; rank costs.num_ranks() is a symbol that does not occur in the index
    _num_ranks = costs.num_ranks() + 1;
    _sub_vals.assign(_len * _num_ranks, apm_forbidden);
    _ins_vals.resize(_len);
    _del_vals.resize(_len + 1);
    _ins_open_in.assign(_len, 0);
    _ins_open_out.assign(_len + 1, 0);
    _del_open.assign(_len + 1, 0);

    for (uint64_t c = 0; c < _len; c++) {
        uint64_t p = pos(c);
        _ins_vals[c] = costs.ins(p);
        for (uint64_t r = 0; r < costs.num_ranks(); r++) _sub_vals[c * _num_ranks + r] = costs.matches(p, r) ? 0 : costs.sub(p, r);
    }

    for (uint64_t c = 0; c <= _len; c++) _del_vals[c] = costs.del(boundary(c));

    // the costs of opening the gaps that do not touch an end of the part
    if (_states == 3) {
        for (uint64_t c = 1; c < _len; c++) {
            if (_reads_right) _ins_open_in[c] = costs.ins_open(pos(c));
            else _ins_open_out[c] = costs.ins_open(pos(c - 1));
            _del_open[c] = costs.del_open(boundary(c));
        }
    }

    // the counters of the operations on the part
    _counts = costs.num_counters() > 0 && _states == 1;

    if (_counts) {
        _counters.assign(num_rows * _band, 0);
        _counters[_cols_before] = counters_in;
        _sub_count_vals.assign(_len * _num_ranks, 0);
        _ins_count_vals.resize(_len);
        _del_count_vals.resize(_len + 1);

        for (uint64_t c = 0; c < _len; c++) {
            uint64_t p = pos(c);
            _ins_count_vals[c] = costs.ins_count(p);
            for (uint64_t r = 0; r < costs.num_ranks(); r++) _sub_count_vals[c * _num_ranks + r] = costs.matches(p, r) ? 0 : costs.sub_count(p);
        }

        for (uint64_t c = 0; c <= _len; c++) _del_count_vals[c] = costs.del_count(boundary(c));
    }

    // row 0: the first symbols of the part are inserted
    const uint64_t q_ins = _states == 3 ? 1 : 0;

    for (uint64_t c = 1; c <= std::min(_len, _cols_after); c++) {
        uint16_t prev = std::min(_rows[(c - 1 + _cols_before) * _states], _rows[(c - 1 + _cols_before) * _states + q_ins]);
        _rows[(c + _cols_before) * _states + q_ins] = add_bounded(prev, _ins_vals[c - 1]);
        if (_counts) _counters[c + _cols_before] = _counters[c - 1 + _cols_before] + _ins_count_vals[c - 1];
    }
}

inline uint64_t apm_cost_band::final_counters(uint64_t i) const
{
    if (!_counts || _len + _cols_before < i || _len > i + _cols_after) return 0;
    return _counters[i * _band + _len + _cols_before - i];
}

inline uint16_t apm_cost_band::final_cost(uint64_t i) const
{
    if (_len + _cols_before < i || _len > i + _cols_after) return apm_forbidden;
    const uint16_t* cell = _rows.data() + (i * _band + _len + _cols_before - i) * _states;
    return _states == 1 ? cell[0] : std::min({cell[0], cell[1], cell[2]});
}

inline uint16_t apm_cost_band::min_partial_cost(uint64_t i) const
{
    const uint16_t* row = _rows.data() + i * _band * _states;
    uint16_t cost = apm_forbidden;

    for (uint64_t c = i > _cols_before ? i - _cols_before : 0; c <= std::min(_len, i + _cols_after) && c < _len; c++) {
        for (uint64_t q = 0; q < _states; q++) cost = std::min(cost, row[(c + _cols_before - i) * _states + q]);
    }

    return cost;
}

inline bool apm_cost_band::compute_row(uint64_t i, uint8_t sym)
{
    const uint64_t r = _costs->rank(sym);
    uint16_t* row = _rows.data() + i * _band * _states;
    const uint16_t* prev = row - _band * _states;
    bool alive = false;

    for (uint64_t c = i > _cols_before ? i - _cols_before : 0; c <= std::min(_len, i + _cols_after); c++) {
        // the cell (i, c) has index x in row i; (i - 1, c - 1) has x and (i - 1, c) has x + 1 in row i - 1, (i, c - 1)
        // has x - 1 in row i
        const uint64_t x = c + _cols_before - i;
        const bool from_above = x + 1 < _band;
        const uint16_t sub = c > 0 ? _sub_vals[(c - 1) * _num_ranks + r] : apm_forbidden;

        if (_states == 1 && _counts) {
            // the counters of the cheapest predecessor (the first one of equal cost)
            uint64_t* cnt = _counters.data() + i * _band;
            const uint64_t* cnt_prev = cnt - _band;
            uint16_t v = from_above ? add_bounded(prev[x + 1], _del_vals[c]) : apm_forbidden;
            uint64_t k = from_above ? cnt_prev[x + 1] + _del_count_vals[c] : 0;

            if (c > 0) {
                const uint16_t v_sub = add_bounded(prev[x], sub);
                if (v_sub < v) { v = v_sub; k = cnt_prev[x] + _sub_count_vals[(c - 1) * _num_ranks + r]; }
                const uint16_t v_ins = x > 0 ? add_bounded(row[x - 1], _ins_vals[c - 1]) : apm_forbidden;
                if (v_ins < v) { v = v_ins; k = cnt[x - 1] + _ins_count_vals[c - 1]; }
            }

            row[x] = v;
            cnt[x] = k;
            alive |= v != apm_forbidden;
            continue;
        }

        if (_states == 1) {
            uint16_t v = from_above ? add_bounded(prev[x + 1], _del_vals[c]) : apm_forbidden;

            if (c > 0) {
                v = std::min(v, add_bounded(prev[x], sub));
                if (x > 0) v = std::min(v, add_bounded(row[x - 1], _ins_vals[c - 1]));
            }

            row[x] = v;
            alive |= v != apm_forbidden;
            continue;
        }

        uint16_t* cur = row + 3 * x;
        cur[0] = cur[1] = cur[2] = apm_forbidden;

        if (c > 0) {
            const uint16_t* d = prev + 3 * x;
            const uint16_t out = _ins_open_out[c - 1];
            cur[0] = std::min({add_bounded(d[0], sub), add_bounded(d[1], sub, out), add_bounded(d[2], sub)});

            if (x > 0) {
                const uint16_t* l = row + 3 * (x - 1);
                const uint16_t ins = _ins_vals[c - 1], open = _ins_open_in[c - 1];
                cur[1] = std::min({add_bounded(l[0], ins, open), add_bounded(l[1], ins), add_bounded(l[2], ins, open)});
            }
        }

        if (from_above) {
            const uint16_t* a = prev + 3 * (x + 1);
            const uint16_t del = _del_vals[c], open = _del_open[c];
            cur[2] = std::min({add_bounded(a[0], del, open), add_bounded(a[1], del, open + _ins_open_out[c]), add_bounded(a[2], del)});
        }

        alive |= std::min({cur[0], cur[1], cur[2]}) != apm_forbidden;
    }

    return alive;
}

inline bool apm_cost_band::restart(uint64_t i, uint16_t cost, uint64_t counters)
{
    assert(0 < i && i + _cols_after < _band);
    uint16_t* row = _rows.data() + i * _band * _states;
    uint64_t* cnt_row = _counts ? _counters.data() + i * _band : nullptr;
    uint16_t& start = row[(_cols_before - i) * _states];

    if (cost < start) {
        start = cost;
        if (_counts) cnt_row[_cols_before - i] = counters;
    }

    // the insertions of the symbols of the part after the new start (no gap is opened at the start of the part)
    for (uint64_t c = 1; c <= std::min(_len, i + _cols_after); c++) {
        const uint64_t x = c + _cols_before - i;
        const uint16_t ins = _ins_vals[c - 1];

        if (_states == 1) {
            const uint16_t v_ins = add_bounded(row[x - 1], ins);

            if (v_ins < row[x]) {
                row[x] = v_ins;
                if (_counts) cnt_row[x] = cnt_row[x - 1] + _ins_count_vals[c - 1];
            }

            continue;
        }

        const uint16_t* l = row + 3 * (x - 1);
        const uint16_t open = _ins_open_in[c - 1];
        row[3 * x + 1] = std::min({row[3 * x + 1], add_bounded(l[0], ins, open), add_bounded(l[1], ins), add_bounded(l[2], ins, open)});
    }

    for (uint64_t x = (_cols_before - i) * _states; x < (std::min(_len, i + _cols_after) + _cols_before - i + 1) * _states; x++) {
        if (row[x] != apm_forbidden) return true;
    }

    return false;
}

#endif // MOVE_R_ALGORITHMS_MOVE_RB_APM_MODEL_TPP
