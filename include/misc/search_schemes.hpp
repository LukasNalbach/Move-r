#pragma once

#include <cstdint>
#include <vector>

enum distance_metric_t : int8_t {
    NO_METRIC = -1,
    HAMMING_DISTANCE = 0,
    EDIT_DISTANCE = 1
};

struct search_step_t {
    uint8_t part;
    uint8_t k_min;
    uint8_t k_max;
};

using search_t = std::vector<search_step_t>;

struct search_scheme_t {
    distance_metric_t dist_metr = NO_METRIC;
    uint8_t k_max = 0;
    uint8_t parts = 1;
    std::vector<search_t> searches = {{{0, 0, 0}}};
};

static search_scheme_t pigeon_hole_scheme(uint8_t k_max, distance_metric_t dist_metr)
{
    uint8_t parts = k_max + 1;
    std::vector<search_t> searches;
    searches.reserve(k_max + 1);

    for (int16_t i = 0; i < k_max + 1; i++) {
        search_t search;
        search.reserve(parts);
        search.emplace_back(search_step_t{.part = i, .k_min = 0, .k_max = 0});

        for (int16_t j = i + 1; j < parts; j++) {
            search.emplace_back(search_step_t{.part = j, .k_min = 0, .k_max = k_max});
        }

        for (int16_t j = i - 1; j >= 0; j--) {
            search.emplace_back(search_step_t{.part = j, .k_min = 0, .k_max = k_max});
        }

        searches.emplace_back(std::move(search));
    }

    return search_scheme_t {
        .dist_metr = dist_metr,
        .k_max = k_max,
        .parts = parts,
        .searches = std::move(searches)
    };
}

static search_scheme_t suffix_filter_scheme(uint8_t k_max, distance_metric_t dist_metr)
{
    uint8_t parts = k_max + 1;
    std::vector<search_t> searches;
    searches.reserve(k_max + 1);

    for (int16_t i = 0; i < k_max + 1; i++) {
        search_t search;
        search.reserve(parts);

        for (int16_t j = i; j < parts; j++) {
            search.emplace_back(search_step_t{.part = j, .k_min = 0, .k_max = j - i});
        }

        for (int16_t j = i - 1; j >= 0; j--) {
            search.emplace_back(search_step_t{.part = j, .k_min = 0, .k_max = k_max});
        }

        searches.emplace_back(std::move(search));
    }

    return search_scheme_t {
        .dist_metr = dist_metr,
        .k_max = k_max,
        .parts = parts,
        .searches = std::move(searches)
    };
}

static search_scheme_t zero_one_scheme(uint8_t k_max, distance_metric_t dist_metr)
{
    uint8_t parts = k_max + 2;
    std::vector<search_t> searches;
    searches.reserve(k_max + 1);

    for (int16_t i = 0; i < k_max + 1; i++) {
        search_t search;
        search.reserve(parts);
        search.emplace_back(search_step_t{.part = i, .k_min = 0, .k_max = 0});
        search.emplace_back(search_step_t{.part = i + 1, .k_min = 0, .k_max = i == k_max ? 0 : 1});

        for (int16_t j = i + 2; j < parts; j++) {
            search.emplace_back(search_step_t{.part = j, .k_min = 0, .k_max = k_max});
        }

        for (int16_t j = i - 1; j >= 0; j--) {
            search.emplace_back(search_step_t{.part = j, .k_min = 0, .k_max = k_max});
        }

        searches.emplace_back(std::move(search));
    }

    return search_scheme_t {
        .dist_metr = dist_metr,
        .k_max = k_max,
        .parts = parts,
        .searches = std::move(searches)
    };
}

static std::vector<uint8_t> parse_bracket(std::string content, int64_t parts)
{
    std::replace(content.begin(), content.end(), ',', ' ');
    
    std::stringstream ss(content);
    std::vector<uint8_t> numbers;
    int64_t num;
    
    while (ss >> num) {
        if (!(0 <= num && num < parts)) {
            print_search_scheme_error();
        }

        numbers.push_back(num);
    }

    return numbers;
}

static search_scheme_t parse_search_scheme(const std::string& str, distance_metric_t dist_metr)
{
    std::stringstream ss_input(str);
    std::string line;
    std::getline(ss_input, line);

    int64_t parts = value_from_key(line, "p=");
    int64_t k_max = value_from_key(line, "k=");

    if (parts == -1 || k_max == -1) {
        print_search_scheme_error();
    }

    std::vector<search_t> searches;

    while (std::getline(ss_input, line)) {
        if (line.empty()) {
            print_search_scheme_error();
        }

        std::vector<std::vector<uint8_t>> bracket_contents;
        uint64_t pos = 0;

        while ((pos = line.find('{', pos)) != std::string::npos) {
            uint64_t end = line.find('}', pos);
            if (end == std::string::npos) break;

            std::string content = line.substr(pos + 1, end - pos - 1);
            bracket_contents.emplace_back(parse_bracket(content, parts));
            
            pos = end + 1;
        }

        if (bracket_contents.size() != 3 ||
            bracket_contents[0].size() != parts ||
            bracket_contents[1].size() != parts ||
            bracket_contents[2].size() != parts
        ) {
            print_search_scheme_error();
        }

        search_t search;
        search.reserve(parts);

        for (int64_t i = 0; i < parts; i++) {
            search.emplace_back(search_step_t{
                .part =  (uint8_t) bracket_contents[0][i],
                .k_min = (uint8_t) bracket_contents[1][i],
                .k_max = (uint8_t) bracket_contents[2][i]
            });
        }

        searches.emplace_back(search);
    }

    return search_scheme_t {
        .dist_metr = dist_metr,
        .k_max = (uint8_t) k_max,
        .parts = (uint8_t) parts,
        .searches = std::move(searches)
    };
}