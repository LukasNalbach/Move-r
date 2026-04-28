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

#pragma once

#include <cstdint>
#include <iostream>
#include <fstream>
#include <set>
#include <vector>

struct Option {
    std::string name;
    std::string value;
};

struct CommandLineArguments {
    std::string command;
    std::vector<std::string> last_param;
    std::vector<Option> value_options;
    std::set<std::string> literal_options;
    bool success;
};

const CommandLineArguments FAILURE = {
    std::string(),
    std::vector<std::string>(),
    std::vector<Option>(),
    std::set<std::string>(),
    false
};

CommandLineArguments parse_args(
    int argc, char** argv,
    std::set<std::string>& allowed_value_options,
    std::set<std::string>& allowed_literal_options,
    int fixed_parameter_count
) {
    if (argc < fixed_parameter_count + 1) {
        return FAILURE;
    }

    std::string command = argv[0];
    std::vector<Option> value_options;
    std::set<std::string> literal_options;
    std::set<std::string> found;
    int i = 1;

    while (i < argc - fixed_parameter_count) { // construct_lzend_of_reverse options
        std::string option = argv[i++];

        if (found.count(option) == 1) {
            std::cout << "Error: duplicate occurence of option " << option << "." << std::endl;
            return FAILURE;
        }

        if (allowed_value_options.count(option) == 1) {
            if (i > argc - 1) {
                std::cout << "Error: missing parameter after the " << option << " option." << std::endl;
                return FAILURE;
            }

            std::string value = argv[i++];
            value_options.push_back({ option, value });
            found.insert(option);
        } else if (allowed_literal_options.count(option) == 1) {
            literal_options.insert(option);
            found.insert(option);
        } else { // unknown option provided
            std::cout << "Error: unknown option " << option << "." << std::endl;
            return FAILURE;
        }
    }

    if (i >= argc) {
        std::cout << "Error: missing parameter" << std::endl;
        return FAILURE;
    }

    std::vector<std::string> last_param;
    last_param.reserve(fixed_parameter_count);

    while (i < argc) {
        last_param.push_back(argv[i++]);
    }

    return { command, last_param, value_options, literal_options, true };
}

std::string get_header_value(std::string header, std::string key)
{
    uint64_t start_pos = header.find(key) + key.length();
    uint64_t end_pos = header.substr(start_pos, std::string::npos).find(" ");
    return header.substr(start_pos).substr(0, end_pos);
}

uint64_t get_pattern_length(std::string header)
{
    std::string len = get_header_value(header, "length=");
    return std::stoll(len);
}

uint64_t get_pattern_count(std::string header)
{
    std::string count = get_header_value(header, "number=");
    return std::stoll(count);
}

std::vector<std::string> load_patterns(std::ifstream& pattern_in, uint64_t pattern_length, uint64_t pattern_count)
{
    std::vector<std::string> patterns;

    for (uint64_t i = 0; i < pattern_count; ++i) {
        std::string pattern;

        for (uint64_t j = 0; j < pattern_length; ++j) {
            char c;
            pattern_in.get(c);
            pattern += c;
        }

        patterns.push_back(pattern);
    }

    return patterns;
}

std::vector<std::vector<uint8_t>> strToUint8Vec(std::vector<std::string>& v)
{
    std::vector<std::vector<uint8_t>> result;
    result.reserve(v.size());
    for (std::string s : v) {
        result.push_back(std::vector<uint8_t>(s.begin(), s.end()));
    }

    return result;
}