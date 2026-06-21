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
#include <filesystem>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <misc/utils.hpp>
#include <misc/log.hpp>
#include <misc/files.hpp>

/** @brief a command-line option of the form "name value" */
struct Option {
    std::string name;
    std::string value;
};

/** @brief the result of parsing the command-line arguments */
struct CommandLineArguments {
    std::string command;
    std::vector<std::string> last_param;
    std::vector<Option> value_options;
    std::set<std::string> literal_options;
    bool success;
};

/** @brief the CommandLineArguments value returned when parsing fails */
const CommandLineArguments FAILURE = {
    std::string(),
    std::vector<std::string>(),
    std::vector<Option>(),
    std::set<std::string>(),
    false
};

/**
 * @brief parses the command-line arguments into options and fixed parameters
 * @param argc the number of command-line arguments
 * @param argv the command-line arguments
 * @param allowed_value_options the set of allowed options that take a value
 * @param allowed_literal_options the set of allowed literal (flag) options
 * @param fixed_parameter_count the number of fixed (positional) parameters expected at the end
 * @return the parsed arguments (with success == false if parsing failed)
 */
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

/**
 * @brief prints an error message about a malformed patterns-file header and exits
 */
[[noreturn]] void cli_header_error()
{
    std::cout << "Error: malformed header in patterns file." << std::endl;
    std::cout << "Take a look here for more info on the file format: "
              << "http://pizzachili.dcc.uchile.cl/experiments.html" << std::endl;
    exit(1);
}

/**
 * @brief returns the value following a key in a header line (e.g. "length=" in "number=100 length=10")
 * @param header the header line
 * @param key the key to look for
 * @return the value following the key, or an empty string if the key is not present
 */
std::string get_header_value(std::string header, std::string key)
{
    uint64_t key_pos = header.find(key);
    if (key_pos == std::string::npos) return "";
    uint64_t start_pos = key_pos + key.length();
    if (start_pos >= header.size()) return "";
    uint64_t end_pos = header.substr(start_pos, std::string::npos).find(" ");
    return header.substr(start_pos).substr(0, end_pos);
}

/**
 * @brief parses a strictly-positive header value, exiting with an error if it is missing or malformed
 * @param header the header line
 * @param key the key to look for
 * @return the parsed value (guaranteed to be >= 1)
 */
uint64_t get_positive_header_value(std::string header, std::string key)
{
    std::string value = get_header_value(header, key);
    if (value.empty()) cli_header_error();

    try {
        int64_t parsed = std::stoll(value);
        if (parsed < 1) cli_header_error();
        return parsed;
    } catch (const std::exception&) {
        cli_header_error();
    }
}

/**
 * @brief reads the pattern length from a patterns-file header
 * @param header the header line
 * @return the pattern length (>= 1; exits with an error if the header is malformed)
 */
uint64_t get_pattern_length(std::string header)
{
    return get_positive_header_value(header, "length=");
}

/**
 * @brief reads the number of patterns from a patterns-file header
 * @param header the header line
 * @return the number of patterns (>= 1; exits with an error if the header is malformed)
 */
uint64_t get_pattern_count(std::string header)
{
    return get_positive_header_value(header, "number=");
}

/**
 * @brief checks that path refers to a readable regular file, printing an error and exiting otherwise
 * @param path the path to check
 */
void require_file(const std::string& path)
{
    if (!std::filesystem::exists(path) || !std::filesystem::is_regular_file(path)) {
        std::cout << "Error: file '" << path << "' does not exist." << std::endl;
        exit(1);
    }

    std::ifstream f(path);
    if (!f.good()) {
        std::cout << "Error: could not open file '" << path << "'." << std::endl;
        exit(1);
    }
}

/**
 * @brief parses the value of a numeric command-line option, exiting with an error if it is not an integer
 * @param value the option value
 * @param option_name the option name (used in the error message)
 * @return the parsed integer
 */
int64_t parse_int_arg(const std::string& value, const std::string& option_name)
{
    try {
        size_t consumed = 0;
        int64_t parsed = std::stoll(value, &consumed);
        if (consumed != value.size()) throw std::invalid_argument("trailing characters");
        return parsed;
    } catch (const std::exception&) {
        std::cout << "Error: the " << option_name << " option expects an integer (got '"
                  << value << "')." << std::endl;
        exit(1);
    }
}

/**
 * @brief loads pattern_count patterns of length pattern_length from the input stream
 * @param pattern_in the input stream (positioned after the header)
 * @param pattern_length the length of each pattern
 * @param pattern_count the number of patterns to read
 * @return the loaded patterns
 */
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

/**
 * @brief converts a vector of strings to a vector of byte vectors
 * @param v a vector of strings
 * @return a vector containing each string as a vector of bytes
 */
std::vector<std::vector<uint8_t>> strToUint8Vec(std::vector<std::string>& v)
{
    std::vector<std::vector<uint8_t>> result;
    result.reserve(v.size());
    for (std::string s : v) {
        result.push_back(std::vector<uint8_t>(s.begin(), s.end()));
    }

    return result;
}