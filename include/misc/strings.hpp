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
#include <string>
#include <random>
#include <fstream>

#include "utils.hpp"

static std::string random_alphanumeric_string(uint64_t length)
{
    static std::string possible_chars = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint8_t> char_idx_distrib(0, possible_chars.size() - 1);
    
    std::string str_rand;
    str_rand.reserve(length);

    for (uint64_t i = 0; i < length; i++) {
        str_rand.push_back(possible_chars[char_idx_distrib(gen)]);
    }

    return str_rand;
}

static uint64_t malloc_count_peak_memory_usage(std::ifstream& leftg_file)
{
    std::string leftg_file_content;
    leftg_file.seekg(0, std::ios::end);
    no_init_resize(leftg_file_content, leftg_file.tellg());
    leftg_file.seekg(0, std::ios::beg);
    leftg_file.read((char*) &leftg_file_content[0], leftg_file_content.size());
    int32_t pos = 0;
    uint64_t cur_peak = 0;
    std::string str_cur_peak;

    while ((pos = leftg_file_content.find(", peak", pos)) != -1) {
        while (!('0' <= leftg_file_content[pos] && leftg_file_content[pos] <= '9')) {
            pos++;
        }

        while (('0' <= leftg_file_content[pos] && leftg_file_content[pos] <= '9') || leftg_file_content[pos] == '.') {
            if (leftg_file_content[pos] != '.') {
                str_cur_peak.push_back(leftg_file_content[pos]);
            }

            pos++;
        }

        cur_peak = std::max(cur_peak, (uint64_t)stol(str_cur_peak));
        str_cur_peak.clear();
    }

    return cur_peak;
}

template <typename inp_t>
static inp_t random_repetitive_input(
    uint64_t min_size, uint64_t max_size,
    typename inp_t::value_type min_sym = std::numeric_limits<typename inp_t::value_type>::min(),
    typename inp_t::value_type max_sym = std::numeric_limits<typename inp_t::value_type>::max()
) {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> prob_distrib(0.0, 1.0);
    std::uniform_int_distribution<uint64_t> input_size_distrib(min_size, max_size);
    std::uniform_int_distribution<typename inp_t::value_type> sym_distrib(min_sym, max_sym);

    uint64_t target_input_size = input_size_distrib(mt);
    enum construction_operation { new_symbol = 0, repetition = 1, run = 2 };
    double repetition_repetitiveness = prob_distrib(mt);
    double run_repetitiveness = prob_distrib(mt);

    std::uniform_int_distribution<uint64_t> repetition_length_distrib(
        1, std::max<double>(1.0, (repetition_repetitiveness * target_input_size) / 100));

    std::uniform_int_distribution<uint64_t> run_length_distrib(
        1, std::max<double>(1.0, ((run_repetitiveness * target_input_size) / 200)));

    std::discrete_distribution<uint8_t> next_operation_distrib({
        2 - (repetition_repetitiveness + run_repetitiveness),
        repetition_repetitiveness,
        run_repetitiveness });

    inp_t input;
    no_init_resize(input, target_input_size);
    input.clear();
    input.push_back(sym_distrib(mt));

    while (input.size() < target_input_size) {
        switch (next_operation_distrib(mt)) {
            case new_symbol: {
                input.push_back(sym_distrib(mt));
                break;
            }
            case repetition: {
                uint64_t repetition_length = std::min<uint64_t>(
                    target_input_size - input.size(), repetition_length_distrib(mt));
                uint64_t repstition_source = std::uniform_int_distribution<uint64_t>(0, input.size() - 1)(mt);

                for (uint64_t i = 0; i < repetition_length; i++)
                    input.push_back(input[repstition_source + i]);

                break;
            }
            case run: {
                uint64_t run_length = std::min<uint64_t>(
                    target_input_size - input.size(), run_length_distrib(mt));
                typename inp_t::value_type run_sym = sym_distrib(mt);

                for (uint64_t i = 0; i < run_length; i++)
                    input.push_back(run_sym);

                break;
            }
        }
    }

    return input;
}