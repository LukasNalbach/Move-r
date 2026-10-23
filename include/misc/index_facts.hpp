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

/**
 * @file index_facts.hpp
 * @brief the facts the info modes print about an index (n, r, k' of the move data structures, the number of phrases
 *        of the RLZ encoding, ...), one per line
 */

#ifndef MOVE_R_MISC_INDEX_FACTS_HPP
#define MOVE_R_MISC_INDEX_FACTS_HPP

#include <cstdint>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include <misc/log.hpp>

/** @brief collects the facts about an index, one per line (see print) */
class index_facts {
    std::vector<std::string> lines; // the facts in the order they were added

public:
    // adds a fact, e.g. "r = 7,992"
    void add(std::string line) { lines.push_back(std::move(line)); }

    // adds the fact "<name> = <value>"
    void add(const std::string& name, uint64_t value) { add(name + " = " + format_number(value)); }

    // prints the facts, one per line
    void print(std::ostream& out) const
    {
        for (const std::string& line : lines) out << line << "\n";
    }
};

#endif // MOVE_R_MISC_INDEX_FACTS_HPP
