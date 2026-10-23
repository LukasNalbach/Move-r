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
 * @file log_structures.hpp
 * @brief helpers of the log_data_structures functions: the limits of the logged entries and trees logged with
 *        indentation
 */

#ifndef MOVE_R_MISC_LOG_STRUCTURES_HPP
#define MOVE_R_MISC_LOG_STRUCTURES_HPP

#include <algorithm>
#include <cstdint>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include <misc/aligned_table.hpp>
#include <misc/log.hpp>

// default depth up to which the log_data_structures functions log the nodes of a tree (a root has depth 0)
inline constexpr uint64_t default_logged_tree_depth = 8;

// default maximum number of nodes or entries a log_data_structures function logs per table or tree
inline constexpr uint64_t default_logged_entries = 256;

/**
 * @brief logs the nodes of a tree, one per line, indented by two spaces per level; nodes deeper than a maximum depth
 *        and nodes beyond a maximum number are left out, which limits() states at the end
 */
class tree_log {
    uint64_t _max_depth; // the nodes up to this depth are logged
    uint64_t _max_nodes; // at most this many nodes are logged
    uint64_t _num_logged = 0; // the number of logged nodes
    bool _below_max_depth = false; // true <=> a node has not been logged because it is deeper than _max_depth
    bool _beyond_max_nodes = false; // true <=> a node has not been logged because _max_nodes nodes have been logged

public:
    tree_log(uint64_t max_depth, uint64_t max_nodes) : _max_depth(max_depth), _max_nodes(max_nodes) {}

    /**
     * @brief logs a node, if it is not deeper than the maximum depth and fewer than the maximum number of nodes have
     *        been logged
     * @param depth the depth of the node
     * @param text the line of the node
     * @return whether the node has been logged (its children are logged only then)
     */
    bool node(uint64_t depth, const std::string& text)
    {
        if (depth > _max_depth) {
            _below_max_depth = true;
            return false;
        }

        if (_num_logged == _max_nodes) {
            _beyond_max_nodes = true;
            return false;
        }

        _num_logged++;
        std::cout << std::string(2 * depth, ' ') << text << std::endl;
        return true;
    }

    // logs the lines that state which nodes were left out
    void limits() const
    {
        if (_below_max_depth) std::cout << "the nodes below depth " << _max_depth << " are not logged" << std::endl;
        if (_beyond_max_nodes)
            std::cout << "only the first " << format_number(_max_nodes) << " nodes are logged" << std::endl;
    }
};

/**
 * @brief logs a table of at most max_entries columns (see aligned_table) and a line, if there are more
 * @param num_entries the number of entries
 * @param max_entries the maximum number of logged entries
 * @param add_rows function called as add_rows(table, number of logged entries), which adds the rows
 */
template <typename add_rows_t>
inline void log_limited_table(uint64_t num_entries, uint64_t max_entries, add_rows_t add_rows)
{
    aligned_table table;
    add_rows(table, std::min(num_entries, max_entries));
    table.print();
    if (num_entries > max_entries) std::cout << "only the first " << format_number(max_entries) << " of "
        << format_number(num_entries) << " entries are logged" << std::endl;
}

#endif // MOVE_R_MISC_LOG_STRUCTURES_HPP
