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

// verify_support_of_iteration for the symbol type TEST_MOVE_R_SYMBOL, whose input container is std::string for char and
// std::vector<TEST_MOVE_R_SYMBOL> otherwise (this file is compiled once per input type, see CMakeLists.txt)

#include "test-move-r-common.hpp"

using instantiated_symbol_t = TEST_MOVE_R_SYMBOL;
using instantiated_input_t = std::conditional_t<std::is_same_v<instantiated_symbol_t, char>, std::string,
    std::vector<instantiated_symbol_t>>;

TEST_MOVE_R_VERIFY_SUPPORT_OF_ITERATION(, instantiated_symbol_t, instantiated_input_t)
