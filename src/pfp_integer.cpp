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

#include <algorithm>
#include <bit>
#include <cstring>
#include <memory>
#include <vector>

#include <misc/files.hpp>
#include <move_r/construction/pfp_integer.hpp>

namespace {

// the text of a pfp_integer_text_t with symbols of type sym_t
template <typename sym_t>
class integer_source : public pfp_rlbwt::text_source<sym_t> {
public:
    integer_source(const pfp_integer_text_t& input, uint64_t offset, bool reverse)
        : input(input), offset(offset), reverse(reverse)
    { }

    uint64_t size() const override { return input.length; }
    std::unique_ptr<typename pfp_rlbwt::text_source<sym_t>::reader> new_reader() const override;

private:
    class reader_t : public pfp_rlbwt::text_source<sym_t>::reader {
    public:
        explicit reader_t(const integer_source& src) : src(src)
        {
            if (src.input.text == nullptr) file.open_file(src.input.file);
        }

        void read(uint64_t pos, uint64_t len, sym_t* dst) override;

    private:
        const integer_source& src;
        direct_reader file; // reads the file (if the text is not in memory)
        std::vector<uint64_t> effective; // the effective symbols of the symbols read from the file
    };

    const pfp_integer_text_t& input;
    uint64_t offset;
    bool reverse;
};

template <typename sym_t>
std::unique_ptr<typename pfp_rlbwt::text_source<sym_t>::reader> integer_source<sym_t>::new_reader() const
{
    return std::make_unique<reader_t>(*this);
}

template <typename sym_t>
void integer_source<sym_t>::reader_t::read(uint64_t pos, uint64_t len, sym_t* dst)
{
    const uint8_t width = src.input.width;
    const uint64_t beg = src.reverse ? src.input.length - pos - len : pos;
    const char* symbols = src.input.text != nullptr ? src.input.text + beg * width
        : (const char*) file.read(beg * width, (beg + len) * width);

    if (src.input.alphabet != nullptr) {
        if (effective.size() < len) effective.resize(len);
        src.input.alphabet->map_symbols(symbols, len, effective.data());
        for (uint64_t i = 0; i < len; i++) dst[i] = sym_t(effective[i] + src.offset);
    } else {
        for (uint64_t i = 0; i < len; i++) {
            uint64_t sym = 0;
            std::memcpy(&sym, symbols + i * width, width);
            dst[i] = sym_t(sym + src.offset);
        }
    }

    if (src.reverse) std::reverse(dst, dst + len);
}

template <typename sym_t>
void build(const pfp_integer_text_t& input, uint64_t offset, bool reverse, const std::string& prefix,
    const pfp_rlbwt::params& prm)
{
    const integer_source<sym_t> text(input, offset, reverse);
    pfp_rlbwt::build(text, prefix, prm);
}

} // namespace

void build_rlbwt_pfp_integer(const pfp_integer_text_t& input, uint64_t offset, bool reverse,
    const std::string& prefix, const pfp_rlbwt::params& prm)
{
    // the symbol type of pfp-rlbwt with the fewest bytes that hold the largest symbol
    switch (byte_width(input.max_symbol + offset)) {
        #define MOVE_R_PFP_SYMBOL_CASE(sym_t) case sizeof(sym_t): build<sym_t>(input, offset, reverse, prefix, prm); break;
        PFP_RLBWT_FOR_EACH_SYMBOL_TYPE(MOVE_R_PFP_SYMBOL_CASE)
        #undef MOVE_R_PFP_SYMBOL_CASE
    }
}
