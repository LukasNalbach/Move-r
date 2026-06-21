/*
 * bi-directional r-index 
 *  the simpler implementation
 */

#ifndef INCLUDED_BR_INDEX_HPP
#define INCLUDED_BR_INDEX_HPP

#include "definitions.hpp"
#include "rle_string.hpp"
#include "sparse_sd_vector.hpp"
#include "permuted_lcp.hpp"
#include "utils.hpp"

namespace bri {

// sample maintained during the search
struct br_sample {
    /*
     * state variables for left_extension & right_extension
     * range: SA range of P
     * j: SA[p]
     * d: offset between starting position of the pattern & j
     * rangeR: correspondents to range in reversed text
     * len: current pattern length
     */
    range_t range, rangeR;
    ulint j, d, len;
    
    br_sample(): range(), rangeR() {}

    br_sample(range_t range_, 
              range_t rangeR_,
              ulint j_,
              ulint d_,
              ulint len_) 
              :
              range(range_),
              rangeR(rangeR_),
              j(j_),
              d(d_),
              len(len_) {}

    void set_values(range_t range_, 
                    range_t rangeR_,
                    ulint j_,
                    ulint d_,
                    ulint len_)
    {
        range = range_;
        rangeR = rangeR_;
        j = j_;
        d = d_;
        len = len_;
    }

    // the pattern does not exist
    bool is_invalid()
    {
        return (range.first > range.second) || (rangeR.first > rangeR.second);
    }

    // range size
    ulint size()
    {
        return range.second + 1 - range.first;
    }
};

template<
    class sparse_bitvector_t = sparse_sd_vector,
    class rle_string_t = rle_string_sd 
>
class br_index {

public:

    using triple = std::tuple<range_t, ulint, ulint>;

    br_index() {}

    /*
     * constructor. 
     * \param input: string on which br-index is built
     * \param sais: flag determining if we use SAIS for suffix sort. 
     *              otherwise we use divsufsort
     */
    br_index(std::string const& input, bool sais = true)
    {
        
        this->sais = sais;

        if (input.size() < 1)
        {

            std::cout << "Error: input string is empty" << std::endl;
            exit(1);

        }

        std::cout << "Text length = " << input.size() << std::endl << std::endl;

        std::cout << "(1/4) Remapping alphabet ... " << std::flush;
        
        // build RLBWT

        // configure & build indexes for sufsort & plcp
        sdsl::cache_config cc;

        // remap alphabet
        remap = std::vector<uchar>(256,0);
        remap_inv = std::vector<uchar>(256,0);
        {
            sigma = 1;
            std::vector<ulint> freqs(256,0);
            for (size_t i = 0; i < input.size(); ++i)
            {
                if (freqs[(uchar)input[i]]++ == 0) sigma++;
                if (sigma >= 255)
                {
                    std::cout << "Error: alphabet cannot be remapped (overflow)" << std::endl;
                    exit(1);
                }
            }
            uchar new_c = 2; // avoid reserved chars
            for (ulint c = 2; c < 256; ++c)
            {
                if (freqs[(uchar)c] != 0)
                {
                    remap[(uchar)c] = new_c;
                    remap_inv[new_c++] = (uchar)c;
                }
            }
        }

        std::cout << "done." << std::endl;
        std::cout << "(2/4) Building BWT, BWT^R, PLCP and computing SA samples";
        if (sais) std::cout << " (SA-SAIS) ... " << std::flush;
        else std::cout << " (DIVSUFSORT) ... " << std::flush;

        // remap input text
        sdsl::int_vector<8> text(input.size());
        for (size_t i = 0; i < input.size(); ++i)
            text[i] = remap[(uchar)input[i]];

        sdsl::append_zero_symbol(text);

        // cache text
        sdsl::store_to_cache(text, sdsl::conf::KEY_TEXT, cc);
        sdsl::construct_config::byte_algo_sa = sais ? sdsl::SE_SAIS : sdsl::LIBDIVSUFSORT;
        
        // cache SA
        sdsl::construct_sa<8>(cc);
        // cache ISA 
        sdsl::construct_isa(cc);

        
        sdsl::int_vector_buffer<> sa(sdsl::cache_file_name(sdsl::conf::KEY_SA, cc));
        last_SA_val = sa[sa.size()-1];
        auto bwt_and_samples = sufsort(text,sa);

        plcp = permuted_lcp<>(cc);

        // remove cache of text and SA
        sdsl::remove(sdsl::cache_file_name(sdsl::conf::KEY_TEXT, cc));
        sdsl::remove(sdsl::cache_file_name(sdsl::conf::KEY_SA, cc));
        sdsl::remove(sdsl::cache_file_name(sdsl::conf::KEY_ISA, cc));




        // configure & build reversed indexes for sufsort
        sdsl::cache_config ccR;

        sdsl::int_vector<8> textR(input.size());
        for (ulint i = 0; i < input.size(); ++i)
            textR[i] = remap[(uchar)input[input.size()-1-i]];

        sdsl::append_zero_symbol(textR);

        // cache textR
        sdsl::store_to_cache(textR, sdsl::conf::KEY_TEXT, ccR);
        sdsl::construct_config::byte_algo_sa = sais ? sdsl::SE_SAIS : sdsl::LIBDIVSUFSORT;
        
        // cache SAR
        sdsl::construct_sa<8>(ccR);
        // cache ISAR
        //sdsl::construct_isa(ccR);

        sdsl::int_vector_buffer<> saR(sdsl::cache_file_name(sdsl::conf::KEY_SA, ccR));
        auto bwt_and_samplesR = sufsort(textR,saR);

        // plcp is not needed in the reversed case

        // remove cache of textR and SAR
        sdsl::remove(sdsl::cache_file_name(sdsl::conf::KEY_TEXT, ccR));
        sdsl::remove(sdsl::cache_file_name(sdsl::conf::KEY_SA, ccR));





        std::string& bwt_s = std::get<0>(bwt_and_samples);
        std::vector<range_t>& samples_first_vec = std::get<1>(bwt_and_samples);
        std::vector<range_t>& samples_last_vec = std::get<2>(bwt_and_samples);

        std::string& bwt_sR = std::get<0>(bwt_and_samplesR);
        std::vector<range_t>& samples_first_vecR = std::get<1>(bwt_and_samplesR);
        std::vector<range_t>& samples_last_vecR = std::get<2>(bwt_and_samplesR);

        std::cout << "done.\n(3/4) Run-length encoding BWT ... " << std::flush;


        // run length compression on BWT and BWTR
        bwt = rle_string_t(bwt_s);
        bwtR = rle_string_t(bwt_sR);

        // build F column (common between text and textR)
        F = std::vector<ulint>(256,0);

        for (uchar c : bwt_s) 
            F[c]++;

        for (ulint i = 255; i > 0; --i) 
            F[i] = F[i-1];

        F[0] = 0;

        for(ulint i = 1; i < 256; ++i) 
            F[i] += F[i-1];


        // remember BWT position of terminator
		for(ulint i = 0; i < bwt_s.size(); ++i)
			if(bwt_s[i]==TERMINATOR)
				terminator_position = i;
        
        for(ulint i = 0; i < bwt_sR.size(); ++i)
			if(bwt_sR[i]==TERMINATOR)
				terminator_positionR = i;

        assert(input.size() + 1 == bwt.size());

        std::cout << "done." << std::endl << std::endl;


        r = bwt.number_of_runs();
        rR = bwtR.number_of_runs();

        assert(samples_first_vec.size() == r);
        assert(samples_last_vec.size() == r);

        assert(samples_first_vecR.size() == rR);
        assert(samples_last_vecR.size() == rR);

        int log_r = bitsize(r);
        int log_rR = bitsize(rR);
        int log_n = bitsize(bwt.size());

        std::cout << "Number of BWT equal-letter runs: r = " << r << std::endl;
		std::cout << "Rate n/r = " << double(bwt.size())/r << std::endl;
		std::cout << "log2(r) = " << std::log2(double(r)) << std::endl;
		std::cout << "log2(n/r) = " << std::log2(double(bwt.size())/r) << std::endl;

        std::cout << "Number of BWT^R equal-letter runs: rR = " << rR << std::endl << std::endl;

        // Phi, Phi inverse is needed only in forward case
        std::cout << "(4/4) Building predecessor for toehold lemma & Phi/Phi^{-1} function ..." << std::flush;

        
        samples_last = sdsl::int_vector<>(r,0,log_n);
        samples_first = sdsl::int_vector<>(r,0,log_n);
        
        samples_firstR = sdsl::int_vector<>(rR,0,log_n);
        samples_lastR = sdsl::int_vector<>(rR,0,log_n);

        for (ulint i = 0; i < r; ++i)
        {
            samples_last[i] = samples_last_vec[i].first;
            samples_first[i] = samples_first_vec[i].first;
        }
        for (ulint i = 0; i < rR; ++i)
        {
            samples_lastR[i] = samples_last_vecR[i].first;
            samples_firstR[i] = samples_first_vecR[i].first;
        }

        // sort samples of first positions in runs according to text position
        std::sort(samples_first_vec.begin(), samples_first_vec.end());
        // sort samples of last positions in runs according to text position
        std::sort(samples_last_vec.begin(), samples_last_vec.end());

        // build Elias-Fano predecessor
        {
            std::vector<bool> first_bv(bwt_s.size(),false);
            for (auto p: samples_first_vec)
            {
                assert(p.first < first_bv.size());
                first_bv[p.first] = true;
            }
            first = sparse_bitvector_t(first_bv);
        }
        {
            std::vector<bool> last_bv(bwt_s.size(),false);
            for (auto p: samples_last_vec)
            {
                assert(p.first < last_bv.size());
                last_bv[p.first] = true;
            }
            last = sparse_bitvector_t(last_bv);
        }

        assert(first.rank(first.size()) == r);
        assert(last.rank(last.size()) == r);

        first_to_run = sdsl::int_vector<>(r,0,log_r);

        last_to_run = sdsl::int_vector<>(r,0,log_r);

        // construct first_to_run
        for (ulint i = 0; i < samples_first_vec.size(); ++i)
        {
            first_to_run[i] = samples_first_vec[i].second;
        }

        // construct last_to_run
        for (ulint i = 0; i < samples_last_vec.size(); ++i)
        {
            last_to_run[i] = samples_last_vec[i].second;
        }

        // release ISA cache
        //sdsl::remove(sdsl::cache_file_name(sdsl::conf::KEY_ISA, cc));
        //sdsl::remove(sdsl::cache_file_name(sdsl::conf::KEY_ISA, ccR));

        std::cout << " done. " << std::endl << std::endl;

    }

    /*
     * get full BWT range
     */
    range_t full_range()
    {
        return {0,bwt_size()-1};
    }

    /*
     * rn: BWT range of a string P
     * c:  remapped character
     * returns: BWT range of cP
     */
    range_t LF(range_t rn, uchar c)
    {

        if ((c == 255 && F[c] == bwt.size()) || F[c] >= F[c+1]) return {1,0};

        ulint c_before = bwt.rank(rn.first, c);

        ulint c_inside = bwt.rank(rn.second+1,c) - c_before;

        if (c_inside == 0) return {1,0};

        ulint lb = F[c] + c_before;

        return {lb, lb + c_inside - 1};
    
    }

    /*
     * rn: BWT^R range of a string P
     * c:  remapped character
     * returns: BWT^R range of cP
     */
    range_t LFR(range_t rn, uchar c)
    {

        if ((c == 255 && F[c] == bwt.size()) || F[c] >= F[c+1]) return {1,0};

        ulint c_before = bwtR.rank(rn.first, c);

        ulint c_inside = bwtR.rank(rn.second+1,c) - c_before;

        if (c_inside == 0) return {1,0};

        ulint lb = F[c] + c_before;

        return {lb, lb + c_inside - 1};

    }

    /*
     * Phi function
     * get SA[i] from SA[i+1]
     */
    ulint Phi(ulint i)
    {
        assert(i != bwt.size() - 1);

        ulint jr = first.predecessor_rank_circular(i);

        assert(jr <= r - 1);

        ulint k = first.select(jr);

        assert(jr < r - 1 || k == bwt.size() - 1);

        // distance from predecessor
        ulint delta = k < i ? i - k : i + 1;

        // check if Phi(SA[0]) is not called
        assert(first_to_run[jr] > 0);

        ulint prev_sample = samples_last[first_to_run[jr]-1];

        return (prev_sample + delta) % bwt.size();
    }
    /*
     * Phi inverse
     * get SA[i] from SA[i-1]
     */
    ulint PhiI(ulint i)
    {
        assert(i != last_SA_val);

        ulint jr = last.predecessor_rank_circular(i);

        assert(jr <= r - 1);

        ulint k = last.select(jr);

        assert(jr < r - 1 || k == bwt.size() - 1);

        // distance from predecessor
        ulint delta = k < i ? i - k : i + 1;

        // check if Phi(SA[0]) is not called
        assert(last_to_run[jr] < r-1);

        ulint prev_sample = samples_first[last_to_run[jr]+1];

        return (prev_sample + delta) % bwt.size();
    }

    ulint LF(ulint i)
    {
        auto c = bwt[i];
        return F[c] + bwt.rank(i,c);
    }

    ulint LFR(ulint i)
    {
        auto c = bwtR[i];
        return F[c] + bwtR.rank(i,c);
    }

    /*
     * inverse of LF (known as Psi)
     */
    ulint FL(ulint i)
    {

        // i-th character in first BWT column F
        auto c = F_at(i);

        // j: occurrences of c before i
        ulint j = i - F[c];

        return bwt.select(j,(uchar)c);

    }

    ulint FLR(ulint i)
    {

        // i-th character in first BWT column F
        auto c = F_at(i);

        // j: occurrences of c before i
        ulint j = i - F[c];

        return bwtR.select(j,(uchar)c);
        
    }

    /*
     * character of position i in column F
     */
    uchar F_at(ulint i)
    {

        ulint c = (std::upper_bound(F.begin(),F.end(),i) - F.begin()) - 1;
        assert(c < 256);
        assert(i >= F[c]);

        return (uchar)c;

    }

    /*
     * return BWT range of original char c (not remapped)
     */
    range_t get_char_range(uchar c)
    {
        // replace c with internal representation
        c = remap[c];

        if ((c == 255 && F[c] == bwt_size()) || F[c] >= F[c+1]) return {1,0};

        ulint lb = F[c];
        ulint rb = bwt_size() - 1;

        if (c < 255) rb = F[c+1] - 1;

        return {lb,rb};

    }

    /*
     * get a sample corresponding to an empty string
     */
    br_sample get_initial_sample(bool right=true)
    {
        if (!right) {
            return br_sample(full_range(), // entire SA range
                            full_range(), // entire SAR range
                            (samples_last[r-1]+1) % bwt.size(), // arbitrary sample
                            0,            // offset 0
                            0);           // null pattern
        }
        else 
        {
            return br_sample(full_range(), // entire SA range
                            full_range(), // entire SAR range
                            (samples_lastR[rR-1]+1) % bwtR.size(), // arbitrary sample
                            0,            // offset 0
                            0);           // null pattern
        }
    }

    /*
     * search the pattern cP (P:the current pattern)
     * returns SA&SAR range corresponding to cP
     * 
     * full version, update ranges for SA and SAR simultaneously
     * 
     * assumes c is original char (not remapped)
     */
    br_sample left_extension(uchar c, br_sample const& prev_sample)
    {
        // replace c with internal representation
        c = remap[c];

        br_sample sample(prev_sample);

        // get SA range of cP
        sample.range = LF(prev_sample.range,c);

        // pattern cP was not found
        if (sample.is_invalid()) return sample;

        // cP and aP occurs for some a s.t. a != c
        if (prev_sample.range.second - prev_sample.range.first != 
            sample.range.second      - sample.range.first)
        {
            
            // accumulated occ of aP (for any a s.t. a < c)
            ulint acc = 0;
            for (ulint a = 1; a < c; ++a)
            {
                range_t smaller_range = LF(prev_sample.range,(uchar)a);
                acc += (smaller_range.second+1) - smaller_range.first;
            }
            // get SAR range of (cP)^R
            sample.rangeR.second = sample.rangeR.first + acc + sample.range.second - sample.range.first;
            sample.rangeR.first = sample.rangeR.first + acc;


            // fint last c in range and get its sample
            // there must be at least one c due to the previous if clause
            ulint rnk = bwt.rank(prev_sample.range.second+1,c);
            assert(rnk > 0);

            // update p by corresponding BWT position
            ulint p = bwt.select(rnk-1,c);
            assert(p >= prev_sample.range.first && p <= prev_sample.range.second);

            // run number of position p
            ulint run_of_p = bwt.run_of_position(p);

            // update j by SA[p]
            if (bwt[prev_sample.range.second] == c)
                sample.j = samples_first[run_of_p];
            else
                sample.j = samples_last[run_of_p];

            // reset d
            sample.d = 0;

        }
        else // only c precedes P 
        {
            // increment offset by 1
            sample.d++;
        }
        sample.len++;
        return sample;
    }

    /*
     * search the pattern Pc (P:the current pattern)
     * return SAR&SA range corresponding to Pc
     * 
     * full version, updates ranges for SAR and SA simultaneously
     * 
     * assumes c is original char (not remapped)
     */
    br_sample right_extension(uchar c, br_sample const& prev_sample)
    {
        // replace c with internal representation
        c = remap[c];

        br_sample sample(prev_sample);

        // get SAR range of Pc
        sample.rangeR = LFR(prev_sample.rangeR,c);

        // pattern Pc was not found
        if (sample.is_invalid()) return sample;

        // Pc and Pa occurs for some a s.t. a != c
        if (prev_sample.rangeR.second - prev_sample.rangeR.first != 
            sample.rangeR.second      - sample.rangeR.first)
        {

            // accumulated occ of Pa (for any a s.t. a < c)
            ulint acc = 0;
            for (ulint a = 1; a < c; ++a)
            {
                range_t smaller_rangeR = LFR(prev_sample.rangeR,(uchar)a);
                acc += (smaller_rangeR.second+1) - smaller_rangeR.first;
            }
            // get SA range of Pc
            sample.range.second = sample.range.first + acc + sample.rangeR.second - sample.rangeR.first; 
            sample.range.first = sample.range.first + acc;


            // fint last c in range and get its sample
            // there must be at least one c due to the previous if clause
            ulint rnk = bwtR.rank(prev_sample.rangeR.second+1,c);
            assert(rnk > 0);

            // update pR by corresponding BWTR position
            ulint pR = bwtR.select(rnk-1,c);
            assert(pR >= prev_sample.rangeR.first && pR <= prev_sample.rangeR.second);

            // run number of position pR
            ulint run_of_pR = bwtR.run_of_position(pR);

            // j = SA[p]
            if (bwtR[prev_sample.rangeR.second] == c)
                sample.j = bwt.size()-2-samples_firstR[run_of_pR];
            else
                sample.j = bwt.size()-2-samples_lastR[run_of_pR];

            // reset d
            sample.d = sample.len;
        }
        /*else // only c follows P
        {
            // increment offset by 1, unncecessary for locate
            sample.dR++; 
        }*/
        sample.len++;
        return sample;
    }

    /*
     * backward search P[left...right]
     */
    br_sample backward_search(std::string const& pattern, ulint left, ulint right, br_sample const& sample)
    {
        br_sample res(sample);
        for (ulint i = right + 1; i-- > left; )
        {
            res = left_extension(pattern[i],res);
            if (res.is_invalid()) return res;
        }
        return res;
    }

    /*
     * forward search P[left...right]
     */
    br_sample forward_search(std::string const& pattern, ulint left, ulint right, br_sample const& sample)
    {
        br_sample res(sample);
        for (ulint i = left; i <= right; ++i)
        {
            res = right_extension(pattern[i],res);
            if (res.is_invalid()) return res;
        }
        return res;
    }

    /*
     * count occurrences of current pattern P
     */
    ulint count_sample(br_sample const& sample)
    {
        return (sample.range.second + 1) - sample.range.first;
    }

    ulint count_samples(std::unordered_map<range_t,br_sample,range_hash> const& samples)
    {
        ulint res = 0;
		for (auto it = samples.begin(); it != samples.end(); ++it)
		{
			res += count_sample(it->second);
		}
		return res;
    }

    /*
     * locate occurrences of current pattern P
     * return them as std::vector
     * (space consuming if result is big)
     */
    std::vector<ulint> locate_sample(br_sample const& sample)
    {
        assert(sample.j >= sample.d);

        ulint sa = sample.j - sample.d;
        ulint pos = sa;

        std::vector<ulint> res;
        res.reserve(sample.range.second + 1 - sample.range.first);

        res.push_back(pos);

        while (plcp[pos] >= sample.len) 
        {
            pos = Phi(pos);
            res.push_back(pos);
        }
        pos = sa;
        while (true)
        {
            if (pos == last_SA_val) break;
            pos = PhiI(pos);
            if (plcp[pos] < sample.len) break;
            res.push_back(pos);
        }

        return res;
    }

    /*
     * locate occurrences of the current pattern P, centered on the toehold occurrence.
     *
     * Returns the suffix-array INDEX (within sample.range) of the toehold occurrence,
     * sets first_pos to its text position (SA value), and fills:
     *   left  = text positions at SA indices center-1, center-2, ..., range.first  (descending)
     *   right = text positions at SA indices center+1, center+2, ..., range.second (ascending)
     * so that the full occurrence set is {reverse(left), first_pos, right}. This mirrors
     * locate_sample but exposes the centered SA-index that move_r's edit-distance
     * approximate-pattern-matching enumeration needs.
     */
    ulint locate_sample_centered(br_sample const& sample, ulint& first_pos,
                                 std::vector<ulint>& left, std::vector<ulint>& right)
    {
        assert(sample.j >= sample.d);

        ulint sa = sample.j - sample.d;
        first_pos = sa;
        left.clear();
        right.clear();

        // walk left (towards smaller SA indices): Phi(SA[i]) = SA[i-1]
        ulint pos = sa;
        while (plcp[pos] >= sample.len)
        {
            pos = Phi(pos);
            left.push_back(pos);
        }

        // the toehold occupies SA index range.first + (number of left steps)
        ulint center = sample.range.first + left.size();

        // walk right (towards larger SA indices): PhiI(SA[i]) = SA[i+1]
        pos = sa;
        while (true)
        {
            if (pos == last_SA_val) break;
            pos = PhiI(pos);
            if (plcp[pos] < sample.len) break;
            right.push_back(pos);
        }

        return center;
    }

    std::vector<ulint> locate_samples(std::vector<br_sample> const& samples)
    {
        std::vector<ulint> res;
		for (auto s: samples)
		{
			auto tmp = locate_sample(s);
			res.insert(res.end(),tmp.begin(),tmp.end());
		}
		return res;
    }

    std::vector<ulint> locate_samples(std::unordered_map<range_t,br_sample,range_hash> const& samples)
    {
        std::vector<ulint> res;
		for (auto it = samples.begin(); it != samples.end(); ++it)
		{
			auto tmp = locate_sample(it->second);
			res.insert(res.end(),tmp.begin(),tmp.end());
		}
		return res;
    }

    /*
     * count the number of a given pattern
     */
    ulint count(std::string const& pattern, bool right=false)
    {
        if (right) 
        {
            br_sample sample(get_initial_sample(true));
            for (size_t i = 0; i < pattern.size(); ++i)
            {
                sample = right_only(pattern[i],sample);
                if (sample.is_invalid()) return {};
            }
            return count_sample(sample);
        }
        else 
        {
            br_sample sample(get_initial_sample());
            for (size_t i = 0; i < pattern.size(); ++i)
            {
                sample = left_only(pattern[pattern.size()-1-i],sample);
                if (sample.is_invalid()) return {};
            }
            return count_sample(sample);
        }
    }

    /*
     * locate occurrences of a given pattern
     */
    std::vector<ulint> locate(std::string const& pattern, bool right=false)
    {
        if (right) 
        {
            br_sample sample(get_initial_sample(true));
            for (size_t i = 0; i < pattern.size(); ++i)
            {
                sample = right_only(pattern[i],sample);
                if (sample.is_invalid()) return {};
            }
            return locate_sample(sample);
        }
        else 
        {
            br_sample sample(get_initial_sample());
            for (size_t i = 0; i < pattern.size(); ++i)
            {
                sample = left_only(pattern[pattern.size()-1-i],sample);
                if (sample.is_invalid()) return {};
            }
            return locate_sample_backward(sample);
        }
    }

    ulint count_with_mismatch(std::string const& pattern, ulint allowed_mis=0)
    {
        auto samples = search_with_mismatch(pattern,allowed_mis);
        ulint res = 0;
        for (auto s: samples)
        {
            res += count_sample(s.second);
        }
        return res;
    }

    std::vector<ulint> locate_with_mismatch(std::string const& pattern, ulint allowed_mis=0)
    {
        auto samples = search_with_mismatch(pattern,allowed_mis);
        return locate_samples(samples);
    }

    std::unordered_map<range_t,br_sample,range_hash> search_with_mismatch(std::string const& pattern, ulint allowed_mis=0)
    {
        std::unordered_map<range_t,br_sample,range_hash> res;
        ulint m = pattern.size();
        br_sample init_sample(get_initial_sample());

        if (allowed_mis == 0)
        {
            br_sample sample(backward_search(pattern,0,m-1,init_sample));
            if (sample.is_invalid()) return res;
            res[sample.range] = sample;
            return res;
        }

        ulint div = allowed_mis + 1;
        // divide pattern into div parts and search each part in advance
        for (ulint part = 0; part < div; ++part)
        {
            br_sample sample(backward_search(pattern,(part*m)/div,((part+1)*m)/div-1,init_sample));
            if (sample.is_invalid()) continue;
            if (part==div-1) 
            {
                backward_dfs(res,pattern,m,allowed_mis,(part*m)/div-1,((part+1)*m)/div-1,0,sample);
            }
            else 
            {
                forward_dfs(res,pattern,m,allowed_mis,(part*m)/div,((part+1)*m)/div,0,sample);
            }
        }

        return res;
    }

    std::unordered_map<range_t,br_sample,range_hash> seed_and_extend(std::string const& pattern, ulint m1, ulint m2, ulint allowed_mis=0)
    {
        // P[0,m1-1], P[m1,m2-1], P[m2,m-1]
        std::unordered_map<range_t,br_sample,range_hash> res;
        ulint m = pattern.size();
        br_sample init_sample(get_initial_sample());

        br_sample sample(backward_search(pattern,m1,m2-1,init_sample));
        if (sample.is_invalid()) return res;

        forward_dfs(res,pattern,m,allowed_mis,m1,m2,0,sample);

        return res;
    }

    void backward_dfs(std::unordered_map<range_t,br_sample,range_hash>& res, std::string const& pattern,
                    ulint m, ulint allowed_mis, ulint left_pos, ulint right_pos, ulint mis, br_sample prev_sample)
    {
        uchar c = remap[pattern[left_pos]];

        //std::cout << mis << " " << allowed_mis << "  ";

        if (mis == allowed_mis)
        {
            ulint acc = 0;

            br_sample sample(prev_sample);

            sample.range = LF(prev_sample.range,c);
            if (sample.is_invalid()) return;

            if (sample.range.second - sample.range.first == 
                prev_sample.range.second - prev_sample.range.first)
            {
                sample.d++;
            } 
            else 
            {
                for (ulint a = 1; a < c; ++a)
                {
                    range_t smaller_range = LF(prev_sample.range,(uchar)a);
                    acc += smaller_range.second + 1 - smaller_range.first;
                }

                sample.rangeR.first = sample.rangeR.first + acc;
                sample.rangeR.second = sample.rangeR.first + sample.range.second - sample.range.first;
                ulint rnk = bwt.rank(prev_sample.range.second+1,c);
                ulint p = bwt.select(rnk-1,c);
                ulint run_of_p = bwt.run_of_position(p);
                if (bwt[prev_sample.range.second] == c)
                    sample.j = samples_first[run_of_p];
                else
                    sample.j = samples_last[run_of_p];
                sample.d = 0;
            }
            sample.len++;

            if (left_pos == 0)
            {
                res[sample.range] = sample;
            }
            else 
            {
                backward_dfs(res,pattern,m,allowed_mis,left_pos-1,right_pos,mis,sample);
            }
            
        } 
        else // mis < allowed_mis 
        {
            ulint acc = 0;
        
            for (ulint a = 1; a < sigma+1; ++a)
            {
                br_sample sample(prev_sample);

                sample.range = LF(prev_sample.range,(uchar)a);
                if (sample.is_invalid()) continue;

                if (a == 1)
                {
                    acc++;
                    continue;
                }
                
                if (sample.range.second - sample.range.first == 
                    prev_sample.range.second - prev_sample.range.first)
                {
                    sample.d++;
                }
                else // sample.size() != prev_sample.size()
                {
                    sample.rangeR.first = sample.rangeR.first + acc;
                    sample.rangeR.second = sample.rangeR.first + sample.range.second - sample.range.first;
                    ulint rnk = bwt.rank(prev_sample.range.second+1,a);
                    ulint p = bwt.select(rnk-1,a);
                    ulint run_of_p = bwt.run_of_position(p);
                    if (bwt[prev_sample.range.second] == a)
                        sample.j = samples_first[run_of_p];
                    else
                        sample.j = samples_last[run_of_p];
                    sample.d = 0;
                }
                sample.len++;
                acc += sample.range.second + 1 - sample.range.first;

                if (c == a)
                {
                    if (left_pos == 0)
                    {
                        res[sample.range] = sample;
                    }
                    else 
                    {
                        backward_dfs(res,pattern,m,allowed_mis,left_pos-1,right_pos,mis,sample);
                    }
                }
                else // c != a
                {
                    if (left_pos == 0)
                    {
                        res[sample.range] = sample;
                    }
                    else
                    {
                        backward_dfs(res,pattern,m,allowed_mis,left_pos-1,right_pos,mis+1,sample);
                    }
                }
            }
        }   
    }

    void forward_dfs(std::unordered_map<range_t,br_sample,range_hash>& res, std::string const& pattern,
                    ulint m, ulint allowed_mis, ulint left_pos, ulint right_pos, ulint mis, br_sample prev_sample)
    {
        uchar c = remap[pattern[right_pos]];

        if (mis == allowed_mis)
        {
            ulint acc = 0;

            br_sample sample(prev_sample);

            sample.rangeR = LFR(prev_sample.rangeR,c);
            if (sample.is_invalid()) return;


            if (sample.rangeR.second - sample.rangeR.first != 
                prev_sample.rangeR.second - prev_sample.rangeR.first)
            {
                for (ulint a = 1; a < c; ++a)
                {
                    range_t smaller_rangeR = LFR(prev_sample.rangeR,(uchar)a);
                    acc += (smaller_rangeR.second + 1) - smaller_rangeR.first;
                }

                sample.range.first = sample.range.first + acc;
                sample.range.second = sample.range.first + sample.rangeR.second - sample.rangeR.first;
                ulint rnk = bwtR.rank(prev_sample.rangeR.second+1,c);
                ulint p = bwtR.select(rnk-1,c);
                ulint run_of_p = bwtR.run_of_position(p);
                if (bwtR[prev_sample.rangeR.second] == c)
                    sample.j = bwt.size()-2-samples_firstR[run_of_p];
                else
                    sample.j = bwt.size()-2-samples_lastR[run_of_p];
                sample.d = sample.len;
            }
            sample.len++;

            if (left_pos == 0 && right_pos >= m - 1)
            {
                res[sample.range] = sample;
            }
            else if (right_pos >= m-1)
            {
                backward_dfs(res,pattern,m,allowed_mis,left_pos-1,right_pos,mis,sample);
            }
            else
            {
                forward_dfs(res,pattern,m,allowed_mis,left_pos,right_pos+1,mis,sample);
            }
        } 
        else // mis < allowed_mis 
        {
            ulint acc = 0;
        
            for (ulint a = 1; a < sigma+1; ++a)
            {
                br_sample sample(prev_sample);

                sample.rangeR = LFR(prev_sample.rangeR,(uchar)a);
                if (sample.is_invalid()) continue;

                if (a == 1)
                {
                    acc++;
                    continue;
                }
                
                if (sample.rangeR.second - sample.rangeR.first != 
                    prev_sample.rangeR.second - prev_sample.rangeR.first)
                {
                    sample.range.first = sample.range.first + acc;
                    sample.range.second = sample.range.first + sample.rangeR.second - sample.rangeR.first;
                    ulint rnk = bwtR.rank(prev_sample.rangeR.second+1,a);
                    ulint p = bwtR.select(rnk-1,a);
                    ulint run_of_p = bwtR.run_of_position(p);
                    if (bwtR[prev_sample.rangeR.second] == a)
                        sample.j = bwt.size()-2-samples_firstR[run_of_p];
                    else
                        sample.j = bwt.size()-2-samples_lastR[run_of_p];
                    sample.d = sample.len;
                }
                sample.len++;
                acc += sample.range.second + 1 - sample.range.first;

                if (c == a)
                {
                    if (left_pos == 0 && right_pos >= m - 1)
                    {
                        res[sample.range] = sample;
                    }
                    else if (right_pos >= m - 1)
                    {
                        backward_dfs(res,pattern,m,allowed_mis,left_pos-1,right_pos,mis,sample);
                    }
                    else
                    {
                        forward_dfs(res,pattern,m,allowed_mis,left_pos,right_pos+1,mis,sample);
                    }
                }
                else // c != a
                {
                    if (left_pos == 0 && right_pos >= m - 1)
                    {
                        res[sample.range] = sample;
                    }
                    else if (right_pos >= m - 1)
                    {
                        backward_dfs(res,pattern,m,allowed_mis,left_pos-1,right_pos,mis+1,sample);
                    }
                    else
                    {
                        forward_dfs(res,pattern,m,allowed_mis,left_pos,right_pos+1,mis+1,sample);
                    }
                }
            }
        }
    }

    /*
     * get BWT[i] or BWT^R[i]
     */
    uchar bwt_at(ulint i, bool reversed = false)
    {
        if (!reversed) return remap_inv[bwt[i]];
        return remap_inv[bwtR[i]];
    }

    /*
     * get number of runs in BWT
     */
    ulint number_of_runs(bool reversed = false)
    {
        if (!reversed) return bwt.number_of_runs();
        return bwtR.number_of_runs();
    }

    /*
     * get position of terminator symbol in BWT
     */
    ulint get_terminator_position(bool reversed = false)
    {
        if (!reversed) return terminator_position;
        return terminator_positionR;
    }

    /*
     * get string representation of BWT
     */
    std::string get_bwt(bool reversed = false)
    {
        if (!reversed)
        {
            std::string res(bwt.to_string());
            for (size_t i = 0; i < res.size(); ++i)
                res[i] = remap_inv[(uchar)res[i]];
            return res;
        } else {
            std::string res(bwtR.to_string());
            for (size_t i = 0; i < res.size(); ++i)
                res[i] = remap_inv[(uchar)res[i]];
            return res;
        }
    }

    uint serialize(std::ostream& out)
    {
        ulint w_bytes = 0;

        out.write((char*)&sigma,sizeof(sigma));

        out.write((char*)remap.data(),256*sizeof(uchar));
        out.write((char*)remap_inv.data(),256*sizeof(uchar));

        out.write((char*)&terminator_position,sizeof(terminator_position));
        out.write((char*)&terminator_positionR,sizeof(terminator_positionR));
        out.write((char*)&last_SA_val,sizeof(last_SA_val));
        out.write((char*)F.data(),256*sizeof(ulint));

        w_bytes += sizeof(sigma)
                   + 256*sizeof(uchar)
                   + 256*sizeof(uchar)
                   + sizeof(terminator_position)
                   + sizeof(terminator_positionR)
                   + sizeof(last_SA_val)
                   + 256*sizeof(ulint);
        
        w_bytes += bwt.serialize(out);
        w_bytes += bwtR.serialize(out);

        w_bytes += samples_first.serialize(out);
        w_bytes += samples_last.serialize(out);

        w_bytes += first.serialize(out);
        w_bytes += first_to_run.serialize(out);

        w_bytes += last.serialize(out);
        w_bytes += last_to_run.serialize(out);

        w_bytes += samples_firstR.serialize(out);
        w_bytes += samples_lastR.serialize(out);

        w_bytes += plcp.serialize(out);

        return w_bytes;
    
    }

    void load(std::istream& in)
    {

        in.read((char*)&sigma,sizeof(sigma));

        remap = std::vector<uchar>(256);
        in.read((char*)remap.data(),256*sizeof(uchar));
        remap_inv = std::vector<uchar>(256);
        in.read((char*)remap_inv.data(),256*sizeof(uchar));
        
        in.read((char*)&terminator_position,sizeof(terminator_position));
        in.read((char*)&terminator_positionR,sizeof(terminator_positionR));
        in.read((char*)&last_SA_val,sizeof(last_SA_val));
        
        F = std::vector<ulint>(256);
        in.read((char*)F.data(),256*sizeof(ulint));

        bwt.load(in);
        bwtR.load(in);
        r = bwt.number_of_runs();
        rR = bwtR.number_of_runs();

        samples_first.load(in);
        samples_last.load(in);

        first.load(in);
        first_to_run.load(in);

        last.load(in);
        last_to_run.load(in);

        samples_firstR.load(in);
        samples_lastR.load(in);

        plcp.load(in);

    }

    /*
     * save index to "{path_prefix}.bri" file
     */
    void save_to_file(std::string const& path_prefix)
    {

        std::string path = path_prefix + ".bri";
        
        std::ofstream out(path);
        serialize(out);
        out.close();
    
    }

    /*
     * load index file from path
     */
    void load_from_file(std::string const& path)
    {

        std::ifstream in(path);
        load(in);
        in.close();

    }

    ulint text_size() { return bwt.size() - 1; }

    ulint bwt_size(bool reversed=false) { return bwt.size(); }

    uchar get_terminator() {
        return TERMINATOR;
    }

    /*
     * alphabet size (including the terminator); internal symbols are 1 (terminator)
     * and 2..sigma for the sigma-1 distinct original characters
     */
    ulint alphabet_size() const { return sigma; }

    /*
     * map an original character to its internal (remapped) symbol (0 if absent)
     */
    uchar char_to_internal(uchar c) const { return remap[c]; }

    /*
     * map an internal (remapped) symbol back to its original character
     */
    uchar internal_to_char(uchar i) const { return remap_inv[i]; }

    /*
     * get statistics
     */
    ulint print_space() 
    {

        std::cout << "text length           : " << bwt.size() << std::endl;
        std::cout << "alphabet size         : " << sigma << std::endl;
        std::cout << "number of runs in bwt : " << bwt.number_of_runs() << std::endl;
        std::cout << "numbef of runs in bwtR: " << bwtR.number_of_runs() << std::endl << std::endl;
        
        ulint tot_bytes = sizeof(sigma)
                        + 256*sizeof(uchar)
                        + 256*sizeof(uchar)
                        + sizeof(terminator_position)
                        + sizeof(terminator_positionR)
                        + sizeof(last_SA_val)
                        + 256*sizeof(ulint);
        
        tot_bytes += bwt.print_space();
        tot_bytes += bwtR.print_space();
        std::cout << "total space for BWT: " << tot_bytes << " bytes" << std::endl << std::endl;

        tot_bytes += plcp.print_space();

        std::ofstream out("/dev/null");

        ulint bytes = 0;

        
        bytes =  samples_first.serialize(out);
        tot_bytes += bytes;
        std::cout << "samples_first: " << bytes << " bytes" << std::endl;

        bytes =  samples_last.serialize(out);
        tot_bytes += bytes;
        std::cout << "samples_last: " << bytes << " bytes" << std::endl;


        bytes =  first.serialize(out);
        tot_bytes += bytes;
        std::cout << "first: " << bytes << " bytes" << std::endl;

        bytes =  first_to_run.serialize(out);
        tot_bytes += bytes;
        std::cout << "first_to_run: " << bytes << " bytes" << std::endl;


        bytes =  last.serialize(out);
        tot_bytes += bytes;
        std::cout << "last: " << bytes << " bytes" << std::endl;

        bytes =  last_to_run.serialize(out);
        tot_bytes += bytes;
        std::cout << "last_to_run: " << bytes << " bytes" << std::endl;


        bytes =  samples_firstR.serialize(out);
        tot_bytes += bytes;
        std::cout << "samples_firstR: " << bytes << " bytes" << std::endl;

        bytes =  samples_lastR.serialize(out);
        tot_bytes += bytes;
        std::cout << "samples_lastR: " << bytes << " bytes" << std::endl;

        std::cout << "<total space of br-index>: " << tot_bytes << " bytes" << std::endl << std::endl;
        std::cout << "<bits/symbol>            : " << (double) tot_bytes * 8 / (double) bwt.size() << std::endl;

        return tot_bytes;

    }

    /*
     * get space complexity
     */
    ulint get_space()
    {

        ulint tot_bytes = sizeof(sigma)
                        + 256*sizeof(uchar)
                        + 256*sizeof(uchar)
                        + sizeof(terminator_position)
                        + sizeof(terminator_positionR)
                        + sizeof(last_SA_val)
                        + 256*sizeof(ulint);

        tot_bytes += bwt.get_space();
        tot_bytes += bwtR.get_space();

        tot_bytes += plcp.get_space();

        std::ofstream out("/dev/null");

        tot_bytes += samples_first.serialize(out);
        tot_bytes += samples_last.serialize(out);

        tot_bytes += first.serialize(out);
        tot_bytes += first_to_run.serialize(out);

        tot_bytes += last.serialize(out);
        tot_bytes += last_to_run.serialize(out);

        tot_bytes += samples_firstR.serialize(out);
        tot_bytes += samples_lastR.serialize(out);

        return tot_bytes;

    }

private:

    /*
     * only updates range for SA
     * use when you only search backward
     * 
     * assumes sample is SA[r] if range is [l,r]
     * assumes c is original char (not remapped)
     */
    br_sample left_only(uchar c, br_sample const& prev_sample)
    {
        // replace c with internal representation
        c = remap[c];

        br_sample sample(prev_sample);

        // get SA range of cP
        sample.range = LF(prev_sample.range,c);

        // pattern cP was not found
        if (sample.is_invalid()) return sample;

        if (bwt[prev_sample.range.second] == c)
        {
            sample.d++;
        } 
        else 
        {
            ulint rnk = bwt.rank(prev_sample.range.second,c);
            rnk--;
            ulint k = bwt.select(rnk,c);
            ulint run_of_k = bwt.run_of_position(k);
            sample.j = samples_last[run_of_k];
            sample.d = 0;
        }
        sample.len++;

        return sample;
    }

    /*
     * only updates range for SAR
     * use when you only search forward
     * 
     * assumes sample is SAR[e] if range is [s,e]
     * assumes c is original char (not remapped)
     */
    br_sample right_only(uchar c, br_sample const& prev_sample)
    {
        // replace c with internal representation
        c = remap[c];

        br_sample sample(prev_sample);

        // get SAR range of Pc
        sample.rangeR = LFR(prev_sample.rangeR,c);

        // pattern Pc was not found
        if (sample.is_invalid()) return sample;

        if (bwtR[prev_sample.rangeR.second] != c)
        {
            ulint rnk = bwtR.rank(prev_sample.rangeR.second,c);
            rnk--;
            ulint k = bwtR.select(rnk,c);
            ulint run_of_k = bwtR.run_of_position(k);
            sample.j = bwt.size()-2-samples_lastR[run_of_k];
            sample.d = sample.len;
        }
        sample.len++;
        return sample;
    }

    /*
     * backward search P[left...right]
     * range for SAR is not updated
     */
    br_sample backward_only(std::string const& pattern, ulint left, ulint right, br_sample const& sample)
    {
        br_sample res(sample);
        for (ulint i = right + 1; i-- > left; )
        {
            res = left_only(pattern[i],res);
            if (res.is_invalid()) return res;
        }
        return res;
    }

    /*
     * forward search P[left...right]
     * range for SA is not updated
     */
    br_sample forward_only(std::string const& pattern, ulint left, ulint right, br_sample const& sample)
    {
        br_sample res(sample);
        for (ulint i = left; i <= right; ++i)
        {
            res = right_only(pattern[i],res);
            if (res.is_invalid()) return res;
        }
        return res;
    }

    /*
     * locate occurrences of current pattern P 
     * use while backward searching only (since it assumes sample is SA[r] when SA range is [l,r])
     * 
     * return them as std::vector
     * (space consuming if result is big)
     */
    std::vector<ulint> locate_sample_backward(br_sample const& sample)
    {
        ulint sa = sample.j - sample.d;
        ulint n_occ = sample.range.second + 1 - sample.range.first;

        std::vector<ulint> res(n_occ);
        if (n_occ > 0)
        {
            res[0] = sa;
            for (ulint i = 1; i < n_occ; ++i)
            {
                sa = Phi(sa);
                res[i] = sa;
            }
        }
        return res;
    }

    /*
     * builds BWT
     */
    std::tuple<std::string, std::vector<range_t>, std::vector<range_t> > 
    sufsort(sdsl::int_vector<8>& text, sdsl::int_vector_buffer<>& sa)
    {
        std::string bwt_s;
        std::vector<range_t> samples_first;
        std::vector<range_t> samples_last;

        {
            for (ulint i = 0; i < sa.size(); ++i)
            {
                auto x = sa[i];

                assert(x <= text.size());

                if (x > 0) 
                    bwt_s.push_back((uchar)text[x-1]);
                else 
                    bwt_s.push_back(TERMINATOR);
                
                // insert samples at beginnings of runs
                if (i > 0)
                {
                    if (i==1 || (i>1 && bwt_s[i-1] != bwt_s[i-2]))
                    {
                        samples_first.push_back({
                            sa[i-1] > 0
                            ? sa[i-1] - 1
                            : sa.size() - 1,
                            samples_first.size()
                        });
                    }
                    if (i==sa.size()-1 && bwt_s[i] != bwt_s[i-1])
                    {
                        samples_first.push_back({
                            sa[i] > 0
                            ? sa[i] - 1
                            : sa.size() - 1,
                            samples_first.size()
                        });
                    }
                }

                // insert samples at ends of runs
                if (i > 0)
                {
                    if (bwt_s[i-1] != bwt_s[i])
                    {
                        samples_last.push_back({
                            sa[i-1] > 0
                            ? sa[i-1] - 1
                            : sa.size() - 1,
                            samples_last.size()
                        });
                    }
                    if (i == sa.size()-1)
                    {
                        samples_last.push_back({
                            sa[i] > 0
                            ? sa[i] - 1
                            : sa.size() - 1,
                            samples_last.size()
                        });
                    }
                }
            }
        }

        return std::tuple<std::string, std::vector<range_t>, std::vector<range_t> >
            (bwt_s, samples_first, samples_last);
    }

    static const uchar TERMINATOR = 1;

    bool sais = true;

    /*
     * sparse RLBWT for text & textR
     */

    // alphabet remapper
    std::vector<uchar> remap;
    std::vector<uchar> remap_inv;
    ulint sigma;

    // accumulated number of characters in lex order
    std::vector<ulint> F;
    
    // RLBWT
    rle_string_t bwt;
    ulint terminator_position = 0;
    ulint last_SA_val = 0;
    ulint r = 0;

    // RLBWT^R
    rle_string_t bwtR;
    ulint terminator_positionR = 0;
    ulint rR = 0;

    // needed for left_extension
    sdsl::int_vector<> samples_first;
    sdsl::int_vector<> samples_last;
    
    // needed for Phi (SA[i] -> SA[i-1])
    sparse_bitvector_t first;
    sdsl::int_vector<> first_to_run;
    
    // needed for Phi^{-1} (SA[i] -> SA[i+1])
    sparse_bitvector_t last;
    sdsl::int_vector<> last_to_run;

    // needed for right_extension
    sdsl::int_vector<> samples_firstR;
    sdsl::int_vector<> samples_lastR;

    // needed for determining the end of locate
    permuted_lcp<> plcp;

};

};

#endif /* INCLUDED_BR_INDEX_HPP */