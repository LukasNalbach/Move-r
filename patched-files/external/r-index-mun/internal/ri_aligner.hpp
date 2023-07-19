#ifndef R_INDEX_S_F_H
#define R_INDEX_S_F_H

#include <queue>
#include <definitions.hpp>
#include <r_index.hpp>

// inexact_(locate|count) routines (D_algo, etc) are lifted from BWA:
// doi:10.1093/bioinformatics/btp698

namespace ri_mun {
    struct ri_opts_t {
        ri_opts_t() {};
        ulint max_hits = (ulint)-1;
        ulint max_range = (ulint)-1;
        uint z = 0;
        std::string algo = "";
    };

    struct partial_aln {
        partial_aln() {}
        uint mm = 0; // mismatches
        uint go = 0; // gapopens
        uint ge = 0; // gapextends
        uint len = 0;
        range_t r, pr; // range, previous range
        ulint k=0,j=0; // sample's bwt pos, SA pos
        uint score = (uint)-1;
        // char CIGAR[512]; // TODO: store this as well?
    };

    // these tuple types are for the stack/recursion routines used for inexact matching
    typedef std::tuple<uint,uint,range_t> iec_t; // inexact_count: i, z, range
    // typedef std::tuple<uchar,uint,uint,range_t,range_t,ulint,ulint> iel_t; // inexact_locate: c,i,z,prev_range,range,k,j
    typedef std::tuple<uchar,uint,uint,partial_aln> iel_t; // inexact_locate: c,i,z,pa
    typedef std::tuple<range_t,ulint,ulint,uint> loc_t;

    /*
    struct cmp_iel_score_lt {
        bool operator()(const iel_t& lhs, const iel_t& rhs) const {
            return std::get<3>(lhs).score < std::get<3>(rhs).score;
        }
    };
    struct cmp_iel_score_gt {
        bool operator()(const iel_t& lhs, const iel_t& rhs) const {
            return std::get<3>(lhs).score > std::get<3>(rhs).score;
        }
    };
    typedef std::priority_queue<iel_t, std::vector<iel_t>, cmp_iel_score_gt> iel_heap;
    */

    class ri_aligner {
        public:
            // ri-aligner(std::string input) {
            //     // create the forward and reverse indexes, save to file?
            // }
            // TODO: add constructor for initializing scoring penalties
            ri_aligner() {}
            // TODO: operator= not nyi for r_index
            // TODO: ri-aligner(r_index f, r_index r) fwd(f), rev(r) {}
            // load fwd and reverse indexes from files
            /*
            ri_aligner(std::string idx_fwd, std::string idx_rev) {
                // use ri-buildfasta to generate fwd,rev idxs
                bool fast;

                ifstream in1(idx_fwd);
                in1.read((char*)&fast,sizeof(fast));
                cerr << "loading from file " << idx_fwd << endl;
                fwd.load(in1);

                ifstream in2(idx_rev);
                in2.read((char*)&fast,sizeof(fast));
                cerr << "loading from file " << idx_rev << endl;
                rev.load(in2);
                rev_loaded = true;
            }
            */

            ri_aligner(std::string prefix) {
                bool fast;
                std::string fwd_fname = prefix + ".ri";
                std::string rev_fname = prefix + ".rev.ri";
                std::string seq_idx = prefix  + ".1.ri";
                if (!load_fwd(fwd_fname)) {
                    cerr << fwd_fname << " does not exist! Exiting...\n";
                    exit(1);
                }
                if (!load_rev(rev_fname)) {
                    // cerr << rev_fname << " does not exist. inexact matching will not be possible\n";
                }
                if (!load_seqidx(seq_idx)) {
                    cerr << seq_idx << " does not exist. output will default to offset into full BWT\n";
                    fai = false;
                } else fai = true;
            }

            bool load_fwd(std::string fname) {
                bool fast;
                ifstream in(fname);
                if (!in.good()) {
                   return false;
                }
                in.read((char*)&fast,sizeof(fast));
                // cerr << "loading fwd index from file " << fname << endl;
                fwd.load(in);
                return true;
            }


            bool load_rev(std::string fname) {
                bool fast;
                ifstream in(fname);
                if (!in.good()) {
                    return false;
                }
                in.read((char*)&fast,sizeof(fast));
                cerr << "loading reverse index from file " << fname << endl;
                rev.load(in);
                rev_loaded = true;
                return true;
            }

            bool load_seqidx(std::string fname) {
                ifstream in(fname);
                if (!in.good()) {
                    return false;
                }
                std::string seq_name;
                std::vector<uint64_t> seq_pos;
                uint64_t prev_pos = 0;
                uint64_t pos;
                int x = 1;
                while (in >> seq_name >> pos) {
                    seq_names.push_back(seq_name);
                    seq_pos.push_back(pos);
                    max_rank = x;
                    x += 1;
                }
                seq_bounds = std::move(sdsl::sd_vector<>(seq_pos.begin(), seq_pos.end()));
                seq_rank = sdsl::sd_vector<>::rank_1_type(&seq_bounds);
                seq_sel = sdsl::sd_vector<>::select_1_type(&seq_bounds);
                return true;
            }

            std::vector<uint> D_algo(const char* P, const uint m,  std::vector<uint>& D) {
                if (!rev_loaded) {
                    std::cerr << "reverse index not loaded, exiting..." << endl;
                    exit(1);
                }
                D.clear();
                D.resize(m);
                range_t range = rev.full_range();
                uint z = 0;
                ulint j = 0;
                ulint k = 0;
                uchar c;
                for (uint i = 0; i < m; i++) {
                    // do LF mapping against rev index
                    c = P[i];
                    range = rev.LF(range, c);
                    if (range.second < range.first) {
                        z++;
                        k,j = 0 ;
                        range = rev.full_range();
                    }
                    D[i] = z;
                }
                return D;
            }

            std::vector<uint> D_algo(const char* P, const uint m ) {
                std::vector<uint> D;
                return D_algo(P, m, D);
            }

            std::vector<uint> D_algo(std::string P) {
                std::vector<uint> D;
                return D_algo(P.c_str(), P.size(), D);
            }

            range_t exact_count(const std::string& P) {
                return fwd.count(P);
            }

            std::pair<range_t, ulint> exact_count_longest_suffix(const std::string& P) {
                auto range = fwd.full_range();
                ulint m = P.size();
                ulint i = 0;
                auto prev_range = range;
                for(i=0; i<m && range.second >= range.first; ++i) {
                    prev_range = range;
                    range = fwd.LF(range,P[m-i-1]);
                }
                if (range.second < range.first)
                    return std::make_pair(prev_range, i - 1);
                return std::make_pair(range, m);
            }

            // iec=inexact_count
            /*
            inline void iec_push_to_stack(const char c, const uint i, const uint z, const range_t cur_range, std::vector<iec_t>& stack) {
                range_t range;
                for (const char b: {'A', 'C', 'G', 'T'}) {
                    range = fwd.LF(cur_range, b);
                    if (range.second >= range.first) {
                        if (b == c) stack.push_back(iec_t(i-1, z, range));
                        else stack.push_back(iec_t(i-1, z-1, range));
                    }
                }
            }
            */

            /*
            std::vector<range_t> inexact_count(const char* P, const uint m, ri_opts_t opts, std::vector<range_t>& ranges) {
                ranges.clear();
                std::vector<uint> D = D_algo(P, m);
                range_t cur_range = fwd.full_range();
                std::vector<iec_t> stack; // i, z, range
                // push back initial things here
                iec_push_to_stack(P[m-1], m-1, opts.z, cur_range, stack);
                while (!stack.empty()) {
#if __cpp_structured_bindings
                    auto [i, z2, cur_range] = stack.back();
#else
                    auto params = stack.back();
                    i = std::get<0>(params);
                    z2 = std::get<1>(params);
                    cur_range = std::get<2>(params);
#endif
                    stack.pop_back();
                    if (i == (uint)-1 and z2 != (uint)-1) {
                        ranges.push_back(cur_range);
                    } else if ( !(z2 == (uint)-1 or z2 < D[i]) ) {
                        iec_push_to_stack(P[i], i, z2, cur_range, stack);
                    }
                }
                return ranges;
            }
            */

            // report locations of the string P within the reference.
            //  use a lambda function to capture ouput (vector<ulint>) 
            //  F(string, vector<ulint>)
            //  if no locations are found, then locs vector will be empty
            /*
            template<typename Func>
            void exact_locate(std::string P, const ri_opts_t& opts, Func f, vector<ulint>& locs) {
                locs.clear();
                auto range = fwd.locate_all(P, opts.max_hits, locs);
                f(P, range, locs);
            }
            */
            template<typename Func>
            void exact_locate(std::string P, const ri_opts_t& opts, Func f, vector<ulint>& locs) {
                // vector<ulint> locs;
                locs.clear();
                ulint k = 0;
                range_t range = fwd.full_range();
                k = fwd.get_last_run_sample();
                ulint m = P.size();
                for (uint i = 0; i < m; ++i) {
                    std::tie(range, k) = fwd.LF_w_loc(range, P[m-i-1], k);
                    if (range.second < range.first) {
                        return;
                    }
                }
                if (range.second-range.first+1 <= opts.max_range) {
                    fwd.locate_range(range.first, range.second, k, opts.max_hits, locs);
                }
                f(P,range,locs);
                return;
            }

            std::tuple<std::string, uint64_t> resolve_offset(uint64_t x) {
                if (!fai) return make_tuple("seq", x);
                size_t offset;
                size_t rank;
                if (x > seq_rank.size()) {
                    rank = max_rank;
                    offset = x - seq_sel(max_rank);
                }
                else {
                    rank = seq_rank(x);
                    offset = x - seq_sel(rank);
                }
                return std::tuple<std::string, uint64_t>(seq_names[rank-1], offset);
            }


            std::vector< std::pair<std::string, uint64_t> > get_ref_names_and_lengths() {
                std::vector<std::pair<std::string, uint64_t>> lengths;
                if (!fai) {
                    lengths.push_back(std::make_pair("seq", 0));
                    return lengths;
                }
                uint64_t prev_pos = 0;
                for (auto i = 0; i < seq_names.size() - 1; ++i) {
                    auto x = seq_sel(i + 2);
                    lengths.push_back(std::make_pair(seq_names[i], x - prev_pos)); 
                    prev_pos = x;
                }
                lengths.push_back(std::make_pair(seq_names[seq_names.size()-1], fwd.text_size() - prev_pos));
                return lengths;
            }

            // iel=inexact_locate
            /*
            template<typename RecurseType=vector<iel_t>>
            inline void iel_push(const char c, const uint i,
                                          const uint z, const partial_aln& pa,
                                          RecurseType& stack) {
                for (const char b: {'A', 'C', 'G', 'T'}) {
                    // create new pa
                    partial_aln npa(pa);
                    npa.pr = npa.r;
                    std::tie(npa.r, npa.k, npa.j) = fwd.LF_w_loc(npa.pr, b, npa.k, npa.j);
                    if (npa.r.second >= npa.r.first) {
                        if (b == c) {
                            stack.push_back(iel_t(b, i-1, z, npa));
                        }
                        else {
                            npa.mm += 1;
                            npa.score = score_fn(npa);
                            stack.push_back(iel_t(b, i-1, z-1, npa));
                        }
                    }
                }
            }

            template<typename RecurseType=vector<iel_t>>
            std::vector<ulint> inexact_locate(const char* P, const uint m, const ri_opts_t& opts, std::vector<ulint>& locs) {
                std::vector<uint> D = D_algo(P, m); // uses reverse index
                // vector<iel_t> stack; // TODO turn into priority queue
                RecurseType stack;
                iel_heap heap;
                uint i, z;
                ulint range_size;
                ulint num_locs = 0;
                uchar b;
                uint hits_thres = opts.max_hits;
                // push back initial things here
                partial_aln pa;
                pa.r = fwd.full_range();
                pa.pr = fwd.full_range();
                iel_push<RecurseType>(P[m-1], m-1, opts.z, pa, stack);
                while (!stack.empty()) {
                    // TODO: switch from stack to priotity queue (plus implement scoring fn)
                    std::tie(b, i, z, pa) = stack.back();
                    stack.pop_back();
                    if (i == (uint)-1 and z != (uint)-1) {
                        // sometimes toehold isn't captured if there are no runs in any LF calls
                        // so we correct that here (and k==j==0 is impossible, right?)
                        if (pa.k == 0 && pa.j == 0) {
                            std::tie(pa.k, pa.j) = fwd.get_sample_from_range(pa.pr, b);
                        }
                        // TODO: do locate queries later. Instead, push range onto priority queue
                        range_size = pa.r.second - pa.r.first + 1;
                        num_locs = locs.size();
                        if (hits_thres < range_size)
                            fwd.locate(pa.r, pa.k, pa.j, locs, hits_thres);
                        else
                            fwd.locate(pa.r, pa.k, pa.j, locs);
                        if (hits_thres >= (locs.size() - num_locs)) { hits_thres -= (locs.size() - num_locs);}
                        else {return locs;}
                    } else if ( !(z == (uint)-1 or z < D[i]) ) {
                        iel_push(P[i], i, z, pa, stack);
                    }
                }
                return locs;
            }
            */


            //  locate until failure, then restart,repeat from next position.
            //  use a lambda function to capture ouput for each "segment" (vector<ulint>):
            //  F(string input, uint segment_start, uint segment_end, vector<ulint> locs)
            template<typename Func>
            void piecewise_locate(std::string P, const ri_opts_t& opts, Func f) {
                ulint k = 0;
                range_t range = fwd.full_range();
                k = fwd.get_last_run_sample();
                ulint m = P.size();
                range_t range1;
                ulint k1;
                uint end = m-1;
                for (uint i = 0; i < m; ++i) {
                    std::tie(range1, k1) = fwd.LF_w_loc(range, P[m-i-1], k);
                    if (range1.second < range1.first) {
                        // locate range
                        vector<ulint> locs;
                        do_locate(range, k, opts, locs);
                        f(P, range, (m-i-1)+1, end ,locs);
                        // reset range, k, end
                        range = fwd.full_range();
                        k = fwd.get_last_run_sample(); 
                        end=m-i-1;
                    } else {
                        range = range1;
                        k = k1;
                    }
                }
                if (range.second >= range.first) { // capture final range if it hasn't already
                    vector<ulint> locs;
                    do_locate(range, k, opts, locs);
                    f(P,range,0, end,locs);
                }
                return;
            }

        private:
            ri::r_index<> fwd, rev;
            std::vector<std::string> seq_names;
            sdsl::sd_vector<> seq_bounds;
            sdsl::sd_vector<>::rank_1_type seq_rank;
            sdsl::sd_vector<>::select_1_type seq_sel;
            bool rev_loaded = false;
            bool fai = false;
            int max_rank = 0;
            uint mms = 1;
            uint gos = 0;
            uint ges = 0;
            inline uint score_fn(const partial_aln& pa) {
                return((mms)*(pa.mm) + (gos)*(pa.go) + (ges)*(pa.ge));
            }
            inline vector<ulint>& do_locate(range_t range, ulint k, ri_opts_t opts, vector<ulint>& locs) {
                if (range.second - range.first + 1 <= opts.max_range) {
                    fwd.locate_range(range.first, range.second, k, opts.max_hits, locs);
                }
                return locs;
            }
    };

}
#endif /*R_INDEX_S_F_H*/
