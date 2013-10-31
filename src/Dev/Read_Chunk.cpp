//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Read_Chunk.hpp"

#include <boost/tuple/tuple.hpp>

#include "Mutation.hpp"
#include "print_seq.hpp"

using namespace std;
using boost::tie;


namespace MAC
{
    bool Read_Chunk::have_mutation(const Mutation* mut_cptr) const
    {
        for (size_t i = 0; i < _mut_ptr_cont.size(); ++i)
        {
            if (_mut_ptr_cont[i] == mut_cptr)
                return true;
        }
        return false;
    }

    void Read_Chunk::cond_add_mutation(const Mutation* m_old_cptr, const Mutation* m_new_cptr)
    {
        size_t i;
        for (i = 0; i < _mut_ptr_cont.size(); ++i)
        {
            if (_mut_ptr_cont[i] == m_old_cptr)
            {
                _mut_ptr_cont.insert(_mut_ptr_cont.begin() + i + 1, m_new_cptr);
                break;
            }
        }
    }

    /** Create a set of mutations to a reference string based on a cigar object.
     * @param rf_start Index in reference of the start of the match (0-based).
     * @param cigar Cigar object describing the match.
     * @return Container of Mutation objects.
     */
    Mutation_Cont make_mutations_from_cigar(Size_Type rf_start, const Cigar& cigar)
    {
        Mutation_Cont res;
        for (size_t i = 0; i < cigar.get_n_ops(); ++i)
        {
            if (not Cigar::is_match_op(cigar.get_op(i)))
            {
                Mutation_Cont::iterator it;
                bool success;
                tie(it, success) = res.insert(Mutation(rf_start + cigar.get_rf_offset(i), cigar.get_rf_op_len(i), cigar.get_qr_op_len(i)));
                assert(success);
            }
        }
        return res;
    }

    vector< pair< Read_Chunk, Mutation_Cont > > Read_Chunk::make_chunks_from_cigar(
        Size_Type r1_start, Size_Type r2_start, const Cigar& cigar)
    {
        // create objects with default constructor
        vector< pair< Read_Chunk, Mutation_Cont > > res;

        // fix lengths and rc flags
        Size_Type r1_len = cigar.get_rf_len();
        Size_Type r2_len = cigar.get_qr_len();
        res[0].first._r_start = r2_start;
        res[0].first._r_len = r2_len;
        res[0].first._c_start = r1_start;
        res[0].first._c_len = r1_len;
        res[0].first._rc = cigar.is_reversed();
        res[1].first._r_start = r1_start;
        res[1].first._r_len = r1_len;
        res[1].first._c_start = r2_start;
        res[1].first._c_len = r2_len;
        res[1].first._rc = cigar.is_reversed();

        // construct r1->r2 mutations
        res[0].second = make_mutations_from_cigar(r1_start, cigar);
        for (auto it = res[0].second.begin(); it != res[0].second.end(); ++it)
            res[0].first._mut_ptr_cont.push_back(&(*it));

        // construct r2->r1 mutations
        res[1].second = make_mutations_from_cigar(r2_start, cigar.complement());
        for (auto it = res[1].second.begin(); it != res[1].second.end(); ++it)
            res[1].first._mut_ptr_cont.push_back(&(*it));

        return res;
    }

    boost::tuple< Size_Type, Size_Type, size_t > Read_Chunk::get_read_split_data(Size_Type r_brk) const
    {
        size_t r_idx;
        size_t c_idx;
        size_t n = _mut_ptr_cont.size();
        Size_Type r_pos = get_r_start();
        Size_Type c_pos = (not _rc? get_c_start() : get_c_end());

        assert(get_r_start() <= r_brk and r_brk <= get_r_end());

        // invariant:
        //  - [_r_start,r_pos) is matched to [_c_start,c_pos) or [c_pos,_c_end) if rc
        //  - for remaining mutations, r_start >= r_pos
        //
        // Hack:
        //   To deal with the last matched portion consistently, we "imagine" an (n+1) empty read mutation.
        //   If not rc, this has (start = end = _c_end; len = 0; seq_len = 0).
        //   If rc, this has (start = end = _c_start; len = 0; seq_len = 0).
        Mutation fake_mut((not _rc? get_c_end() : get_c_start()), 0, 0);
        const Mutation* mut_cptr;
        for (r_idx = 0; r_idx <= n; ++r_idx)
        {
            // r_idx-th mutation in read == c_idx-th mutation in contig
            // Note: this is invalid when r_idx == n; must protect access
            c_idx = (not _rc ? r_idx : n - 1 - r_idx);
            mut_cptr = (r_idx < n? _mut_ptr_cont[c_idx] : &fake_mut);

            // before we consider r_idx-th read mutation,
            // we must account for the matched (non-mutated) region preceding it
            // read may not observe overlapping or out-of-range mutations in the reference
            Size_Type match_len;
            if (not _rc)
            {
                assert(c_pos <= mut_cptr->get_start());
                match_len = mut_cptr->get_start() - c_pos;
            }
            else
            {
                assert(mut_cptr->get_end() <= c_pos);
                match_len = c_pos - mut_cptr->get_end();
            }

            // compute read start&end of next mutation
            Size_Type m_r_start = r_pos + match_len;
            Size_Type m_r_end = m_r_start + mut_cptr->get_seq_len();

            if (m_r_start >= r_brk)
            {
                // break occurs inside matched portion, or immediately after it
                Size_Type r_leftover = r_brk - r_pos;
                r_pos += r_leftover;
                c_pos = (not _rc? c_pos + r_leftover : c_pos - r_leftover);
                break;
            }
            else // m_r_start < r_brk
            {
                // definitely go past matched portion
                r_pos += match_len;
                c_pos = (not _rc? c_pos + match_len : c_pos - match_len);

                // now consider advancing past r_idx-th read mutation
                if (m_r_end > r_brk)
                {
                    // we found a mutation with m_r_start < r_brk < m_r_end
                    break;
                }
                else // m_r_end <= r_brk
                {
                    // go past it as well
                    r_pos += mut_cptr->get_seq_len();
                    c_pos = (not _rc? c_pos + mut_cptr->get_len() : c_pos - mut_cptr->get_len());
                }
            }
        }

        // loop should exit with m_r_start >= r_brk before the last increment
        assert(r_idx < n + 1);

        // basic bounds
        assert(get_r_start() <= r_pos and r_pos <= get_r_end());
        assert(get_c_start() <= c_pos and c_pos <= get_c_end());

        // either r_pos < r_brk or r_pos == r_brk
        assert(r_pos <= r_brk);

        // final adjustment: we found the first mutation in read order after the break;
        // if break is not inside mutation, we must return the first mutation after the break in contig order
        if (_rc and r_pos == r_brk)
        {
            if (r_idx == n or _mut_ptr_cont[c_idx]->get_end() < c_pos)
            {
                --r_idx;
                c_idx = n - 1 - r_idx;
            }
            else
            {
                while (c_idx > 0 and _mut_ptr_cont[c_idx]->get_end() == c_pos)
                {
                    ++r_idx;
                    --c_idx;
                }
            }
        }

        // if r_pos < r_brk:
        //   r_idx < n, and
        //   r_brk < r_pos + m_seq_len, and
        //   c_pos == m_start (no _rc)
        //   c_pos == m_end (rc)
        assert(not (r_pos < r_brk) or r_idx < n);
        assert(not (r_pos < r_brk) or r_brk < r_pos + _mut_ptr_cont[c_idx]->get_seq_len());
        assert(not (r_pos < r_brk) or _rc or c_pos == _mut_ptr_cont[c_idx]->get_start());
        assert(not (r_pos < r_brk) or not _rc or _mut_ptr_cont[c_idx]->get_end() == c_pos);

        // if r_pos == r_brk:
        //   if r_idx < n:
        //     c_pos <= m_start (no _rc)
        //     m_end <= c_pos (rc)
        assert(not (r_pos == r_brk) or not (r_idx < n) or _rc or c_pos <= _mut_ptr_cont[c_idx]->get_start());
        assert(not (r_pos == r_brk) or not (r_idx < n) or not _rc or _mut_ptr_cont[c_idx]->get_end() <= c_pos);

        return boost::make_tuple(r_pos, c_pos, c_idx);
    }

    boost::tuple< const Mutation*, Size_Type, Size_Type > Read_Chunk::find_mutation_to_cut(Size_Type r_brk, const set< Size_Type >& c_brk_cont) const
    {
        Size_Type r_pos;
        Size_Type c_pos;
        size_t i;
        boost::tie(r_pos, c_pos, i) = get_read_split_data(r_brk);

        // if we reached r_pos, there is no mutation to be cut
        if (r_pos == r_brk)
            return boost::make_tuple< const Mutation*, Size_Type, Size_Type >(NULL, 0, 0);

        // if not, we must cut i-th mutation
        assert(i < _mut_ptr_cont.size());
        assert(r_pos < r_brk and r_brk < r_pos + _mut_ptr_cont[i]->get_seq_len());
        assert(_rc or c_pos == _mut_ptr_cont[i]->get_start());
        assert(not _rc or c_pos == _mut_ptr_cont[i]->get_end());

        Size_Type c_brk;
        // range for possible contig cuts is [m_c_start...m_c_end]
        // look for preferred contig breakpoint
        set< Size_Type >::iterator it = c_brk_cont.lower_bound(_mut_ptr_cont[i]->get_start());
        if (it != c_brk_cont.end() and *it <= _mut_ptr_cont[i]->get_end())
            c_brk = *it;
        else
            c_brk = c_pos;

        Size_Type r_offset = (not _rc? r_brk - r_pos : r_pos + _mut_ptr_cont[i]->get_seq_len() - r_brk);
        return boost::make_tuple< const Mutation*, Size_Type, Size_Type >(_mut_ptr_cont[i], c_brk - _mut_ptr_cont[i]->get_start(), r_offset);
    }

    /*
    pair< Size_Type, Size_Type > Read_Chunk::get_contig_brk_range(Size_Type pos)
    {
        pair< Size_Type, Size_Type > res;

        assert(_r_start <= pos and pos <= _r_start + _r_len);

        boost::tuple< size_t, Size_Type, Size_Type > t = count_mutations_left_of_brk(pos);

        size_t j = (not _rc? t.get<0>() : _mut_ptr_cont.size() - 1 - t.get<0>());
        Size_Type r_pos = _r_start + t.get<2>();
        Size_Type c_pos;
        if (not _rc)
        {
            c_pos = _c_start + t.get<1>();
        }
        else
        {
            c_pos = _c_start + _c_len - t.get<1>();
        }

        if (r_pos < pos)
        {
            // read breakpoint is strictly inside of a mutation with non-empty alternate sequence (insertion/multisnp)
            assert(pos < r_pos + _mut_ptr_cont[j]->get_seq_len());
            assert(rc or c_pos == _mut_ptr_cont[j]->get_start());
            assert(not _rc or _mut_ptr_cont->get_end() == c_pos);
            res.first = c_pos;
            res.second = _mut_ptr_cont[j]->get_end();
        }
        else
        {
            assert(r_pos == pos);
            //TODO:
        }
        
        assert(_r_start + t.get<2>() <= pos);
        if (t.get<0>() == _mut_ptr_cont.size())
        {
            // all mutations are left of breakpoint
            Size_Type r_leftover = pos - (_r_start + t.get<2>());
            if (not _rc)
            {
                res.first = _c_start + t.get<1>() + r_leftover;
                res.second = res.first;
            }
            else
            {
                res.first = _c_start + _c_len - t.get<1>() - r_leftover;
                res.second = res.first;
            }
        }
        else
        {
            // there are mutations at or right of read breakpoint
            size_t j = (not _rc? t.get<0>() : _mut_ptr_cont.size() - 1 - t.get<0>());
            Size_Type c_match;
            if (not _rc)
            {
                assert(_c_start + t.get<1>() <= _mut_ptr_cont[j]->get_start());
                c_match = _mut_ptr_cont[j]->get_start() - (_c_start + t.get<1>());
            }
            else
            {
                assert(_mut_ptr_cont[j]->get_end() <= _c_start + _c_len - t.get<1>());
                c_match = (_c_start + _c_len - t.get<1>()) - _mut_ptr_cont[j]->get_end();
            }

            // does the breakpoint occur in the matched portion before the mutation?
            if (_r_start + t.get<2>() + c_match < pos)
            {
                // yes
                Size_Type r_leftover = pos - (_r_start + t.get<2>() + c_match);
                if (not _rc)
                {
                    res.first = _c_start + t.get<1>() + r_leftover;
                    res.second = res.first;
                }
                else
                {
                    res.first = _c_start + _c_len - t.get<1>() - r_leftover;
                    res.second = res.first;
                }
            }
            else
            {
                // no; breakpoint exactly after matched portion
                assert(_r_start + t.get<2>() + c_match == pos);
                
            }
        }
    }
    */

    ostream& operator << (ostream& os, const Read_Chunk& rhs)
    {
        os << "(r_start=" << rhs.get_r_start() << ",r_len=" << rhs.get_r_len()
        << ",c_start=" << rhs.get_c_start() << ",c_len=" << rhs.get_c_len() << ",rc=" << (int)rhs.get_rc() << ",mut_list=(";
        print_ptr_seq(os, rhs._mut_ptr_cont);
        os << "))";
        return os;
    }
}
