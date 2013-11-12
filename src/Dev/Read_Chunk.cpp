//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Read_Chunk.hpp"

#include <boost/tuple/tuple.hpp>

#include "Mutation.hpp"
#include "print_seq.hpp"
#include "indent.hpp"

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
     *
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
        vector< pair< Read_Chunk, Mutation_Cont > > res(2);

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
    */

    boost::tuple< Size_Type, Size_Type, size_t > Read_Chunk::get_cut_data(Size_Type brk, bool is_contig_brk) const
    {
        assert(not is_contig_brk or (get_c_start() <= brk and brk <= get_c_end()));
        assert(is_contig_brk or (get_r_start() <= brk and brk <= get_r_end()));

        Size_Type c_pos = get_c_start();
        Size_Type r_pos = (not _rc? get_r_start() : get_r_end());
        size_t i;

        // we use a fake Mutation to deal with end matching region smoothly
        Mutation fake_mut(get_c_end(), 0, 0);
        const Mutation* mut_cptr;
        for (i = 0; i <= _mut_ptr_cont.size(); ++i)
        {
            mut_cptr = (i < _mut_ptr_cont.size()? _mut_ptr_cont[i] : &fake_mut);

            // mutations are expected in contig order
            assert(c_pos <= mut_cptr->get_start());

            // progress
            assert(not is_contig_brk or c_pos <= brk);
            assert(is_contig_brk or _rc or r_pos <= brk);
            assert(is_contig_brk or not _rc or brk <= r_pos);

            Size_Type match_len = mut_cptr->get_start() - c_pos;
            // stop if breakpoint is inside the next matched region, or right after it
            if (is_contig_brk and brk <= c_pos + match_len)
            {
                match_len = brk - c_pos;
                c_pos = brk;
                r_pos = (not _rc? r_pos + match_len : r_pos - match_len);
                break;
            }
            else if (not is_contig_brk
                     and ((not _rc and brk <= r_pos + match_len)
                          or (_rc and r_pos - match_len <= brk)))
            {
                c_pos += (not _rc? brk - r_pos : r_pos - brk);
                r_pos = brk;
                break;
            }

            // advance past matched region
            c_pos += match_len;
            r_pos = (not _rc? r_pos + match_len : r_pos - match_len);

            // we do not get here if breakpoint is now on the edge of the last matching region
            assert(not is_contig_brk or c_pos < brk);
            assert(is_contig_brk or _rc or r_pos < brk);
            assert(is_contig_brk or not _rc or brk < r_pos);

            // stop if next mutation completely spans read breakpoint
            if ((is_contig_brk and brk < c_pos + mut_cptr->get_len())
                or (not is_contig_brk and ((not _rc and brk < r_pos + mut_cptr->get_seq_len())
                                           or (_rc and r_pos - mut_cptr->get_seq_len() < brk))))
            {
                break;
            }

            // otherwise, advance past next mutation
            c_pos += mut_cptr->get_len();
            r_pos = (not _rc? r_pos + mut_cptr->get_seq_len() : r_pos - mut_cptr->get_seq_len());
        }

        assert(get_c_start() <= c_pos and c_pos <= get_c_end());
        assert(get_r_start() <= r_pos and r_pos <= get_r_end());
        assert(not is_contig_brk or c_pos <= brk);
        assert(is_contig_brk or _rc or r_pos <= brk);
        assert(is_contig_brk or not _rc or brk <= r_pos);

        return boost::make_tuple< Size_Type, Size_Type, size_t >(c_pos, r_pos, i);
    }

    boost::tuple< const Mutation*, Size_Type, Size_Type > Read_Chunk::get_mutation_to_cut(
        Size_Type brk, bool is_contig_brk, const set< Size_Type >& brk_cont) const
    {
        Size_Type c_pos;
        Size_Type r_pos;
        size_t i;
        tie(c_pos, r_pos, i) = get_cut_data(brk, is_contig_brk);

        // if we reached brk, there is no mutation to be cut
        if ((is_contig_brk and c_pos == brk)
            or (not is_contig_brk and r_pos == brk))
            return boost::make_tuple< const Mutation*, Size_Type, Size_Type >(NULL, 0, 0);

        // if not, we must cut i-th mutation
        assert(i < _mut_ptr_cont.size());
        assert(c_pos == _mut_ptr_cont[i]->get_start());

        assert(not is_contig_brk or (c_pos < brk and brk < _mut_ptr_cont[i]->get_end()));
        assert(is_contig_brk or _rc or (r_pos < brk and brk < r_pos + _mut_ptr_cont[i]->get_seq_len()));
        assert(is_contig_brk or not _rc or (r_pos > brk and brk > r_pos - _mut_ptr_cont[i]->get_seq_len()));

        // range for possible contig cuts
        Size_Type range_start = (is_contig_brk?
                                 (not _rc? r_pos : r_pos - _mut_ptr_cont[i]->get_seq_len())
                                 : c_pos);
        Size_Type range_end = (is_contig_brk?
                               (not _rc? r_pos + _mut_ptr_cont[i]->get_seq_len() : r_pos)
                               : _mut_ptr_cont[i]->get_end());
        Size_Type sel_brk;

        // look for preferred breakpoint
        set< Size_Type >::iterator it = brk_cont.lower_bound(range_start);
        if (it != brk_cont.end() and *it <= range_end)
            sel_brk = *it;
        else
            sel_brk = range_start;

        if (is_contig_brk)
            return boost::make_tuple< const Mutation*, Size_Type, Size_Type >(_mut_ptr_cont[i], brk - c_pos,
                                                                       (not _rc? sel_brk - r_pos : r_pos - sel_brk));
        else
            return boost::make_tuple< const Mutation*, Size_Type, Size_Type >(_mut_ptr_cont[i], sel_brk - c_pos,
                                                                       (not _rc? brk - r_pos : r_pos - brk));
    }

    boost::tuple< bool, shared_ptr< Read_Chunk > > Read_Chunk::apply_contig_split(
        Size_Type c_brk, const map< const Mutation*, const Mutation* >& mut_cptr_map, const Contig_Entry* ce_cptr)
    {
        Size_Type c_pos = c_brk;
        Size_Type r_pos;
        size_t i;
        if (c_brk < get_c_start())
        {
            r_pos = (not _rc? get_r_start() : get_r_end());
            i = 0;
        }
        else if (get_c_end() < c_brk)
        {
            r_pos = (not _rc? get_r_end() : get_r_start());
            i = _mut_ptr_cont.size();
        }
        else
        {
            tie(c_pos, r_pos, i) = get_cut_data(c_brk, true);
            assert(c_pos == c_brk);
            // if there exists an insertion at c_brk that is not in the map, it will stay on the left side of the break
            if (i < _mut_ptr_cont.size()
                and _mut_ptr_cont[i]->get_start() == c_brk
                and _mut_ptr_cont[i]->is_ins()
                and mut_cptr_map.count(_mut_ptr_cont[i]) == 0)
            {
                r_pos = (not _rc? r_pos + _mut_ptr_cont[i]->get_seq_len() : r_pos - _mut_ptr_cont[i]->get_seq_len());
                ++i;
            }
            assert(i == _mut_ptr_cont.size() or mut_cptr_map.count(_mut_ptr_cont[i]) > 0);
        }

        if (r_pos == get_r_start() or r_pos == get_r_end())
        {
            // the chunk stays in one piece, but ce_ptr & mutation pointers might change
            bool move_to_rhs = false;
            assert(i == 0 or i == _mut_ptr_cont.size());
            if (c_brk <= _c_start)
            {
                move_to_rhs = true;
                // fix contig start coordinate
                _c_start -= c_brk;
                // fix mutation pointers
                vector< const Mutation* > tmp = _mut_ptr_cont;
                _mut_ptr_cont.clear();
                for (size_t j = 0; j < tmp.size(); ++j)
                {
                    assert(mut_cptr_map.count(tmp[j]) > 0);
                    _mut_ptr_cont.push_back(mut_cptr_map.find(tmp[j])->second);
                }
                // fix ce_ptr
                _ce_ptr = ce_cptr;
            }
            return boost::make_tuple< bool, shared_ptr< Read_Chunk > >(move_to_rhs, NULL);
        }
        else
        {
            // read chunk gets split in 2
            assert(get_r_start() < r_pos and r_pos < get_r_end());

            // create new read chunk for second part of the contig
            Read_Chunk* rc_ptr = new Read_Chunk();
            rc_ptr->_rc = _rc;
            rc_ptr->_re_ptr = _re_ptr;

            rc_ptr->_ce_ptr = ce_cptr;
            rc_ptr->_c_start = 0;
            rc_ptr->_c_len = _c_start + _c_len - c_brk;
            _c_len -= rc_ptr->_c_len;

            if (not _rc)
            {
                rc_ptr->_r_start = r_pos;
                rc_ptr->_r_len = _r_start + _r_len - r_pos;
                _r_len -= rc_ptr->_r_len;
            }
            else
            {
                rc_ptr->_r_start = _r_start;
                rc_ptr->_r_len = r_pos - _r_start;
                _r_start = r_pos;
                _r_len -= rc_ptr->_r_len;
            }

            // transfer mutations starting at i (to values under mapping) to read chunk for contig rhs
            for (size_t j = i; j < _mut_ptr_cont.size(); ++j)
            {
                assert(mut_cptr_map.count(_mut_ptr_cont[j]) > 0);
                rc_ptr->_mut_ptr_cont.push_back(mut_cptr_map.find(_mut_ptr_cont[j])->second);
            }

            // drop transfered mutations from this chunk
            _mut_ptr_cont.resize(i);

            return boost::make_tuple< bool, shared_ptr< Read_Chunk > >(false, shared_ptr< Read_Chunk >(rc_ptr));
        }
    }

    vector< boost::tuple< bool, Read_Chunk_CPtr, Size_Type, Size_Type > > Read_Chunk::collapse_mutations(
        const Read_Chunk& rc1, const Mutation_Extra_Cont& rc1_me_cont, const Read_Chunk& rc2)
    {
        assert(rc1.get_r_start() == rc2.get_c_start() and rc1.get_r_len() == rc2.get_c_len());

        vector< boost::tuple< bool, Read_Chunk_CPtr, Size_Type, Size_Type > > res;

        Size_Type c_pos = rc1.get_c_start();
        Size_Type r1_pos = (not rc1.get_rc()? rc1.get_r_start() : rc1.get_r_end());
        Size_Type r2_pos = (rc1.get_rc() == rc2.get_rc()? rc2.get_r_start() : rc2.get_r_end());
        size_t i1 = 0; // number of mutations completely passed on r1
        size_t k1 = 0; // length of prefix of mutation i1 that was already considered
        size_t i2 = 0; // number of mutations passed on r2 (counting from the end if rc1._rc)

        while (i1 < rc1._mut_ptr_cont.size() or i2 < rc2._mut_ptr_cont.size())
        {
            // iteration may not stop in the middle of a read mutation
            assert(not rc1.get_rc()? r1_pos <= rc1.get_r_end() : rc1.get_r_start() <= r1_pos);
            // i2 < n2 and r1 not reversed => before i2 mutation
            assert(not (i2 < rc2._mut_ptr_cont.size() and not rc1.get_rc()) or r1_pos <= rc2._mut_ptr_cont[i2]->get_start());
            // 0 <= i2 - 1 < n2 and r1 not reversed => after i2 - 1 mutation
            assert(not (0 < i2 and i2 < rc2._mut_ptr_cont.size() + 1 and not rc1.get_rc()) or rc2._mut_ptr_cont[i2 - 1]->get_end());
            // i2 < n2 and r1 reversed => after i2 mutation
            assert(not (i2 < rc2._mut_ptr_cont.size() and rc1.get_rc()) or rc2._mut_ptr_cont[i2]->get_end() <= r1_pos);
            // 0 <= i2 - 1 < n2 and r1 reversed => before i2 - 1 mutation
            assert(not (0 < i2 and i2 < rc2._mut_ptr_cont.size() + 1 and rc1.get_rc()) or r1_pos <= rc2._mut_ptr_cont[i2 - 1]->get_start());

            assert(c_pos <= rc1.get_c_end());
            // if k1 != 0, we must be in the middle of a contig mutation
            assert(not k1 != 0 or (i1 < rc1._mut_ptr_cont.size() and c_pos == rc1._mut_ptr_cont[i1]->get_start() + k1));
            // if k1 == 0 and i1 not at the end, the next contig mutation is at the right
            assert(not (k1 == 0 and i1 < rc1._mut_ptr_cont.size()) or c_pos <= rc1._mut_ptr_cont[i1]->get_start());

            // if we are not on the boundary of either a contig mutation or a read mutation, advance until next boundary
            if (k1 == 0
                and (i1 == rc1._mut_ptr_cont.size()
                     or c_pos == rc1._mut_ptr_cont[i1]->get_start())
                and (i2 == rc2._mut_ptr_cont.size()
                     or r1_pos == (not rc1.get_rc()? rc2._mut_ptr_cont[i2]->get_start() : rc2._mut_ptr_cont[i2]->get_end())))
            {
                Size_Type c_match_len = (i1 < rc1._mut_ptr_cont.size()? rc1._mut_ptr_cont[i1]->get_start() : rc1.get_c_end()) - c_pos;
                Size_Type r1_match_len = (not rc1.get_rc()?
                                          (i2 < rc2._mut_ptr_cont.size()? rc2._mut_ptr_cont[i2]->get_start() : rc1.get_r_end()) - r1_pos
                                          : r1_pos - (i2 < rc2._mut_ptr_cont.size()? rc2._mut_ptr_cont[i2]->get_end() : rc1.get_r_start()));
                Size_Type skip_len = min(c_match_len, r1_match_len);
                c_pos += skip_len;
                r1_pos = (not rc1.get_rc()? r1_pos + skip_len : r1_pos - skip_len);
                r2_pos = (rc1.get_rc() == rc2.get_rc()? r2_pos + skip_len : r2_pos - skip_len);
                continue;
            }

            // if we are on the boundary of a read mutation, consume all contig mutations that span it
            if (i2 < rc2._mut_ptr_cont.size()
                and r1_pos == (not rc1.get_rc()? rc2._mut_ptr_cont[i2]->get_start() : rc2._mut_ptr_cont[i2]->get_end()))
            {
                // advance past contig mutations until they cover the current read mutation
                Size_Type c_span = 0;
                Size_Type r1_span = 0;
                while (r1_span < rc2._mut_ptr_cont[i2]->get_len())
                {
                    if (c_pos + c_span < (i1 < rc1._mut_ptr_cont.size()? rc1._mut_ptr_cont[i1]->get_start() + k1 : rc1.get_c_end()))
                    {
                        // in the middle of a matched region
                        assert(k1 == 0);
                        Size_Type match_len = (i1 < rc1._mut_ptr_cont.size()? rc1._mut_ptr_cont[i1]->get_start() : rc1.get_c_end()) - (c_pos + c_span);
                        Size_Type skip_len = min(match_len, rc2._mut_ptr_cont[i2]->get_len() - r1_span);
                        c_span += skip_len;
                        r1_span += skip_len;
                    }
                    else
                    {
                        // on the border or inside a contig mutation
                        // could not have finished the contig mutations at this point
                        // because there are read bases to be consumed
                        assert(i1 < rc1._mut_ptr_cont.size());
                        assert(c_pos + c_span == rc1._mut_ptr_cont[i1]->get_start() + k1);

                        Size_Type r1_leftover = (rc1._mut_ptr_cont[i1]->get_seq_len() >= k1? rc1._mut_ptr_cont[i1]->get_seq_len() - k1 : 0);
                        Size_Type c_leftover = (rc1._mut_ptr_cont[i1]->get_len() >= k1? rc1._mut_ptr_cont[i1]->get_len() - k1 : 0);
                        Size_Type r1_skip_len = min(r1_leftover, rc2._mut_ptr_cont[i2]->get_len() - r1_span);
                        c_span += min(r1_skip_len, c_leftover);
                        r1_span += r1_skip_len;
                        k1 += r1_skip_len;
                        if (rc1._mut_ptr_cont[i1]->get_seq_len() <= k1)
                        {
                            // consumed entire read span of this mutation; advance past it
                            c_span += (k1 < rc1._mut_ptr_cont[i1]->get_len()? rc1._mut_ptr_cont[i1]->get_len() - k1 : 0);
                            k1 = 0;
                            ++i1;
                        }
                        else
                        {
                            // still in the middle of the same mutation
                            assert(r1_span == rc2._mut_ptr_cont[i2]->get_len());
                        }
                    }
                }
                // construct new mutation to replace this read mutation
                Mutation new_mut(c_pos, c_span, rc2._mut_ptr_cont[i2]->get_seq());
                c_pos += c_span;
                r1_pos = (not rc1.get_rc()? r1_pos + r1_span : r1_pos - r1_span);
                r2_pos = (rc1.get_rc() == rc2.get_rc()? r2_pos + rc2._mut_ptr_cont[i2]->get_seq_len() : r2_pos - rc2._mut_ptr_cont[i2]->get_seq_len());
                ++i2;
            }

            //TODO

            // consume deleted contig bases; if any
            if (i1 < rc1._mut_ptr_cont.size()
                and c_pos == rc1._mut_ptr_cont[i1]->get_start() + k1
                and rc1._mut_ptr_cont[i1]->get_seq_len() <= k1)
            {
                assert(rc1._mut_ptr_cont[i1]->get_len() > k1);
                res.push_back(true, _mut_ptr_cont[i1], k1, _mut_ptr_cont[i1]->get_len());
                c_pos += _mut_ptr_cont[i1]->get_len();
                k1 = 0;
                ++i1;
                continue;
            }

            

            // consider the mapped strech preceding the next contig mutation
            Size_Type match_len = (i1 < rc1._mut_ptr_cont.size()? rc1._mut_ptr_cont[i1]->get_start() : rc1.get_c_end()) - c_pos;
            // find the first read1 mutation
            
        }

        return res;
    }

    Mutation_Extra_Cont Read_Chunk::get_me_cont() const
    {
        Mutation_Extra_Cont res;

        Size_Type c_pos = _c_start;
        Size_Type r_pos = (not _rc? _r_start : _r_start + _r_len);
        for (size_t i = 0; i < _mut_ptr_cont.size(); ++i)
        {
            assert(c_pos <= _mut_ptr_cont[i]->get_start());
            Size_Type match_len = _mut_ptr_cont[i]->get_start() - c_pos;
            c_pos += match_len;
            r_pos = (not _rc? r_pos + match_len : r_pos - match_len);
            Mutation_Extra tmp;
            tmp.mut_cptr = _mut_ptr_cont[i];
            tmp.alt_start = r_pos;
            res.insert(tmp);
        }
        return res;
    }

    ostream& operator << (ostream& os, const Read_Chunk& rhs)
    {
        os << "(Read_Chunk &=" << (void*)&rhs
           << indent::inc << indent::nl << "re_cptr=" << (void*)rhs._re_ptr
           << ",ce_cptr=" << (void*)rhs._ce_ptr
           << indent::nl << "r_start=" << rhs.get_r_start()
           << ",r_len=" << rhs.get_r_len()
           << ",c_start=" << rhs.get_c_start()
           << ",c_len=" << rhs.get_c_len()
           << ",rc=" << (int)rhs.get_rc()
           << indent::nl << "mut_cptr_cont="
           << indent::inc;
        print_seq(os, rhs._mut_ptr_cont, indent::nl, indent::nl);
        os << indent::dec << indent::dec << indent::nl << ")";
        return os;
    }
}
