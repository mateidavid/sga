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
#include "../Util/Util.h"

using namespace std;
using boost::tie;


namespace MAC
{
    ostream& operator << (ostream& os, const Read_Chunk_Pos& pos)
    {
        os << "(c_pos=" << (size_t)pos.c_pos << ",r_pos=" << (size_t)pos.r_pos << ",mut_idx=" << pos.mut_idx << ",mut_offset=" << pos.mut_offset << ")";
        return os;
    }

    inline bool Read_Chunk::check_pos(const Read_Chunk_Pos& pos) const
    {
        assert(get_c_start() <= pos.c_pos and pos.c_pos <= get_c_end());
        assert(get_r_start() <= pos.r_pos and pos.r_pos <= get_r_end());
        assert(pos.mut_idx <= _mut_ptr_cont.size());
        assert(not pos.mut_offset == 0
               or pos.c_pos <= (pos.mut_idx < _mut_ptr_cont.size()? _mut_ptr_cont[pos.mut_idx]->get_start() : get_c_end()));
        assert(not pos.mut_offset != 0
               or (pos.mut_idx < _mut_ptr_cont.size()
                   and pos.c_pos == _mut_ptr_cont[pos.mut_idx]->get_start() + min(pos.mut_offset, _mut_ptr_cont[pos.mut_idx]->get_len())));
        return true;
    }

    inline Size_Type Read_Chunk::get_match_len_from_pos(const Read_Chunk_Pos& pos, bool forward) const
    {
        assert(check_pos(pos));
        if (pos.mut_offset != 0)
        {
            return 0;
        }
        if (forward)
        {
            return (pos.mut_idx < _mut_ptr_cont.size()? _mut_ptr_cont[pos.mut_idx]->get_start() : get_c_end()) - pos.c_pos;
        }
        else
        {
            return pos.c_pos - (pos.mut_idx > 0? _mut_ptr_cont[pos.mut_idx - 1]->get_end() : get_c_start());
        }
    }

    void Read_Chunk::increment_pos(Read_Chunk_Pos& pos, Size_Type brk, bool on_contig) const
    {
        assert(check_pos(pos));
        if (brk == 0)
        {
            // by default, use an unreachable breakpoint
            on_contig = true;
            brk = get_c_end() + 1;
        }
        assert(on_contig? pos.c_pos < brk : (not _rc? pos.r_pos < brk : brk < pos.r_pos));
        if (pos.mut_idx < _mut_ptr_cont.size()
            and pos.c_pos == _mut_ptr_cont[pos.mut_idx]->get_start() + min(pos.mut_offset, _mut_ptr_cont[pos.mut_idx]->get_len()))
        {
            // at the start or inside a mutation
            Size_Type c_leftover = (_mut_ptr_cont[pos.mut_idx]->get_len() > pos.mut_offset? _mut_ptr_cont[pos.mut_idx]->get_len() - pos.mut_offset : 0);
            Size_Type r_leftover = (_mut_ptr_cont[pos.mut_idx]->get_seq_len() > pos.mut_offset? _mut_ptr_cont[pos.mut_idx]->get_seq_len() - pos.mut_offset : 0);
            Size_Type skip_len = 0;
            if (on_contig and brk < pos.c_pos + c_leftover)
            {
                skip_len = brk - pos.c_pos;
            }
            else if (not on_contig and (not _rc? brk < pos.r_pos + r_leftover : pos.r_pos - r_leftover < brk))
            {
                skip_len = (not _rc? brk - pos.r_pos : pos.r_pos - brk);
            }
            else
            {
                // breakpoint does not occur later in this mutation: we advance past it entirely
                skip_len = max(c_leftover, r_leftover);
            }
            pos.c_pos += min(c_leftover, skip_len);
            pos.r_pos = (not _rc? pos.r_pos + min(r_leftover, skip_len) : pos.r_pos - min(r_leftover, skip_len));
            pos.mut_offset += skip_len;
            if (pos.mut_offset >= _mut_ptr_cont[pos.mut_idx]->get_len() and pos.mut_offset >= _mut_ptr_cont[pos.mut_idx]->get_seq_len())
            {
                pos.mut_offset = 0;
                ++pos.mut_idx;
            }
        }
        else
        {
            // a match stretch follows
            Size_Type match_len = get_match_len_from_pos(pos, true);
            assert(match_len > 0);
            Size_Type skip_len = 0;
            if (on_contig and brk < pos.c_pos + match_len)
            {
                skip_len = brk - pos.c_pos;
            }
            else if (not on_contig and (not _rc? brk < pos.r_pos + match_len : pos.r_pos - match_len < brk))
            {
                skip_len = (not _rc? brk - pos.r_pos : pos.r_pos - brk);
            }
            else
            {
                skip_len = match_len;
            }
            pos.c_pos += skip_len;
            pos.r_pos = (not _rc? pos.r_pos + skip_len : pos.r_pos - skip_len);
        }
    }

    void Read_Chunk::decrement_pos(Read_Chunk_Pos& pos, Size_Type brk, bool on_contig) const
    {
        assert(check_pos(pos));
        if (brk == 0)
        {
            on_contig = true;
        }
        assert(on_contig? brk <= pos.c_pos : (not _rc? brk < pos.r_pos : pos.r_pos < brk));
        if (pos.mut_offset != 0
            or (pos.mut_offset == 0 and pos.mut_idx > 0
                and pos.c_pos == _mut_ptr_cont[pos.mut_idx - 1]->get_end()))
        {
            // at the end or inside a mutation
            if (pos.mut_offset == 0)
            {
                --pos.mut_idx;
                pos.mut_offset = max(_mut_ptr_cont[pos.mut_idx]->get_len(), _mut_ptr_cont[pos.mut_idx]->get_seq_len());
            }
            assert(pos.mut_offset > 0);
            Size_Type c_leftover = min(_mut_ptr_cont[pos.mut_idx]->get_len(), pos.mut_offset);
            Size_Type r_leftover = min(_mut_ptr_cont[pos.mut_idx]->get_seq_len(), pos.mut_offset);
            Size_Type c_skip_len = 0;
            Size_Type r_skip_len = 0;
            if (on_contig and pos.c_pos - c_leftover < brk)
            {
                c_skip_len = pos.c_pos - brk;
                r_skip_len = (r_leftover > c_leftover - c_skip_len? r_leftover - (c_leftover - c_skip_len) : 0);
            }
            else if (not on_contig and (not _rc? pos.r_pos - r_leftover < brk : brk < pos.r_pos + r_leftover))
            {
                r_skip_len = (not _rc? pos.r_pos - brk : brk - pos.r_pos);
                c_skip_len = (c_leftover > r_leftover - r_skip_len? c_leftover - (r_leftover - r_skip_len) : 0);
            }
            else
            {
                // breakpoint does not occur earlier in this mutation: we advance past it entirely
                c_skip_len = c_leftover;
                r_skip_len = r_leftover;
            }
            pos.c_pos -= c_skip_len;
            pos.r_pos = (not _rc? pos.r_pos - r_skip_len : pos.r_pos + r_skip_len);
            pos.mut_offset = max(c_leftover - c_skip_len, r_leftover - r_skip_len);
        }
        else
        {
            // a match stretch follows
            Size_Type match_len = get_match_len_from_pos(pos, false);
            assert(match_len > 0);
            Size_Type skip_len = 0;
            if (on_contig and pos.c_pos - match_len < brk)
            {
                skip_len = pos.c_pos - brk;
            }
            else if (not on_contig and (not _rc? pos.r_pos - match_len < brk : brk < pos.r_pos + match_len))
            {
                skip_len = (not _rc? pos.r_pos - brk : brk - pos.r_pos);
            }
            else
            {
                skip_len = match_len;
            }
            pos.c_pos -= skip_len;
            pos.r_pos = (not _rc? pos.r_pos - skip_len : pos.r_pos + skip_len);
        }
    }

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
     * @param cigar Cigar object describing the match.
     * @param qr Query string; optional: if not given, Mutation objects contain alternate string lengths only.
     * @return Container of Mutation objects.
     */
    shared_ptr< Mutation_Cont > make_mutations_from_cigar(const Cigar& cigar, const string& qr)
    {
        shared_ptr< Mutation_Cont > res(new Mutation_Cont());
        for (size_t i = 0; i < cigar.get_n_ops(); ++i)
        {
            // disallow ambigous 'M' operation
            assert(cigar.get_op(i) != 'M');
            if (cigar.get_op(i) != '=')
            {
                Mutation mut;
                Mutation_Cont::iterator it;
                bool success;
                if (qr.size() > 0)
                {
                    // accept either the matched part of the query (length = cigar.qr_len)
                    // or entire query (length > cigar.qr_len)
                    assert(qr.size() >= cigar.get_qr_len());
                    Size_Type qr_offset = (not cigar.is_reversed()? cigar.get_qr_offset(i) : cigar.get_qr_offset(i) - cigar.get_qr_op_len(i));
                    if (qr.size() == cigar.get_qr_len())
                    {
                        // assume the query sequence before qr_start is missing
                        qr_offset -= cigar.get_qr_start();
                    }
                    mut = Mutation(cigar.get_rf_offset(i), cigar.get_rf_op_len(i),
                                   (not cigar.is_reversed()?
                                    qr.substr(qr_offset, cigar.get_qr_op_len(i))
                                    : reverseComplement(qr.substr(qr_offset, cigar.get_qr_op_len(i)))));
                }
                else
                {
                    mut = Mutation(cigar.get_rf_offset(i), cigar.get_rf_op_len(i), cigar.get_qr_op_len(i));
                }
                tie(it, success) = res->insert(mut);
                assert(success);
            }
        }
        return res;
    }

    boost::tuple< Read_Chunk, shared_ptr< Mutation_Cont > > Read_Chunk::make_chunk_from_cigar(const Cigar& cigar, const string& qr)
    {
        // create objects with default constructor
        boost::tuple< Read_Chunk, shared_ptr< Mutation_Cont > > res;

        // fix lengths and rc flags
        res.get<0>()._c_start = cigar.get_rf_start();
        res.get<0>()._c_len = cigar.get_rf_len();
        res.get<0>()._r_start = cigar.get_qr_start();
        res.get<0>()._r_len = cigar.get_qr_len();
        res.get<0>()._rc = cigar.is_reversed();

        // construct mutations
        res.get<1>() = make_mutations_from_cigar(cigar, qr);

        // store pointers in the Read_Chunk object
        for (auto it = res.get<1>()->begin(); it != res.get<1>()->end(); ++it)
        {
            res.get<0>()._mut_ptr_cont.push_back(&*it);
        }

        return res;
    }

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

    boost::tuple< shared_ptr< map< Mutation_CPtr, Mutation_CPtr > >, shared_ptr< Mutation_Cont > >
    Read_Chunk::lift_read_mutations(const Mutation_Cont& mut_cont) const
    {
        shared_ptr< map< Mutation_CPtr, Mutation_CPtr > > mut_map_sptr(new map< Mutation_CPtr, Mutation_CPtr >());
        shared_ptr< Mutation_Cont > new_mut_cont_sptr(new Mutation_Cont());

        Read_Chunk_Pos pos = (not _rc? get_start_pos() : get_end_pos());
        Read_Chunk_Pos pos_end = (not _rc? get_end_pos() : get_start_pos());
        for (auto old_mut_it = mut_cont.begin(); old_mut_it != mut_cont.end(); ++old_mut_it)
        {
            const Mutation& old_mut = *old_mut_it;
            Mutation m;
            Size_Type c_start;
            Size_Type c_end;

            // look for the start of read mutation
            while (pos != pos_end and pos.r_pos != old_mut.get_start())
            {
                advance_pos(pos, not _rc, old_mut.get_start(), false);
            }
            assert(pos.r_pos == old_mut.get_start());
            // NOTE: Here we could detect a contig deletion perfectly matched with a read insertion
            // in order to enforce a certain order.
            // By not doing anything special, the order is: if not _rc, DI; if _rc, ID.

            // Continue advancing past any deletions.
            Read_Chunk_Pos next_pos = pos;
            while (next_pos.r_pos == old_mut.get_start())
            {
                pos = next_pos;
                if (next_pos == pos_end)
                    break;
                advance_pos(next_pos, not _rc);
            }
            assert(pos.r_pos == old_mut.get_start());

            if (not _rc)
                c_start = pos.c_pos;
            else
                c_end = pos.c_pos;

            // look for end of read mutation
            next_pos = pos;
            while (next_pos != pos_end and next_pos.r_pos != old_mut.get_end())
            {
                advance_pos(next_pos, not _rc, old_mut.get_end(), false);
            }
            assert(pos.r_pos == old_mut.get_end());

            if (not _rc)
                c_end = next_pos.c_pos;
            else
                c_start = next_pos.c_pos;

            if (old_mut.have_seq())
                m = Mutation(c_start, c_end, not _rc? old_mut.get_seq() : reverseComplement(old_mut.get_seq()));
            else
                m = Mutation(c_start, c_end, old_mut.get_seq_len());

            Mutation_Cont::iterator new_mut_it;
            bool success;
            tie(new_mut_it, success) = new_mut_cont_sptr->insert(m);
            assert(success);
            auto it = mut_map_sptr->begin();
            tie(it, success) = mut_map_sptr->insert(pair< Mutation_CPtr, Mutation_CPtr >(&*old_mut_it, &*new_mut_it));
            assert(success);
        }

        return boost::make_tuple(mut_map_sptr, new_mut_cont_sptr);
    }

    /*
    shared_ptr< vector< boost::tuple< Mutation_CPtr, Size_Type, Size_Type > > >
    Read_Chunk::get_old_mutations_in_mapping(const vector< Mutation_CPtr >& new_mut_cptr_cont, const map< Mutation_CPtr, Mutation_CPtr >& mut_map) const
    {
        shared_ptr< vector< boost::tuple< Mutation_CPtr, Size_Type, Size_Type > > > res(new vector< boost::tuple< Mutation_CPtr, Size_Type, Size_Type > >());

        Read_Chunk_Pos pos = get_start_pos();
        size_t old_idx = 0; // next old mutation to look for
        size_t new_idx = 0; // next new mutation to look for
        while (old_idx < _mut_ptr_cont.size())
        {
            const Mutation& old_mut = *_mut_ptr_cont[old_idx];
            Size_Type brk = old_mut.get_start();
            if (new_idx < new_mut_cptr_cont.size())
            {
                brk = min(brk, new_mut_cptr_cont[new_idx]->get_start());
            }
            
        }

        return res;
    }
    */

    /*
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
    */

    vector< Size_Type > Read_Chunk::get_mut_pos() const
    {
        vector< Size_Type > res;
        for (Read_Chunk_Pos pos = get_start_pos(); not (pos == get_end_pos()); increment_pos(pos))
        {
            // if no breakpoints are used, we should never stop in the middle of a mutation
            assert(pos.mut_offset == 0);
            if (get_match_len_from_pos(pos) == 0)
            {
                // at the start of a mutation
                res.push_back(pos.r_pos);
            }
        }
        assert(res.size() == _mut_ptr_cont.size());
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
