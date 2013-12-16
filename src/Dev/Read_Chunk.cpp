//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Read_Chunk.hpp"

#include "Mutation.hpp"
#include "Contig_Entry.hpp"
#include "print_seq.hpp"
#include "indent.hpp"
#include "../Util/Util.h"

using namespace std;


namespace MAC
{
    ostream& operator << (ostream& os, const Read_Chunk_Pos& pos)
    {
        os << "(c_pos=" << (size_t)pos.c_pos << ",r_pos=" << (size_t)pos.r_pos << ",mut_idx=" << pos.mut_idx << ",mut_offset=" << pos.mut_offset << ")";
        return os;
    }

    Seq_Type Read_Chunk::get_seq() const
    {
        Seq_Type res;
        Read_Chunk_Pos pos = (not _rc? get_start_pos() : get_end_pos());
        while (pos != (not _rc? get_end_pos() : get_start_pos()))
        {
            Read_Chunk_Pos pos_next = pos;
            advance_pos(pos_next, not _rc);
            Size_Type match_len = get_match_len_from_pos(pos, not _rc);
            if (match_len > 0)
            {
                string tmp = _ce_ptr->substr((not _rc? pos.c_pos : pos_next.c_pos) - _ce_ptr->get_seq_offset(), match_len);
                res += (not _rc? tmp : reverseComplement(tmp));
            }
            else if (pos.r_pos != pos_next.r_pos)
            {
                if (not _rc)
                {
                    res += _mut_ptr_cont[pos.mut_idx]->get_seq().substr(pos.mut_offset, pos_next.r_pos - pos.r_pos);
                }
                else
                {
                    res += reverseComplement(_mut_ptr_cont[pos_next.mut_idx]->get_seq().substr(pos_next.mut_offset, pos.r_pos - pos_next.r_pos));
                }
            }
            pos = pos_next;
        }
        return res;
    }

    Seq_Type Read_Chunk::substr(Size_Type start, Size_Type len) const
    {
        assert(start >= _r_start and start + len <= _r_start + _r_len);
        return get_seq().substr(start -_r_start, len);
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

    bool Read_Chunk::advance_pos_til_mut(Read_Chunk_Pos& pos, const Mutation& mut, bool forward) const
    {
        Size_Type mut_first;
        Size_Type mut_second;
        if (forward == not _rc)
        {
            mut_first = mut.get_start();
            mut_second = mut.get_end();
            assert(pos.r_pos <= mut_first);
        }
        else
        {
            mut_first = mut.get_end();
            mut_second = mut.get_start();
            assert(mut_first <= pos.r_pos);
        }
        if (pos.r_pos != mut_first)
        {
            // there are still mapping stretches before the read Mutation
            // advance using mut_first as read breakpoint
            advance_pos(pos, forward, mut_first, false);
            return false;
        }
        else // pos.r_pos == mut_first
        {
            // we produce any remaining contig deletions at current read position,
            // except for the case when a contig deletion matches a read insertion;
            // in that case we produce the deletion first only when moving forward
            if ((forward or not mut.is_ins())
                and ((forward and pos != get_end_pos()) or (not forward and pos != get_start_pos())))
            {
                Read_Chunk_Pos pos_next = pos;
                advance_pos(pos_next, forward);
                if (pos_next.r_pos == pos.r_pos)
                {
                    pos = pos_next;
                    return false;
                }
            }
            // at this point, we produce the stretch corresponding to the read Mutation
            while (pos.r_pos != mut_second)
            {
                advance_pos(pos, forward, mut_second, false);
            }
            return true;
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

    std::tuple< shared_ptr< Read_Chunk >, shared_ptr< Contig_Entry > >
    Read_Chunk::make_chunk_from_cigar(const Cigar& cigar, Seq_Type* rf_ptr, const Seq_Type& qr)
    {
        assert(rf_ptr->size() == cigar.get_rf_len() or rf_ptr->size() >= cigar.get_rf_start() + cigar.get_rf_len());
        assert(qr.size() == cigar.get_qr_len() or qr.size() >= cigar.get_qr_start() + cigar.get_qr_len());

        // create objects with default constructor
        shared_ptr< Read_Chunk > chunk_sptr(new Read_Chunk());
        shared_ptr< Contig_Entry > ce_sptr(new Contig_Entry(rf_ptr, (rf_ptr->size() == cigar.get_rf_len()? cigar.get_rf_start() : 0)));

        // fix lengths and rc flags
        chunk_sptr->_c_start = cigar.get_rf_start();
        chunk_sptr->_c_len = cigar.get_rf_len();
        chunk_sptr->_r_start = cigar.get_qr_start();
        chunk_sptr->_r_len = cigar.get_qr_len();
        chunk_sptr->_rc = cigar.is_reversed();

        // fix cross-pointers
        chunk_sptr->set_ce_ptr(ce_sptr.get());
        ce_sptr->add_chunk(chunk_sptr.get());

        // construct mutations and store them in Contig_Entry object
        ce_sptr->mut_cont() = *(make_mutations_from_cigar(cigar, qr));

        // store pointers in the Read_Chunk object
        for (auto it = ce_sptr->get_mut_cont().begin(); it != ce_sptr->get_mut_cont().end(); ++it)
        {
            chunk_sptr->_mut_ptr_cont.push_back(&*it);
        }

        return std::make_tuple(chunk_sptr, ce_sptr);
    }

    tuple< shared_ptr< Read_Chunk >, shared_ptr< Contig_Entry > >
    Read_Chunk::make_chunk_from_cigar_and_chunks(const Cigar& cigar, const Read_Chunk& rc1, const Read_Chunk& rc2)
    {
        shared_ptr< Read_Chunk > chunk_sptr;
        shared_ptr< Contig_Entry > ce_sptr;
        std::tie(chunk_sptr, ce_sptr) = make_chunk_from_cigar(cigar, new Seq_Type(rc1.get_seq()), rc2.get_seq());
        // fix Read_Entry pointer
        chunk_sptr->_re_ptr = rc2._re_ptr;
        return std::make_tuple(chunk_sptr, ce_sptr);
    }

    std::tuple< Size_Type, Size_Type, size_t > Read_Chunk::get_cut_data(Size_Type brk, bool is_contig_brk) const
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

        return make_tuple(c_pos, r_pos, i);
    }

    std::tuple< const Mutation*, Size_Type, Size_Type > Read_Chunk::get_mutation_to_cut(
        Size_Type brk, bool is_contig_brk, const set< Size_Type >& brk_cont) const
    {
        Size_Type c_pos;
        Size_Type r_pos;
        size_t i;
        tie(c_pos, r_pos, i) = get_cut_data(brk, is_contig_brk);

        // if we reached brk, there is no mutation to be cut
        if ((is_contig_brk and c_pos == brk)
            or (not is_contig_brk and r_pos == brk))
            return make_tuple< const Mutation*, Size_Type, Size_Type >(NULL, 0, 0);

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
            return make_tuple(_mut_ptr_cont[i], brk - c_pos, (not _rc? sel_brk - r_pos : r_pos - sel_brk));
        else
            return make_tuple(_mut_ptr_cont[i], sel_brk - c_pos, (not _rc? brk - r_pos : r_pos - brk));
    }

    std::tuple< bool, shared_ptr< Read_Chunk > > Read_Chunk::apply_contig_split(
        Size_Type c_brk, const map< const Mutation*, const Mutation* >& mut_cptr_map, const Contig_Entry* ce_cptr)
    {
        // compute read chunk position corresponding to a contig cut at c_brk
        Read_Chunk_Pos pos;
        if (get_c_end() < c_brk)
        {
            pos = get_end_pos();
        }
        else // c_brk <= c_end
        {
            pos = get_start_pos();
            while (pos.c_pos < c_brk)
            {
                increment_pos(pos, c_brk, true);
            }
            assert(pos.c_pos == c_brk);
            // any mutations spanning the break must be cut prior to calling this method
            assert(pos.mut_offset == 0);
            // if there exists an insertion at c_brk that is not in the map, it will stay on the left side of the break
            if (pos.mut_idx < _mut_ptr_cont.size()
                and _mut_ptr_cont[pos.mut_idx]->get_start() == c_brk
                and _mut_ptr_cont[pos.mut_idx]->is_ins()
                and mut_cptr_map.count(_mut_ptr_cont[pos.mut_idx]) == 0)
            {
                increment_pos(pos);
                assert(pos.c_pos == c_brk);
            }
            assert(pos.mut_idx == _mut_ptr_cont.size() or mut_cptr_map.count(_mut_ptr_cont[pos.mut_idx]) > 0);
        }

        if (pos.r_pos == get_r_start() or pos.r_pos == get_r_end())
        {
            // the chunk stays in one piece, but ce_ptr & mutation pointers might change
            // also, drop any deletion mapped to an empty read chunk
            bool move_to_rhs = false;
            //assert(pos.mut_idx == 0 or pos.mut_idx == _mut_ptr_cont.size());
            if ((not _rc and pos.r_pos == get_r_start())
                or (_rc and pos.r_pos == get_r_end()))
            {
                assert((c_brk <= _c_start and pos.mut_idx == 0)
                       or (_c_start < c_brk and pos.mut_idx == 1 and _mut_ptr_cont[0]->is_del()));
                move_to_rhs = true;
                // fix contig coordinates
                _c_start = (c_brk < _c_start? _c_start - c_brk : 0);
                _c_len = (c_brk <= _c_start? _c_len : _c_len - (c_brk - _c_start));
                // fix mutation pointers
                vector< const Mutation* > tmp(_mut_ptr_cont.begin() + pos.mut_idx, _mut_ptr_cont.end()); // drop initial deletion, if any
                _mut_ptr_cont.clear();
                for (size_t j = 0; j < tmp.size(); ++j)
                {
                    assert(mut_cptr_map.count(tmp[j]) > 0);
                    _mut_ptr_cont.push_back(mut_cptr_map.find(tmp[j])->second);
                }
                // fix ce_ptr
                _ce_ptr = ce_cptr;
            }
            else
            {
                assert((not _rc and pos.r_pos == get_r_end())
                       or (_rc and pos.r_pos == get_r_start()));
                assert((get_c_end() <= c_brk and pos.mut_idx == _mut_ptr_cont.size())
                       or (c_brk < get_c_end() and pos.mut_idx == _mut_ptr_cont.size() - 1 and _mut_ptr_cont[pos.mut_idx]->is_del()));
                if (pos.mut_idx == _mut_ptr_cont.size() - 1)
                {
                    _mut_ptr_cont.resize(pos.mut_idx); // drop final deletion, if any
                }
            }
            return std::make_tuple(move_to_rhs, shared_ptr< Read_Chunk >(NULL));
        }
        else
        {
            // read chunk gets split in 2
            assert(get_r_start() < pos.r_pos and pos.r_pos < get_r_end());
            assert(get_c_start() <= pos.c_pos and pos.c_pos <= get_c_end());

            // create new read chunk for second part of the contig
            shared_ptr< Read_Chunk > rc_sptr(new Read_Chunk());
            rc_sptr->_rc = _rc;
            rc_sptr->_re_ptr = _re_ptr;

            rc_sptr->_ce_ptr = ce_cptr;
            rc_sptr->_c_start = 0;
            rc_sptr->_c_len = _c_start + _c_len - c_brk;
            _c_len -= rc_sptr->_c_len;

            if (not _rc)
            {
                rc_sptr->_r_start = pos.r_pos;
                rc_sptr->_r_len = _r_start + _r_len - pos.r_pos;
                _r_len -= rc_sptr->_r_len;
            }
            else
            {
                rc_sptr->_r_start = _r_start;
                rc_sptr->_r_len = pos.r_pos - _r_start;
                _r_start = pos.r_pos;
                _r_len -= rc_sptr->_r_len;
            }

            // transfer mutations starting at i (to values under mapping) to read chunk for contig rhs
            for (size_t j = pos.mut_idx; j < _mut_ptr_cont.size(); ++j)
            {
                assert(mut_cptr_map.count(_mut_ptr_cont[j]) > 0);
                rc_sptr->_mut_ptr_cont.push_back(mut_cptr_map.find(_mut_ptr_cont[j])->second);
            }

            // drop transfered mutations from this chunk
            _mut_ptr_cont.resize(pos.mut_idx);

            return std::make_tuple(false, rc_sptr);
        }
    }

    std::tuple< shared_ptr< Read_Chunk >, shared_ptr< Contig_Entry >, shared_ptr< Mutation_Trans_Cont > >
    Read_Chunk::reverse() const
    {
        shared_ptr< Read_Chunk > rc_sptr(new Read_Chunk());
        shared_ptr< Contig_Entry > ce_sptr(new Contig_Entry(new string(get_seq()), _r_start));
        shared_ptr< Mutation_Trans_Cont > mut_trans_cont_sptr(new Mutation_Trans_Cont());

        rc_sptr->_r_start = _c_start;
        rc_sptr->_r_len = _c_len;
        rc_sptr->_c_start = _r_start;
        rc_sptr->_c_len = _r_len;
        rc_sptr->_rc = _rc;

        // add cross-pointers
        rc_sptr->set_ce_ptr(ce_sptr.get());
        ce_sptr->add_chunk(rc_sptr.get());

        Read_Chunk_Pos pos = get_start_pos();
        while (pos != get_end_pos())
        {
            // ignore matched stretches
            while (get_match_len_from_pos(pos) > 0)
                increment_pos(pos);
            if (pos == get_end_pos())
                break;

            Read_Chunk_Pos pos_next = pos;
            increment_pos(pos_next);

            // create reversed Mutation
            Mutation rev_mut((not _rc? pos.r_pos : pos_next.r_pos),
                             (not _rc? pos_next.r_pos - pos.r_pos : pos.r_pos - pos_next.r_pos),
                             (not _rc? _ce_ptr->substr(pos.c_pos, pos_next.c_pos - pos.c_pos)
                              : reverseComplement(_ce_ptr->substr(pos.c_pos, pos_next.c_pos - pos.c_pos))));

            // save it in Mutation container
            Mutation_Cont::iterator rev_mut_it;
            bool success;
            tie(rev_mut_it, success) = ce_sptr->mut_cont().insert(rev_mut);
            assert(success);

            // save the pair in the Mutation translation container
            Mutation_Trans trans;
            trans.old_mut_cptr = _mut_ptr_cont[pos.mut_idx];
            trans.new_mut_cptr = &*rev_mut_it;
            Mutation_Trans_Cont::iterator it;
            tie(it, success) = mut_trans_cont_sptr->insert(trans);
            assert(success);
            pos = pos_next;
        }

        // save mutation pointers in the Read_Chunk object
        for (Mutation_Cont::iterator rev_mut_it = ce_sptr->get_mut_cont().begin(); rev_mut_it != ce_sptr->get_mut_cont().end(); ++rev_mut_it)
        {
            rc_sptr->_mut_ptr_cont.push_back(&*rev_mut_it);
        }

        return std::make_tuple(rc_sptr, ce_sptr, mut_trans_cont_sptr);
    }

    /*
    std::tuple< shared_ptr< Mutation_Trans_Cont >, shared_ptr< Mutation_Cont > >
    Read_Chunk::lift_read_mutations(const Mutation_Cont& r_mut_cont) const
    {
        shared_ptr< Mutation_Trans_Cont > mut_trans_cont_sptr(new Mutation_Trans_Cont());
        shared_ptr< Mutation_Cont > new_mut_cont_sptr(new Mutation_Cont());

        Read_Chunk_Pos pos = get_start_pos();
        Mutation_Cont::iterator r_mut_it = r_mut_cont.begin();
        Mutation_Cont::reverse_iterator r_mut_rit = r_mut_cont.rbegin();
        while (not _rc? r_mut_it != r_mut_cont.end() : r_mut_rit != r_mut_cont.rend())
        {
            const Mutation& r_mut = (not _rc? *r_mut_it : *r_mut_rit);

            // find stretch corresponding to read Mutation r_mut
            Read_Chunk_Pos pos_next;
            for (pos_next = pos; not advance_pos_til_mut(pos_next, r_mut); pos = pos_next);

            assert(_rc or (pos.r_pos == r_mut.get_start() and pos_next.r_pos == r_mut.get_end()));
            assert(not _rc or (pos_next.r_pos == r_mut.get_start() and pos.r_pos == r_mut.get_end()));

            Mutation m;
            if (r_mut.have_seq())
                m = Mutation(pos.c_pos, pos_next.c_pos, not _rc? r_mut.get_seq() : reverseComplement(r_mut.get_seq()));
            else
                m = Mutation(pos.c_pos, pos_next.c_pos, r_mut.get_seq_len());

            // insert translated Mutation into container
            Mutation_Cont::iterator new_mut_it;
            bool success;
            tie(new_mut_it, success) = new_mut_cont_sptr->insert(m);
            assert(success);

            // prepare translation entry
            Mutation_Trans trans;
            trans.old_mut_cptr = &r_mut;
            trans.new_mut_cptr = &*new_mut_it;

            // backtrack, and save the original contig mutation slices observed by this read mutation
            Read_Chunk_Pos pos_1 = pos;
            while (pos_1 != pos_next)
            {
                assert(pos_1.r_pos != pos_next.r_pos);
                Read_Chunk_Pos pos_2 = pos;
                increment_pos(pos_2, pos_next.r_pos, false);

                if (get_match_len_from_pos(pos_1) == 0)
                {
                    // a non-matched stretch follows
                    trans.new_mut_rev_list.push_back(make_tuple(
                        _mut_ptr_cont[pos_1.mut_idx], pos_1.mut_offset, (pos_2.mut_idx == pos_1.mut_idx? pos_2.mut_offset : 0)));
                }

                pos_1 = pos_2;
            }

            // insert new Mutation translation
            Mutation_Trans_Cont::iterator it;
            tie(it, success) = mut_trans_cont_sptr->insert(trans);
            assert(success);

            // we do not advance pos to pos_next, because read Mutations can be overlapping
            // next iteration starts again at pos, looking for the next read Mutation
            if (not _rc)
                ++r_mut_it;
            else
                ++r_mut_rit;
        }
        return make_tuple(mut_trans_cont_sptr, new_mut_cont_sptr);
    }
    */

    shared_ptr< vector< std::tuple< Mutation_CPtr, Size_Type, Size_Type, bool > > >
    Read_Chunk::get_mutations_under_mapping(const vector< Mutation_CPtr >& r_mut_cptr_cont,
                                            const Mutation_Trans_Cont& mut_map,
                                            shared_ptr< Mutation_Cont > new_mut_cont_sptr) const
    {
        (void)new_mut_cont_sptr;
        shared_ptr< vector< std::tuple< Mutation_CPtr, Size_Type, Size_Type, bool > > > res(
            new vector< std::tuple< Mutation_CPtr, Size_Type, Size_Type, bool > >());

        Read_Chunk_Pos pos = get_start_pos();
        size_t r_mut_cnt = 0;
        Mutation fake_mut((not _rc? get_r_end() : get_r_start()), (not _rc? get_r_end() : get_r_start()), 0);
        while (true)
        {
            size_t r_mut_idx = (not _rc? r_mut_cnt : r_mut_cptr_cont.size() - 1 - r_mut_cnt); // next new mutation to look for
            const Mutation& r_mut = (r_mut_cnt < r_mut_cptr_cont.size()? *r_mut_cptr_cont[r_mut_idx] : fake_mut);

            Read_Chunk_Pos pos_next;
            bool got_r_mut = advance_pos_til_mut(pos_next, r_mut);

            if (got_r_mut and r_mut_cnt == r_mut_cptr_cont.size())
                break;

            if (got_r_mut)
            {
                assert(mut_map.count(&r_mut) == 1);
                auto it = mut_map.find(&r_mut);
                assert(it->new_mut_cptr->get_start() == pos.c_pos and it->new_mut_cptr->get_end() == pos_next.c_pos);
                res->push_back(std::make_tuple(it->new_mut_cptr, 0, 0, true));
                ++r_mut_cnt;
            }
            else
            {
                res->push_back(std::make_tuple(
                    _mut_ptr_cont[pos.mut_idx], pos.mut_offset, (pos_next.mut_idx == pos.mut_idx? pos_next.mut_offset : 0), false));
            }
            pos = pos_next;
        }
        assert(pos == get_end_pos());

        return res;
    }

    shared_ptr< Read_Chunk > Read_Chunk::collapse_mapping(const Read_Chunk& rc2, Mutation_Cont& extra_mut_cont) const
    {
        assert(rc2.get_c_start() == get_r_start() and rc2.get_c_end() == get_r_end());
        shared_ptr< Read_Chunk > res(new Read_Chunk());

        res->_c_start = get_c_start();
        res->_c_len = get_c_len();
        res->_r_start = rc2.get_r_start();
        res->_r_len = rc2.get_r_len();
        res->_rc = (get_rc() != rc2.get_rc());
        res->_re_ptr = rc2._re_ptr;
        res->_ce_ptr = _ce_ptr;

        // traverse c1 left-to-right; traverse rc1 left-to-right if not _rc, r-to-l ow;
        Read_Chunk_Pos pos = get_start_pos();
        size_t r_mut_cnt = 0;
        Mutation fake_mut((not get_rc()? get_r_end() : get_r_start()), 0);
        Mutation m;
        vector< Mutation_CPtr > c_muts;
        vector< Mutation_CPtr > r_muts;
        while (true)
        {
            size_t r_mut_idx = (not get_rc()? r_mut_cnt : rc2.get_mut_ptr_cont().size() - 1 - r_mut_cnt); // next new mutation to look for
            const Mutation& r_mut = (r_mut_cnt < rc2.get_mut_ptr_cont().size()? *rc2.get_mut_ptr_cont()[r_mut_idx] : fake_mut);

            Read_Chunk_Pos pos_next = pos;
            bool got_r_mut = advance_pos_til_mut(pos_next, r_mut);

            if ((got_r_mut and r_mut_cnt == rc2.get_mut_ptr_cont().size())
                or (not got_r_mut and get_match_len_from_pos(pos) > 0))
            {
                // the stretch pos->pos_next is a match, or we are at the end;
                // consolidate any outstanding mutations
                m.simplify(get_ce_ptr()->substr(m.get_start(), m.get_len()));
                if (not m.is_empty())
                {
                    if (r_muts.size() == 0)
                    {
                        // no adjacent contig mutations
                        assert(c_muts.size() == 1);
                        res->mut_ptr_cont().push_back(c_muts[0]);
                    }
                    else
                    {
                        // new metamutation
                        // add it to extra Mutation container if it doesn't exist
                        Mutation_CPtr new_mut = add_mut_to_cont(extra_mut_cont, m);
                        res->_mut_ptr_cont.push_back(new_mut);
                        // NOTE: here we could save the metamutation composition
                    }
                }
                m = Mutation(pos_next.c_pos, 0);
                c_muts.clear();
                r_muts.clear();

                if (got_r_mut and r_mut_cnt == rc2.get_mut_ptr_cont().size())
                    break;
            }
            else
            {
                // the stretch pos->post_next is not a match and we are not at the end;
                // add mutation slice to m
                if (got_r_mut)
                {
                    // this was an entire read mutation
                    m.merge(Mutation(pos.c_pos, pos_next.c_pos - pos.c_pos,
                                     (not get_rc()? r_mut.get_seq() : reverseComplement(r_mut.get_seq()))));
                    r_muts.push_back(&r_mut);
                    ++r_mut_cnt;
                }
                else
                {
                    // this was a (possibly sliced) contig mutation
                    assert((pos_next.mut_idx == pos.mut_idx and pos_next.mut_offset > pos.mut_offset)
                           or (pos_next.mut_idx == pos.mut_idx + 1 and pos_next.mut_offset == 0));
                    const Mutation& c_mut = *_mut_ptr_cont[pos.mut_idx];
                    m.merge(Mutation(pos.c_pos, pos_next.c_pos - pos.c_pos,
                                     c_mut.get_seq().substr(
                                         min(pos.mut_offset, c_mut.get_seq_len()),
                                         min((pos_next.mut_offset == 0? string::npos : pos_next.mut_offset - pos.mut_offset), c_mut.get_seq_len()))));
                    c_muts.push_back(&c_mut);
                }
            }

            pos = pos_next;
        }
        assert(pos == get_end_pos());

        return res;
    }

    /*
    vector< std::tuple< bool, Read_Chunk_CPtr, Size_Type, Size_Type > > Read_Chunk::collapse_mutations(
        const Read_Chunk& rc1, const Mutation_Extra_Cont& rc1_me_cont, const Read_Chunk& rc2)
    {
        assert(rc1.get_r_start() == rc2.get_c_start() and rc1.get_r_len() == rc2.get_c_len());

        vector< std::tuple< bool, Read_Chunk_CPtr, Size_Type, Size_Type > > res;

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

    void Read_Chunk::reverse()
    {
        _rc = not _rc;
        _c_start = _ce_ptr->get_len() - (_c_start + _c_len);
        for (size_t i = 0; i < _mut_ptr_cont.size() / 2; ++i)
        {
            swap(_mut_ptr_cont[i], _mut_ptr_cont[_mut_ptr_cont.size() - 1 - i]);
        }
    }

    void Read_Chunk::merge_next(Read_Chunk_CPtr rc_next_cptr)
    {
        assert(rc_next_cptr != NULL);
        assert(rc_next_cptr->get_re_ptr() == get_re_ptr());
        assert(rc_next_cptr->get_ce_ptr() == get_ce_ptr());
        assert(rc_next_cptr->get_rc() == get_rc());
        assert(rc_next_cptr->get_r_start() == get_r_end());
        assert(get_rc() or rc_next_cptr->get_c_start() == get_c_end());
        assert(not get_rc() or rc_next_cptr->get_c_end() == get_c_start());

        // fix coordinates
        if (get_rc())
        {
            _c_start = rc_next_cptr->_c_start;
        }
        _c_len += rc_next_cptr->_c_len;
        _r_len += rc_next_cptr->_r_len;
    }

    void Read_Chunk::rebase(const Contig_Entry* ce_cptr, Size_Type prefix_len, const Mutation_Trans_Cont& mut_map)
    {
        _ce_ptr = ce_cptr;
        _c_start += prefix_len;
        vector< Mutation_CPtr > old_mut_cptr_cont = _mut_ptr_cont;
        _mut_ptr_cont.clear();
        for (auto old_mut_cptr_it = old_mut_cptr_cont.begin(); old_mut_cptr_it != old_mut_cptr_cont.end(); ++old_mut_cptr_it)
        {
            Mutation_CPtr old_mut_cptr = *old_mut_cptr_it;
            assert(mut_map.count(old_mut_cptr) == 1);
            Mutation_Trans_Cont::const_iterator it = mut_map.find(old_mut_cptr);
            _mut_ptr_cont.push_back(it->new_mut_cptr);
        }
    }

    ostream& operator << (ostream& os, const Read_Chunk& rhs)
    {
        os << indent::tab << "(Read_Chunk &=" << (void*)&rhs
           << indent::inc << indent::nl << "re_cptr=" << (void*)rhs._re_ptr
           << ",ce_cptr=" << (void*)rhs._ce_ptr
           << indent::nl << "r_start=" << rhs.get_r_start()
           << ",r_len=" << rhs.get_r_len()
           << ",c_start=" << rhs.get_c_start()
           << ",c_len=" << rhs.get_c_len()
           << ",rc=" << (int)rhs.get_rc()
           << indent::nl << "mut_cptr_cont:"
           << indent::inc << '\n';
        print_seq(os, rhs._mut_ptr_cont, indent::nl, indent::tab, '\n');
        os << indent::dec << indent::dec << indent::tab << ")\n";
        return os;
    }
}
