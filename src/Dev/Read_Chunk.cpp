//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Read_Chunk.hpp"

#include "Mutation.hpp"
#include "Read_Entry.hpp"
#include "Contig_Entry.hpp"
#include "print_seq.hpp"
#include "indent.hpp"
#include "../Util/Util.h"

using namespace std;


namespace MAC
{
    bool Read_Chunk::Pos::check() const
    {
        ASSERT(rc_cptr != NULL);
        ASSERT(rc_cptr->get_c_start() <= c_pos and c_pos <= rc_cptr->get_c_end());
        ASSERT(rc_cptr->get_r_start() <= r_pos and r_pos <= rc_cptr->get_r_end());
        ASSERT(mut_idx <= rc_cptr->_mut_ptr_cont.size());
        ASSERT(not mut_offset == 0
               or c_pos <= (mut_idx < rc_cptr->_mut_ptr_cont.size()?
                            rc_cptr->_mut_ptr_cont[mut_idx]->get_start()
                            : rc_cptr->get_c_end()));
        ASSERT(not mut_offset != 0
               or (mut_idx < rc_cptr->_mut_ptr_cont.size()
                   and c_pos == rc_cptr->_mut_ptr_cont[mut_idx]->get_start() + min(mut_offset, rc_cptr->_mut_ptr_cont[mut_idx]->get_len())));
        return true;
    }

    inline Size_Type Read_Chunk::Pos::get_match_len(bool forward) const
    {
        ASSERT(check());
        if (mut_offset != 0)
        {
            return 0;
        }
        if (forward)
        {
            return (mut_idx < rc_cptr->_mut_ptr_cont.size()? rc_cptr->_mut_ptr_cont[mut_idx]->get_start() : rc_cptr->get_c_end()) - c_pos;
        }
        else
        {
            return c_pos - (mut_idx > 0? rc_cptr->_mut_ptr_cont[mut_idx - 1]->get_end() : rc_cptr->get_c_start());
        }
    }

    void Read_Chunk::Pos::increment(Size_Type brk, bool on_contig)
    {
        ASSERT(check());
        ASSERT(*this != rc_cptr->get_end_pos());
        if (brk == 0)
        {
            // by default, use an unreachable breakpoint
            on_contig = true;
            brk = rc_cptr->get_c_end() + 1;
        }
        ASSERT(not on_contig or c_pos < brk);
        ASSERT(on_contig or (not rc_cptr->get_rc()? r_pos < brk : brk < r_pos));
        if (mut_idx < rc_cptr->_mut_ptr_cont.size()
            and c_pos == rc_cptr->_mut_ptr_cont[mut_idx]->get_start() + min(mut_offset, rc_cptr->_mut_ptr_cont[mut_idx]->get_len()))
        {
            // at the start or inside a mutation
            Size_Type c_leftover = (rc_cptr->_mut_ptr_cont[mut_idx]->get_len() > mut_offset?
                                    rc_cptr->_mut_ptr_cont[mut_idx]->get_len() - mut_offset
                                    : 0);
            Size_Type r_leftover = (rc_cptr->_mut_ptr_cont[mut_idx]->get_seq_len() > mut_offset?
                                    rc_cptr->_mut_ptr_cont[mut_idx]->get_seq_len() - mut_offset
                                    : 0);
            Size_Type skip_len = 0;
            if (on_contig and brk < c_pos + c_leftover)
            {
                skip_len = brk - c_pos;
            }
            else if (not on_contig and (not rc_cptr->get_rc()? brk < r_pos + r_leftover : r_pos - r_leftover < brk))
            {
                skip_len = (not rc_cptr->get_rc()? brk - r_pos : r_pos - brk);
            }
            else
            {
                // breakpoint does not occur later in this mutation: we advance past it entirely
                skip_len = max(c_leftover, r_leftover);
            }
            c_pos += min(c_leftover, skip_len);
            r_pos = (not rc_cptr->get_rc()? r_pos + min(r_leftover, skip_len) : r_pos - min(r_leftover, skip_len));
            mut_offset += skip_len;
            if (mut_offset >= rc_cptr->_mut_ptr_cont[mut_idx]->get_len() and mut_offset >= rc_cptr->_mut_ptr_cont[mut_idx]->get_seq_len())
            {
                mut_offset = 0;
                ++mut_idx;
            }
        }
        else
        {
            // a match stretch follows
            Size_Type match_len = get_match_len(true);
            ASSERT(match_len > 0);
            Size_Type skip_len = 0;
            if (on_contig and brk < c_pos + match_len)
            {
                skip_len = brk - c_pos;
            }
            else if (not on_contig and (not rc_cptr->get_rc()? brk < r_pos + match_len : r_pos - match_len < brk))
            {
                skip_len = (not rc_cptr->get_rc()? brk - r_pos : r_pos - brk);
            }
            else
            {
                skip_len = match_len;
            }
            c_pos += skip_len;
            r_pos = (not rc_cptr->get_rc()? r_pos + skip_len : r_pos - skip_len);
        }
    }

    void Read_Chunk::Pos::decrement(Size_Type brk, bool on_contig)
    {
        ASSERT(check());
        ASSERT(*this != rc_cptr->get_start_pos());
        if (brk == 0)
        {
            on_contig = true;
        }
        ASSERT(on_contig? brk <= c_pos : (not rc_cptr->get_rc()? brk < r_pos : r_pos < brk));
        if (mut_offset != 0
            or (mut_offset == 0 and mut_idx > 0
                and c_pos == rc_cptr->_mut_ptr_cont[mut_idx - 1]->get_end()))
        {
            // at the end or inside a mutation
            if (mut_offset == 0)
            {
                --mut_idx;
                mut_offset = max(rc_cptr->_mut_ptr_cont[mut_idx]->get_len(), rc_cptr->_mut_ptr_cont[mut_idx]->get_seq_len());
            }
            ASSERT(mut_offset > 0);
            Size_Type c_leftover = min(rc_cptr->_mut_ptr_cont[mut_idx]->get_len(), mut_offset);
            Size_Type r_leftover = min(rc_cptr->_mut_ptr_cont[mut_idx]->get_seq_len(), mut_offset);
            Size_Type c_skip_len = 0;
            Size_Type r_skip_len = 0;
            if (on_contig and c_pos - c_leftover < brk)
            {
                c_skip_len = c_pos - brk;
                r_skip_len = (r_leftover > c_leftover - c_skip_len? r_leftover - (c_leftover - c_skip_len) : 0);
            }
            else if (not on_contig and (not rc_cptr->get_rc()? r_pos - r_leftover < brk : brk < r_pos + r_leftover))
            {
                r_skip_len = (not rc_cptr->get_rc()? r_pos - brk : brk - r_pos);
                c_skip_len = (c_leftover > r_leftover - r_skip_len? c_leftover - (r_leftover - r_skip_len) : 0);
            }
            else
            {
                // breakpoint does not occur earlier in this mutation: we advance past it entirely
                c_skip_len = c_leftover;
                r_skip_len = r_leftover;
            }
            c_pos -= c_skip_len;
            r_pos = (not rc_cptr->get_rc()? r_pos - r_skip_len : r_pos + r_skip_len);
            mut_offset = max(c_leftover - c_skip_len, r_leftover - r_skip_len);
        }
        else
        {
            // a match stretch follows
            Size_Type match_len = get_match_len(false);
            ASSERT(match_len > 0);
            Size_Type skip_len = 0;
            if (on_contig and c_pos - match_len < brk)
            {
                skip_len = c_pos - brk;
            }
            else if (not on_contig and (not rc_cptr->get_rc()? r_pos - match_len < brk : brk < r_pos + match_len))
            {
                skip_len = (not rc_cptr->get_rc()? r_pos - brk : brk - r_pos);
            }
            else
            {
                skip_len = match_len;
            }
            c_pos -= skip_len;
            r_pos = (not rc_cptr->get_rc()? r_pos - skip_len : r_pos + skip_len);
        }
    }

    void Read_Chunk::Pos::jump_to_brk(Size_Type brk, bool on_contig)
    {
        ASSERT(check());
        ASSERT(not on_contig or (rc_cptr->get_c_start() <= brk and brk <= rc_cptr->get_c_end()));
        ASSERT(on_contig or (rc_cptr->get_r_start() <= brk and brk <= rc_cptr->get_r_end()));
        bool forward = (on_contig?
                        c_pos <= brk
                        : (not rc_cptr->get_rc()) == (r_pos <= brk));
        while ((on_contig and c_pos != brk) or (not on_contig and r_pos != brk))
        {
            advance(forward, brk, on_contig);
        }
    }

    bool Read_Chunk::Pos::advance_til_mut(const Mutation& mut, bool forward)
    {
        Size_Type mut_first;
        Size_Type mut_second;
        if (forward == not rc_cptr->get_rc())
        {
            mut_first = mut.get_start();
            mut_second = mut.get_end();
            ASSERT(r_pos <= mut_first);
        }
        else
        {
            mut_first = mut.get_end();
            mut_second = mut.get_start();
            ASSERT(mut_first <= r_pos);
        }
        if (r_pos != mut_first)
        {
            // there are still mapping stretches before the read Mutation
            // advance using mut_first as read breakpoint
            advance(forward, mut_first, false);
            return false;
        }
        else
        {
            // r_pos == mut_first
            // at this point, we produce the stretch corresponding to the read Mutation
            jump_to_brk(mut_second, false);
            return true;
        }
    }

    bool Read_Chunk::have_mutation(const Mutation* mut_cptr) const
    {
        for (size_t i = 0; i < _mut_ptr_cont.size(); ++i)
        {
            if (_mut_ptr_cont[i] == mut_cptr)
            {
                return true;
            }
        }
        return false;
    }

    void Read_Chunk::cond_add_mutation(const Mutation* old_mut_cptr, const Mutation* new_mut_cptr)
    {
        size_t i;
        for (i = 0; i < _mut_ptr_cont.size(); ++i)
        {
            if (_mut_ptr_cont[i] == old_mut_cptr)
            {
                _mut_ptr_cont.insert(_mut_ptr_cont.begin() + i + 1, new_mut_cptr);
                break;
            }
        }
    }

    std::tuple< shared_ptr< Read_Chunk >, shared_ptr< Contig_Entry > >
    Read_Chunk::make_chunk_from_cigar(const Cigar& cigar, Seq_Type* rf_ptr, const Seq_Type& qr)
    {
        ASSERT(rf_ptr->size() == cigar.get_rf_len() or rf_ptr->size() >= cigar.get_rf_start() + cigar.get_rf_len());
        ASSERT(qr.size() == cigar.get_qr_len() or qr.size() >= cigar.get_qr_start() + cigar.get_qr_len());

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

    /*
    std::tuple< Size_Type, Size_Type, size_t > Read_Chunk::get_cut_data(Size_Type brk, bool is_contig_brk) const
    {
        ASSERT(not is_contig_brk or (get_c_start() <= brk and brk <= get_c_end()));
        ASSERT(is_contig_brk or (get_r_start() <= brk and brk <= get_r_end()));

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
            ASSERT(c_pos <= mut_cptr->get_start());

            // progress
            ASSERT(not is_contig_brk or c_pos <= brk);
            ASSERT(is_contig_brk or _rc or r_pos <= brk);
            ASSERT(is_contig_brk or not _rc or brk <= r_pos);

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
            ASSERT(not is_contig_brk or c_pos < brk);
            ASSERT(is_contig_brk or _rc or r_pos < brk);
            ASSERT(is_contig_brk or not _rc or brk < r_pos);

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

        ASSERT(get_c_start() <= c_pos and c_pos <= get_c_end());
        ASSERT(get_r_start() <= r_pos and r_pos <= get_r_end());
        ASSERT(not is_contig_brk or c_pos <= brk);
        ASSERT(is_contig_brk or _rc or r_pos <= brk);
        ASSERT(is_contig_brk or not _rc or brk <= r_pos);

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
        ASSERT(i < _mut_ptr_cont.size());
        ASSERT(c_pos == _mut_ptr_cont[i]->get_start());

        ASSERT(not is_contig_brk or (c_pos < brk and brk < _mut_ptr_cont[i]->get_end()));
        ASSERT(is_contig_brk or _rc or (r_pos < brk and brk < r_pos + _mut_ptr_cont[i]->get_seq_len()));
        ASSERT(is_contig_brk or not _rc or (r_pos > brk and brk > r_pos - _mut_ptr_cont[i]->get_seq_len()));

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
    */

    std::tuple< bool, shared_ptr< Read_Chunk > > Read_Chunk::split(
        Size_Type c_brk, const map< const Mutation*, const Mutation* >& mut_cptr_map, const Contig_Entry* ce_cptr)
    {
        // compute read chunk position corresponding to a contig cut at c_brk
        Pos pos;
        if (get_c_end() < c_brk)
        {
            pos = get_end_pos();
        }
        else if (c_brk < get_c_start())
        {
            pos = get_start_pos();
        }
        else // c_start <= c_brk <= c_end
        {
            pos = get_start_pos();
            pos.jump_to_brk(c_brk, true);
            ASSERT(pos.c_pos == c_brk);
            // any mutations spanning the break must be cut prior to calling this method
            ASSERT(pos.mut_offset == 0);
            // if there exists an insertion at c_brk that is not in the map, it will stay on the left side of the break
            while (pos.mut_idx < _mut_ptr_cont.size()
                and _mut_ptr_cont[pos.mut_idx]->get_start() == c_brk
                and _mut_ptr_cont[pos.mut_idx]->is_ins()
                and mut_cptr_map.count(_mut_ptr_cont[pos.mut_idx]) == 0)
            {
                pos.increment();
                ASSERT(pos.c_pos == c_brk);
            }
            ASSERT(pos.mut_idx == _mut_ptr_cont.size() or mut_cptr_map.count(_mut_ptr_cont[pos.mut_idx]) > 0);
        }

        if (pos.r_pos == get_r_start() or pos.r_pos == get_r_end())
        {
            // the chunk stays in one piece, but ce_ptr & mutation pointers might change
            // also, drop any deletion mapped to an empty read chunk
            bool move_to_rhs = false;
            //ASSERT(pos.mut_idx == 0 or pos.mut_idx == _mut_ptr_cont.size());
            if ((not _rc and pos.r_pos == get_r_start())
                or (_rc and pos.r_pos == get_r_end()))
            {
                ASSERT((c_brk <= _c_start and pos.mut_idx == 0)
                       or (_c_start < c_brk and pos.mut_idx == 1 and _mut_ptr_cont[0]->is_del()));
                move_to_rhs = true;
                // fix contig coordinates
                _c_len = (c_brk <= _c_start? _c_len : _c_len - (c_brk - _c_start));
                _c_start = (c_brk <= _c_start? _c_start - c_brk : 0);
                /*
                if (pos.mut_idx == 1)
                {
                    // skip initial deletion
                    _c_start += _mut_ptr_cont[0]->get_len();
                    _c_len -= _mut_ptr_cont[0]->get_len();
                }
                */
                // fix mutation pointers
                vector< const Mutation* > tmp(_mut_ptr_cont.begin() + pos.mut_idx, _mut_ptr_cont.end()); // drop initial deletion, if any
                _mut_ptr_cont.clear();
                for (size_t j = 0; j < tmp.size(); ++j)
                {
                    ASSERT(mut_cptr_map.count(tmp[j]) == 1);
                    _mut_ptr_cont.push_back(mut_cptr_map.find(tmp[j])->second);
                }
                // fix ce_ptr
                _ce_ptr = ce_cptr;
            }
            else
            {
                ASSERT((not _rc and pos.r_pos == get_r_end())
                       or (_rc and pos.r_pos == get_r_start()));
                ASSERT((get_c_end() <= c_brk and pos.mut_idx == _mut_ptr_cont.size())
                       or (c_brk < get_c_end() and pos.mut_idx == _mut_ptr_cont.size() - 1 and _mut_ptr_cont[pos.mut_idx]->is_del()));
                if (pos.mut_idx == _mut_ptr_cont.size() - 1)
                {
                    // drop final deletion
                    _c_len -= _mut_ptr_cont[pos.mut_idx]->get_len();
                    _mut_ptr_cont.resize(pos.mut_idx);
                }
            }
            return std::make_tuple(move_to_rhs, shared_ptr< Read_Chunk >(NULL));
        }
        else
        {
            // read chunk gets split in 2
            ASSERT(get_r_start() < pos.r_pos and pos.r_pos < get_r_end());
            ASSERT(get_c_start() <= pos.c_pos and pos.c_pos <= get_c_end());

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
                ASSERT(mut_cptr_map.count(_mut_ptr_cont[j]) > 0);
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

        Pos pos = get_start_pos();
        while (pos != get_end_pos())
        {
            // ignore matched stretches
            while (pos.get_match_len() > 0)
            {
                pos.increment();
            }
            if (pos == get_end_pos())
            {
                break;
            }

            Pos pos_next = pos;
            pos_next.increment();

            // create reversed Mutation
            Mutation rev_mut((not _rc? pos.r_pos : pos_next.r_pos),
                             (not _rc? pos_next.r_pos - pos.r_pos : pos.r_pos - pos_next.r_pos),
                             (not _rc? _ce_ptr->substr(pos.c_pos, pos_next.c_pos - pos.c_pos)
                              : reverseComplement(_ce_ptr->substr(pos.c_pos, pos_next.c_pos - pos.c_pos))));

            // save it in Mutation container
            Mutation_Cont::iterator rev_mut_it;
            bool success;
            tie(rev_mut_it, success) = ce_sptr->mut_cont().insert(rev_mut);
            ASSERT(success);

            // save the pair in the Mutation translation container
            Mutation_Trans trans;
            trans.old_mut_cptr = _mut_ptr_cont[pos.mut_idx];
            trans.new_mut_cptr = &*rev_mut_it;
            Mutation_Trans_Cont::iterator it;
            tie(it, success) = mut_trans_cont_sptr->insert(trans);
            ASSERT(success);
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

            ASSERT(_rc or (pos.r_pos == r_mut.get_start() and pos_next.r_pos == r_mut.get_end()));
            ASSERT(not _rc or (pos_next.r_pos == r_mut.get_start() and pos.r_pos == r_mut.get_end()));

            Mutation m;
            if (r_mut.have_seq())
                m = Mutation(pos.c_pos, pos_next.c_pos, not _rc? r_mut.get_seq() : reverseComplement(r_mut.get_seq()));
            else
                m = Mutation(pos.c_pos, pos_next.c_pos, r_mut.get_seq_len());

            // insert translated Mutation into container
            Mutation_Cont::iterator new_mut_it;
            bool success;
            tie(new_mut_it, success) = new_mut_cont_sptr->insert(m);
            ASSERT(success);

            // prepare translation entry
            Mutation_Trans trans;
            trans.old_mut_cptr = &r_mut;
            trans.new_mut_cptr = &*new_mut_it;

            // backtrack, and save the original contig mutation slices observed by this read mutation
            Read_Chunk_Pos pos_1 = pos;
            while (pos_1 != pos_next)
            {
                ASSERT(pos_1.r_pos != pos_next.r_pos);
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
            ASSERT(success);

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

    /*
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
                ASSERT(mut_map.count(&r_mut) == 1);
                auto it = mut_map.find(&r_mut);
                ASSERT(it->new_mut_cptr->get_start() == pos.c_pos and it->new_mut_cptr->get_end() == pos_next.c_pos);
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
        ASSERT(pos == get_end_pos());

        return res;
    }
    */

    shared_ptr< Read_Chunk > Read_Chunk::collapse_mapping(const Read_Chunk& rc2, Mutation_Cont& extra_mut_cont) const
    {
        ASSERT(get_r_start() <= rc2.get_c_start() and rc2.get_c_end() <= get_r_end());
        shared_ptr< Read_Chunk > res(new Read_Chunk());

        res->_r_start = rc2.get_r_start();
        res->_r_len = rc2.get_r_len();
        res->_rc = (get_rc() != rc2.get_rc());
        res->_re_ptr = rc2._re_ptr;
        res->_ce_ptr = _ce_ptr;

        // traverse c1 left-to-right; traverse rc1 left-to-right if not _rc, r-to-l ow;
        Pos pos = get_start_pos();
        size_t r_mut_cnt = 0;
        Mutation fake_mut_start((not get_rc()? rc2.get_c_start() : rc2.get_c_end()), 0);
        Mutation fake_mut_end((not get_rc()? rc2.get_c_end() : rc2.get_c_start()), 0);
        bool past_start = false;
        Mutation m;
        vector< Mutation_CPtr > c_muts;
        vector< Mutation_CPtr > r_muts;
        while (true)
        {
            size_t r_mut_idx = (not get_rc()? r_mut_cnt : rc2.get_mut_ptr_cont().size() - 1 - r_mut_cnt); // next new mutation to look for
            const Mutation& r_mut = (r_mut_cnt == 0 and not past_start?
                                     fake_mut_start :
                                     r_mut_cnt < rc2.get_mut_ptr_cont().size()? *rc2.get_mut_ptr_cont()[r_mut_idx] : fake_mut_end);

            Pos pos_next = pos;
            bool got_r_mut = pos_next.advance_til_mut(r_mut);

            if (not past_start)
            {
                if (got_r_mut)
                {
                    ASSERT(pos_next == pos);
                    res->_c_start = pos.c_pos;
                    past_start = true;
                    if (pos == get_start_pos())
                    {
                        // if rc2 starts on contig break, incorporate initial deletions
                        Read_Chunk::Pos tmp_pos = pos_next;
                        while (tmp_pos != get_end_pos())
                        {
                            tmp_pos.increment();
                            if (tmp_pos.r_pos != pos_next.r_pos)
                            {
                                break;
                            }
                            ASSERT(tmp_pos.r_pos == pos_next.r_pos);
                            ASSERT(tmp_pos.mut_offset == 0);
                            ASSERT(tmp_pos.mut_idx == pos_next.mut_idx + 1);
                            // initial deletions are always whole (offset==0 in pos_next and tmp_pos)
                            Mutation_CPtr c_mut_cptr = _mut_ptr_cont[pos_next.mut_idx];
                            m.merge(*c_mut_cptr);
                            c_muts.push_back(c_mut_cptr);
                            pos_next = tmp_pos;
                        }
                    }
                    else
                    {
                        // if rc2 doesn't start on contig break, skip initial deletions
                        Read_Chunk::Pos tmp_pos = pos_next;
                        tmp_pos.advance_past_del();
                        // FIX: notorious situation:
                        // do not skip deletions if we can reach the contig end on deletins only,
                        // and read is aligned on contig end
                        if (tmp_pos != get_end_pos())
                        {
                            pos_next = tmp_pos;
                            res->_c_start = pos_next.c_pos;
                        }
                        else
                        {
                            // rc2 contains insertions only, and is mapped to the contig end
                            // in this case, we leave the deletion, as it is needed
                            // to make the result aligned on the contig end
                            ASSERT(rc2.get_c_start() == rc2.get_c_end());
                        }
                    }
                }
            }
            else if ((not got_r_mut and pos.get_match_len() > 0)
                or (got_r_mut and r_mut_cnt == rc2.get_mut_ptr_cont().size()))
            {
                ASSERT(not got_r_mut or pos_next == pos);
                if (got_r_mut and pos_next.r_pos == get_end_pos().r_pos)
                {
                    // if rc2 ends on contig end, incorporate remaining deletions
                    Read_Chunk::Pos tmp_pos = pos_next;
                    while (tmp_pos != get_end_pos())
                    {
                        tmp_pos.increment();
                        ASSERT(tmp_pos.r_pos == pos_next.r_pos);
                        ASSERT(tmp_pos.mut_offset == 0);
                        ASSERT(tmp_pos.mut_idx == pos_next.mut_idx + 1);
                        // final deletions are whole, or we have used a read mutation
                        ASSERT(pos_next.mut_offset == 0 or r_muts.size() > 0);
                        Mutation_CPtr c_mut_cptr = _mut_ptr_cont[pos_next.mut_idx];
                        m.merge(Mutation(pos_next.c_pos, tmp_pos.c_pos - pos_next.c_pos, 0));
                        c_muts.push_back(c_mut_cptr);
                        pos_next = tmp_pos;
                    }

                }
                // the stretch pos->pos_next is a match, or we are at the end;
                // consolidate any outstanding mutations
                m.simplify(get_ce_ptr()->substr(m.get_start(), m.get_len()));
                if (not m.is_empty())
                {
                    /*
                    if (r_muts.size() == 0)
                    {
                        // no adjacent contig mutations
                        ASSERT(c_muts.size() == 1);
                        res->mut_ptr_cont().push_back(c_muts[0]);
                    }
                    else
                    */
                    {
                        // new metamutation
                        // add it to extra Mutation container if it doesn't exist
                        Mutation_CPtr new_mut = add_mut_to_cont(extra_mut_cont, m);
                        res->_mut_ptr_cont.push_back(new_mut);
                        // NOTE: here we could save the metamutation composition
                    }
                }
                m = Mutation();
                c_muts.clear();
                r_muts.clear();

                if (got_r_mut)
                {
                    res->_c_len = pos_next.c_pos - res->_c_start;
                    break;
                }
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
                    ASSERT((pos_next.mut_idx == pos.mut_idx and pos_next.mut_offset > pos.mut_offset)
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
        ASSERT(not (rc2.get_c_start() == get_r_start())
               or (not get_rc()? res->get_c_start() == get_c_start() : res->get_c_end() == get_c_end()));
        ASSERT(not (rc2.get_c_end() == get_r_end())
               or (not get_rc()? res->get_c_end() == get_c_end() : res->get_c_start() == get_c_start()));
        return res;
    }

    /*
    vector< std::tuple< bool, Read_Chunk_CPtr, Size_Type, Size_Type > > Read_Chunk::collapse_mutations(
        const Read_Chunk& rc1, const Mutation_Extra_Cont& rc1_me_cont, const Read_Chunk& rc2)
    {
        ASSERT(rc1.get_r_start() == rc2.get_c_start() and rc1.get_r_len() == rc2.get_c_len());

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
            ASSERT(not rc1.get_rc()? r1_pos <= rc1.get_r_end() : rc1.get_r_start() <= r1_pos);
            // i2 < n2 and r1 not reversed => before i2 mutation
            ASSERT(not (i2 < rc2._mut_ptr_cont.size() and not rc1.get_rc()) or r1_pos <= rc2._mut_ptr_cont[i2]->get_start());
            // 0 <= i2 - 1 < n2 and r1 not reversed => after i2 - 1 mutation
            ASSERT(not (0 < i2 and i2 < rc2._mut_ptr_cont.size() + 1 and not rc1.get_rc()) or rc2._mut_ptr_cont[i2 - 1]->get_end());
            // i2 < n2 and r1 reversed => after i2 mutation
            ASSERT(not (i2 < rc2._mut_ptr_cont.size() and rc1.get_rc()) or rc2._mut_ptr_cont[i2]->get_end() <= r1_pos);
            // 0 <= i2 - 1 < n2 and r1 reversed => before i2 - 1 mutation
            ASSERT(not (0 < i2 and i2 < rc2._mut_ptr_cont.size() + 1 and rc1.get_rc()) or r1_pos <= rc2._mut_ptr_cont[i2 - 1]->get_start());

            ASSERT(c_pos <= rc1.get_c_end());
            // if k1 != 0, we must be in the middle of a contig mutation
            ASSERT(not k1 != 0 or (i1 < rc1._mut_ptr_cont.size() and c_pos == rc1._mut_ptr_cont[i1]->get_start() + k1));
            // if k1 == 0 and i1 not at the end, the next contig mutation is at the right
            ASSERT(not (k1 == 0 and i1 < rc1._mut_ptr_cont.size()) or c_pos <= rc1._mut_ptr_cont[i1]->get_start());

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
                        ASSERT(k1 == 0);
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
                        ASSERT(i1 < rc1._mut_ptr_cont.size());
                        ASSERT(c_pos + c_span == rc1._mut_ptr_cont[i1]->get_start() + k1);

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
                            ASSERT(r1_span == rc2._mut_ptr_cont[i2]->get_len());
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
                ASSERT(rc1._mut_ptr_cont[i1]->get_len() > k1);
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

    /*
    vector< Size_Type > Read_Chunk::get_mut_pos() const
    {
        vector< Size_Type > res;
        for (Pos pos = get_start_pos(); not (pos == get_end_pos()); increment_pos(pos))
        {
            // if no breakpoints are used, we should never stop in the middle of a mutation
            ASSERT(pos.mut_offset == 0);
            if (get_match_len_from_pos(pos) == 0)
            {
                // at the start of a mutation
                res.push_back(pos.r_pos);
            }
        }
        ASSERT(res.size() == _mut_ptr_cont.size());
        return res;
    }
    */

    void Read_Chunk::reverse()
    {
        _rc = not _rc;
        _c_start = _ce_ptr->get_len() - (_c_start + _c_len);
        for (size_t i = 0; i < _mut_ptr_cont.size() / 2; ++i)
        {
            swap(_mut_ptr_cont[i], _mut_ptr_cont[_mut_ptr_cont.size() - 1 - i]);
        }
    }

    void Read_Chunk::merge_next(Read_Chunk_CPtr rc_next_cptr, Mutation::add_mut_mod_type add_mut_mod)
    {
        ASSERT(rc_next_cptr != NULL);
        ASSERT(rc_next_cptr->get_re_ptr() == get_re_ptr());
        ASSERT(rc_next_cptr->get_ce_ptr() == get_ce_ptr());
        ASSERT(rc_next_cptr->get_rc() == get_rc());
        ASSERT(rc_next_cptr->get_r_start() == get_r_end());
        ASSERT(get_rc() or rc_next_cptr->get_c_start() == get_c_end());
        ASSERT(not get_rc() or rc_next_cptr->get_c_end() == get_c_start());

        // fix coordinates
        if (get_rc())
        {
            _c_start = rc_next_cptr->_c_start;
        }
        _c_len += rc_next_cptr->_c_len;
        _r_len += rc_next_cptr->_r_len;
        // aquire mutations
        // if there are touching mutations across the break, merge them
        if (not get_rc())
        {
            size_t lhs_start = 0;
            if (_mut_ptr_cont.size() > 0 and rc_next_cptr->_mut_ptr_cont.size() > 0 and
                _mut_ptr_cont.back()->get_end() == rc_next_cptr->_mut_ptr_cont.front()->get_start())
            {
                Mutation m(*_mut_ptr_cont.back());
                m.merge(*rc_next_cptr->_mut_ptr_cont.front());
                Mutation_CPtr m_cptr = add_mut_mod(m);
                _mut_ptr_cont.back() = m_cptr;
                ++lhs_start;
            }
            _mut_ptr_cont.insert(_mut_ptr_cont.end(), rc_next_cptr->_mut_ptr_cont.begin() + lhs_start, rc_next_cptr->_mut_ptr_cont.end());
        }
        else // _rc
        {
            size_t rhs_size = _mut_ptr_cont.size();
            size_t lhs_size = rc_next_cptr->_mut_ptr_cont.size();
            if (lhs_size > 0)
            {
                // mutations from next chunk must be inserted before the ones in here
                if (rhs_size > 0
                    and rc_next_cptr->_mut_ptr_cont.back()->get_end() == _mut_ptr_cont.front()->get_start())
                {
                    Mutation m(*rc_next_cptr->_mut_ptr_cont.back());
                    m.merge(*_mut_ptr_cont.front());
                    Mutation_CPtr m_cptr = add_mut_mod(m);
                    _mut_ptr_cont.front() = m_cptr;
                    --lhs_size;
                }
                _mut_ptr_cont.resize(lhs_size + rhs_size, NULL);
                for (size_t i = 0; i < rhs_size; ++i)
                    _mut_ptr_cont[lhs_size + rhs_size - 1 - i] = _mut_ptr_cont[rhs_size - 1 - i];
                for (size_t i = 0; i <  lhs_size; ++i)
                    _mut_ptr_cont[i] = rc_next_cptr->_mut_ptr_cont[i];
            }
        }
    }

    void Read_Chunk::rebase(const Contig_Entry* ce_cptr, const Mutation_Trans_Cont& mut_map, Size_Type prefix_len)
    {
        _ce_ptr = ce_cptr;
        _c_start += prefix_len;
        vector< Mutation_CPtr > old_mut_cptr_cont = _mut_ptr_cont;
        _mut_ptr_cont.clear();
        for (auto old_mut_cptr_it = old_mut_cptr_cont.begin(); old_mut_cptr_it != old_mut_cptr_cont.end(); ++old_mut_cptr_it)
        {
            Mutation_CPtr old_mut_cptr = *old_mut_cptr_it;
            ASSERT(mut_map.count(old_mut_cptr) == 1);
            Mutation_Trans_Cont::const_iterator it = mut_map.find(old_mut_cptr);
            _mut_ptr_cont.push_back(it->new_mut_cptr);
        }
    }

    Seq_Type Read_Chunk::get_seq() const
    {
        Seq_Type res;
        Pos pos = (not _rc? get_start_pos() : get_end_pos());
        while (pos != (not _rc? get_end_pos() : get_start_pos()))
        {
            Pos next_pos = pos;
            next_pos.advance(not _rc);
            Size_Type match_len = pos.get_match_len(not _rc);
            if (match_len > 0)
            {
                string tmp = _ce_ptr->substr((not _rc? pos.c_pos : next_pos.c_pos) - _ce_ptr->get_seq_offset(), match_len);
                res += (not _rc? tmp : reverseComplement(tmp));
            }
            else if (pos.r_pos != next_pos.r_pos)
            {
                if (not _rc)
                {
                    res += _mut_ptr_cont[pos.mut_idx]->get_seq().substr(pos.mut_offset, next_pos.r_pos - pos.r_pos);
                }
                else
                {
                    res += reverseComplement(_mut_ptr_cont[next_pos.mut_idx]->get_seq().substr(next_pos.mut_offset, pos.r_pos - next_pos.r_pos));
                }
            }
            pos = next_pos;
        }
        return res;
    }

    Seq_Type Read_Chunk::substr(Size_Type start, Size_Type len) const
    {
        ASSERT(start >= _r_start and start + len <= _r_start + _r_len);
        return get_seq().substr(start - _r_start, len);
    }

    bool Read_Chunk::check() const
    {
        // no empty chunks
        ASSERT(get_r_len() > 0);
        // contigs coordinates
        ASSERT(get_c_start() <= get_c_end());
        ASSERT(get_c_start() <= get_ce_ptr()->get_seq_offset() + get_ce_ptr()->get_len());
        ASSERT(get_c_end() <= get_ce_ptr()->get_seq_offset() + get_ce_ptr()->get_len());
        // mapped length
        Size_Type c_len = get_c_end() - get_c_start();
        Size_Type r_len = get_r_end() - get_r_start();
        long long delta = 0;
        for (size_t i = 0; i < get_mut_ptr_cont().size(); ++i)
        {
            // no empty mutations
            ASSERT(not get_mut_ptr_cont()[i]->is_empty());
            // mutations must be in contig order
            ASSERT(i == 0 or get_mut_ptr_cont()[i - 1]->get_end() <= get_mut_ptr_cont()[i]->get_start());
#ifndef ALLOW_CONSECUTIVE_MUTATIONS
            ASSERT(i == 0 or get_mut_ptr_cont()[i - 1]->get_end() < get_mut_ptr_cont()[i]->get_start());
#endif
            delta += (long long)get_mut_ptr_cont()[i]->get_seq_len() - (long long)get_mut_ptr_cont()[i]->get_len();
        }
        ASSERT((long long)c_len + delta == (long long)r_len);
#ifndef ALLOW_PART_MAPPED_INNER_CHUNKS
        // chunks must end on contig breaks except for first and last
        ASSERT(get_r_start() == 0
               or (not get_rc()?
                   get_c_start() == 0
                   : get_c_end() == get_ce_ptr()->get_seq_offset() + get_ce_ptr()->get_len()));
        ASSERT(get_r_end() == get_re_ptr()->get_len()
               or (not get_rc()?
                   get_c_end() == get_ce_ptr()->get_seq_offset() + get_ce_ptr()->get_len()
                   : get_c_start() == 0));
#endif
        return true;
    }

    ostream& operator << (ostream& os, const Read_Chunk::Pos& pos)
    {
        os << "(c_pos=" << (size_t)pos.c_pos << ",r_pos=" << (size_t)pos.r_pos
           << ",mut_idx=" << pos.mut_idx << ",mut_offset=" << (size_t)pos.mut_offset << ")";
        return os;
    }

    ostream& operator << (ostream& os, const Read_Chunk& rhs)
    {
        os << indent::tab << "(Read_Chunk &=" << (void*)&rhs
           << indent::inc << indent::nl << "re_cptr=" << (void*)rhs._re_ptr
           << ",ce_cptr=" << (void*)rhs._ce_ptr
           << ",is_unmappable=" << (int)rhs._is_unmappable
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
