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

bool operator == (const Read_Chunk_Pos& lhs, const Read_Chunk_Pos& rhs)
{
    return (lhs.c_pos == rhs.c_pos
            and lhs.r_pos == rhs.r_pos
            and lhs.mca_cit == rhs.mca_cit
            and lhs.mut_offset == rhs.mut_offset);
}
bool operator != (const Read_Chunk_Pos& lhs, const Read_Chunk_Pos& rhs)
{
    return !(lhs == rhs);
}

/** Check if position is past the last mutation in the read chunk. */
bool Read_Chunk_Pos::past_last_mut() const
{
    ASSERT(rc_cptr);
    return mca_cit == rc_cptr->mut_ptr_cont().end();
}

/** Check if position is past first mutation in the read chunk. */
bool Read_Chunk_Pos::past_first_mut() const
{
    ASSERT(rc_cptr);
    return mca_cit != rc_cptr->mut_ptr_cont().begin();
}

bool Read_Chunk_Pos::check() const
{
    ASSERT(rc_cptr != NULL);
    ASSERT(rc_cptr->get_c_start() <= c_pos and c_pos <= rc_cptr->get_c_end());
    ASSERT(rc_cptr->get_r_start() <= r_pos and r_pos <= rc_cptr->get_r_end());
    // if offset == 0; c_pos <= mut_start (or c_end if no more mutations)
    ASSERT(not (mut_offset == 0)
           or c_pos <= (not past_last_mut()?
                        mut().get_start()
                        : rc_cptr->get_c_end()));
    // if offset != 0; c_pos == mut_start + min(mut_offset, mut_len)
    ASSERT(not (mut_offset != 0)
           or (not past_last_mut()
               and c_pos == mut().get_start() + min(mut_offset, mut().get_len())));
    return true;
}

inline Size_Type Read_Chunk_Pos::get_match_len(bool forward) const
{
    ASSERT(check());
    if (mut_offset != 0)
    {
        return 0;
    }
    if (forward)
    {
        return (not past_last_mut()? mut().get_start() : rc_cptr->get_c_end()) - c_pos;
    }
    else
    {
        return c_pos - (past_first_mut()? prev_mut().get_end() : rc_cptr->get_c_start());
    }
}

void Read_Chunk_Pos::increment(Size_Type brk, bool on_contig)
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
    if (not past_last_mut()
            and c_pos == mut().get_start() + min(mut_offset, mut().get_len()))
    {
        // at the start or inside a mutation
        Size_Type c_leftover = (mut().get_len() > mut_offset?
                                mut().get_len() - mut_offset
                                : 0);
        Size_Type r_leftover = (mut().get_seq_len() > mut_offset?
                                mut().get_seq_len() - mut_offset
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
        if (mut_offset >= mut().get_len() and mut_offset >= mut().get_seq_len())
        {
            mut_offset = 0;
            ++mca_cit;
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

void Read_Chunk_Pos::decrement(Size_Type brk, bool on_contig)
{
    ASSERT(check());
    ASSERT(*this != rc_cptr->get_start_pos());
    if (brk == 0)
    {
        on_contig = true;
    }
    ASSERT(on_contig? brk <= c_pos : (not rc_cptr->get_rc()? brk < r_pos : r_pos < brk));
    if (mut_offset != 0
            or (mut_offset == 0 and past_first_mut()
                and c_pos == prev_mut().get_end()))
    {
        // at the end or inside a mutation
        if (mut_offset == 0)
        {
            --mca_cit;
            mut_offset = max(mut().get_len(), mut().get_seq_len());
        }
        ASSERT(mut_offset > 0);
        Size_Type c_leftover = min(mut().get_len(), mut_offset);
        Size_Type r_leftover = min(mut().get_seq_len(), mut_offset);
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

void Read_Chunk_Pos::jump_to_brk(Size_Type brk, bool on_contig)
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

bool Read_Chunk_Pos::advance_til_mut(const Mutation& mut, bool forward)
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

void Read_Chunk_Pos::advance_past_del(bool forward)
{
    Read_Chunk_Pos tmp_pos = *this;
    while ((forward and *this != rc_cptr->get_end_pos()) or (not forward and *this != rc_cptr->get_start_pos()))
    {
        tmp_pos.advance(forward);
        if (tmp_pos.r_pos != r_pos)
        {
            break;
        }
        *this = tmp_pos;
    }
}

Size_Type Read_Chunk::get_read_len() const
{
    ASSERT(re_bptr());
    return re_bptr()->get_len();
}

Read_Chunk_BPtr Read_Chunk::make_chunk_from_cigar(const Cigar& cigar, Seq_Type&& rf, const Seq_Type& qr)
{
    ASSERT(rf.size() == cigar.get_rf_len() or rf.size() >= cigar.get_rf_start() + cigar.get_rf_len());
    ASSERT(qr.size() == cigar.get_qr_len() or qr.size() >= cigar.get_qr_start() + cigar.get_qr_len());

    // create objects with default constructor
    //shared_ptr< Read_Chunk > chunk_sptr(new Read_Chunk());
    //shared_ptr< Contig_Entry > ce_sptr(new Contig_Entry(rf_ptr, (rf_ptr->size() == cigar.get_rf_len()? cigar.get_rf_start() : 0)));
    Read_Chunk_BPtr rc_bptr = Read_Chunk_Fact::new_elem();
    Contig_Entry_BPtr ce_bptr = Contig_Entry_Fact::new_elem(std::move(rf),
                                                            (rf.size() == cigar.get_rf_len()? cigar.get_rf_start() : 0));

    // fix lengths and rc flags
    rc_bptr->_c_start = cigar.get_rf_start();
    rc_bptr->_c_len = cigar.get_rf_len();
    rc_bptr->_r_start = cigar.get_qr_start();
    rc_bptr->_r_len = cigar.get_qr_len();
    rc_bptr->_rc = cigar.is_reversed();

    // fix cross-pointers
    rc_bptr->_ce_bptr = ce_bptr;
    ce_bptr->chunk_cont().insert(rc_bptr);

    // construct mutations and store them in Contig_Entry object
    ce_bptr->mut_cont() = Mutation_Cont(cigar, qr);

    // store pointers in the Read_Chunk object
    rc_bptr->_mut_ptr_cont = Mutation_Ptr_Cont(ce_bptr->mut_cont(), rc_bptr);

    return rc_bptr;
}

/*
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
*/

std::tuple< Read_Chunk_BPtr, Read_Chunk_BPtr >
Read_Chunk::split(Read_Chunk_BPtr rc_bptr, Size_Type c_brk, Mutation_CBPtr mut_left_cbptr)
{
    ASSERT(rc_bptr->is_unlinked());
    ASSERT(mut_left_cbptr == nullptr or (mut_left_cbptr->get_start() == c_brk and mut_left_cbptr->is_ins()));
    Read_Chunk_BPtr left_rc_bptr = nullptr;
    Read_Chunk_BPtr right_rc_bptr = nullptr;

    // chunk stays intact on the lhs of the cut if:
    if (// endpoint is before c_brk
        rc_bptr->get_c_end() < c_brk
        // or at c_brk, and:
        or (rc_bptr->get_c_end() == c_brk
            and (// there are no mutations
                 rc_bptr->mut_ptr_cont().size() == 0
                 // or last mutation is not an insertion
                 or not rc_bptr->mut_ptr_cont().rbegin()->mut_cbptr()->is_ins()
                 // or last insertion is not at c_brk
                 or rc_bptr->mut_ptr_cont().rbegin()->mut_cbptr()->get_start() < c_brk
                 // or last insertion at c_brk is selected to stay on lhs
                 or rc_bptr->mut_ptr_cont().rbegin()->mut_cbptr() == mut_left_cbptr)))
    {
        return std::make_tuple(rc_bptr, right_rc_bptr);
    }

    // chunk stays intact on the rhs of the cut if:
    else if (// startpoint is after c_brk
             c_brk < rc_bptr->get_c_start()
             // or at c_brk, and:
             or (c_brk == rc_bptr->get_c_start()
                 and (// there are no mutations
                      rc_bptr->mut_ptr_cont().size() == 0
                      // or first mutation is not the insertion selected to stay on lhs
                      or rc_bptr->mut_ptr_cont().begin()->mut_cbptr() != mut_left_cbptr)))
    {
        return std::make_tuple(left_rc_bptr, rc_bptr);
    }

    else
    {
        // the chunk gets altered
        ASSERT(rc_bptr->get_c_start() <= c_brk and c_brk <= rc_bptr->get_c_end());
        // first compute the cut position
        Pos pos = rc_bptr->get_start_pos();
        pos.jump_to_brk(c_brk, true);
        // we are guaranteed no mutations span the break
        ASSERT(pos.c_pos == c_brk);
        ASSERT(pos.mut_offset == 0);
        // advance past an insertion at c_brk that is selected to stay on the lhs
        if (not pos.past_last_mut() and pos.mca_cit->mut_cbptr() == mut_left_cbptr)
        {
            pos.increment();
        }

        if (pos.r_pos == rc_bptr->get_r_start() or pos.r_pos == rc_bptr->get_r_end())
        {
            // the chunk stays in one piece
            if ((not rc_bptr->get_rc() and pos.r_pos == rc_bptr->get_r_start())
                or (rc_bptr->get_rc() and pos.r_pos == rc_bptr->get_r_end()))
            {
                // chunk goes on the rhs
                ASSERT((rc_bptr->get_c_start() == c_brk
                        and pos.mca_cit == rc_bptr->mut_ptr_cont().begin())
                       or (rc_bptr->get_c_start() < c_brk
                           and pos.mca_cit == ++(rc_bptr->mut_ptr_cont().begin())
                           and rc_bptr->mut_ptr_cont().begin()->mut_cbptr()->is_del()
                           and rc_bptr->mut_ptr_cont().begin()->mut_cbptr()->get_start() == rc_bptr->get_c_start()));
                // remove initial deletion, if any
                if (rc_bptr->get_c_start() < c_brk)
                {
                    Mutation_Chunk_Adapter_BPtr mca_bptr = &*rc_bptr->mut_ptr_cont().begin();
                    Mutation_BPtr mut_bptr = mca_bptr->mut_cbptr().unconst();
                    rc_bptr->mut_ptr_cont().erase(mca_bptr);
                    mut_bptr->chunk_ptr_cont().erase(mca_bptr);
                    Mutation_Chunk_Adapter_Fact::del_elem(mca_bptr);
                }
                return std::make_tuple(left_rc_bptr, rc_bptr);
            }
            else
            {
                // chunk goes on the lhs
                ASSERT((c_brk == rc_bptr->get_c_end()
                        and pos.mca_cit == rc_bptr->mut_ptr_cont().end())
                       or (c_brk < rc_bptr->get_c_end()
                           and pos.mca_cit == --(rc_bptr->mut_ptr_cont().end())
                           and pos.mca_cit->mut_cbptr()->is_del()
                           and pos.mca_cit->mut_cbptr()->get_end() == rc_bptr->get_c_end()));
                // remove final deletion, if any
                if (c_brk < rc_bptr->get_c_end())
                {
                    Mutation_Chunk_Adapter_BPtr mca_bptr = &*rc_bptr->mut_ptr_cont().rbegin();
                    Mutation_BPtr mut_bptr = mca_bptr->mut_cbptr().unconst();
                    rc_bptr->mut_ptr_cont().erase(mca_bptr);
                    mut_bptr->chunk_ptr_cont().erase(mca_bptr);
                    Mutation_Chunk_Adapter_Fact::del_elem(mca_bptr);
                }
                return std::make_tuple(rc_bptr, right_rc_bptr);
            }
        }
        else
        {
            // chunk gets cut
            ASSERT(rc_bptr->get_r_start() < pos.r_pos and pos.r_pos < rc_bptr->get_r_end());
            right_rc_bptr = Read_Chunk_Fact::new_elem();
            // fix inner fields of both chunks
            right_rc_bptr->_rc = rc_bptr->_rc;
            right_rc_bptr->_re_bptr = rc_bptr->_re_bptr;
            right_rc_bptr->_ce_bptr = rc_bptr->_ce_bptr;
            right_rc_bptr->_c_start = 0;
            right_rc_bptr->_c_len = rc_bptr->_c_start + rc_bptr->_c_len - c_brk;
            rc_bptr->_c_len -= right_rc_bptr->_c_len;
            if (not rc_bptr->_rc)
            {
                right_rc_bptr->_r_start = pos.r_pos;
                right_rc_bptr->_r_len = rc_bptr->_r_start + rc_bptr->_r_len - pos.r_pos;
            }
            else
            {
                right_rc_bptr->_r_start = rc_bptr->_r_start;
                right_rc_bptr->_r_len = pos.r_pos - rc_bptr->_r_start;
                rc_bptr->_r_start = pos.r_pos;
            }
            rc_bptr->_r_len -= right_rc_bptr->_r_len;
            // transfer mutations beyond breakpoint from left to right chunk
            right_rc_bptr->mut_ptr_cont() = rc_bptr->mut_ptr_cont().split(pos.mca_cit, right_rc_bptr);
            return std::make_tuple(rc_bptr, right_rc_bptr);
        }
    }
}

/*
std::tuple< bool, shared_ptr< Read_Chunk > > Read_Chunk::split(
    Size_Type c_brk, const map< const Mutation*, const Mutation* >& mut_cptr_map, const Contig_Entry* ce_cptr)
{
    // compute read chunk position corresponding to a contig cut at c_brk

    else
    {
        // read chunk gets split in 2
        ASSERT(get_r_start() < pos.r_pos and pos.r_pos < get_r_end());
        ASSERT(get_c_start() <= pos.c_pos and pos.c_pos <= get_c_end());

        // create new read chunk for second part of the contig
        shared_ptr< Read_Chunk > rc_sptr(new Read_Chunk());

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
        size_t r_mut_idx = (not get_rc()? r_mut_cnt : rc2.mut_ptr_cont().size() - 1 - r_mut_cnt); // next new mutation to look for
        const Mutation& r_mut = (r_mut_cnt == 0 and not past_start?
                                 fake_mut_start :
                                 r_mut_cnt < rc2.mut_ptr_cont().size()? *rc2.mut_ptr_cont()[r_mut_idx] : fake_mut_end);

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
                    Pos tmp_pos = pos_next;
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
                    Pos tmp_pos = pos_next;
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
            or (got_r_mut and r_mut_cnt == rc2.mut_ptr_cont().size()))
        {
            ASSERT(not got_r_mut or pos_next == pos);
            if (got_r_mut and pos_next.r_pos == get_end_pos().r_pos)
            {
                // if rc2 ends on contig end, incorporate remaining deletions
                Pos tmp_pos = pos_next;
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
            m.simplify(ce_bptr()->substr(m.get_start(), m.get_len()));
            if (not m.is_empty())
            {
                / *
                if (r_muts.size() == 0)
                {
                    // no adjacent contig mutations
                    ASSERT(c_muts.size() == 1);
                    res->mut_ptr_cont().push_back(c_muts[0]);
                }
                else
                * /
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
*/

void Read_Chunk::reverse()
{
    _rc = not _rc;
    _c_start = _ce_bptr->get_len() - (_c_start + _c_len);
}

/*
void Read_Chunk::merge_next(Read_Chunk_CPtr rc_next_cptr, Mutation::add_mut_mod_type add_mut_mod)
{
    ASSERT(rc_next_cptr != NULL);
    ASSERT(rc_next_cptr->get_re_ptr() == get_re_ptr());
    ASSERT(rc_next_cptr->ce_bptr() == ce_bptr());
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
*/

Seq_Type Read_Chunk::get_seq() const
{
    Seq_Type res;
    // start with start_pos iff not rc
    Pos pos = (not _rc? get_start_pos() : get_end_pos());
    while (pos != (not _rc? get_end_pos() : get_start_pos()))
    {
        Pos next_pos = pos;
        next_pos.advance(not _rc);
        Size_Type match_len = pos.get_match_len(not _rc);
        if (match_len > 0)
        {
            // match stretch follows
            string tmp = _ce_bptr->substr((not _rc? pos.c_pos : next_pos.c_pos) - _ce_bptr->get_seq_offset(), match_len);
            res += (not _rc? tmp : reverseComplement(tmp));
        }
        else if (pos.r_pos != next_pos.r_pos)
        {
            // mutation follows
            if (not _rc)
            {
                ASSERT(not pos.past_last_mut());
                res += pos.mut().get_seq().substr(pos.mut_offset, next_pos.r_pos - pos.r_pos);
            }
            else
            {
                res += reverseComplement(pos.mut().get_seq().substr(next_pos.mut_offset, pos.r_pos - next_pos.r_pos));
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
    ASSERT(get_c_start() <= ce_bptr()->get_seq_offset() + ce_bptr()->get_len());
    ASSERT(get_c_end() <= ce_bptr()->get_seq_offset() + ce_bptr()->get_len());
    // mapped length
    Size_Type c_len = get_c_end() - get_c_start();
    Size_Type r_len = get_r_end() - get_r_start();
    long long delta = 0;
    Mutation_CBPtr last_mut_cbptr = nullptr;
    for (const auto& mca_cbref : _mut_ptr_cont)
    {
        Mutation_CBPtr mut_cbptr = mca_cbref.raw().mut_cbptr();
        // no empty mutations
        ASSERT(not mut_cbptr->is_empty());
        // mutations must be in contig order
        //ASSERT(not last_mut_cbptr or last_mut_cbptr->get_end() <= mut_cbptr->get_start());
#ifndef ALLOW_CONSECUTIVE_MUTATIONS
        ASSERT(not last_mut_cbptr or last_mut_cbptr->get_end() < mut_cbptr->get_start());
#endif
        delta += (long long)mut_cbptr->get_seq_len() - (long long)mut_cbptr->get_len();
    }
    ASSERT((long long)c_len + delta == (long long)r_len);
#ifndef ALLOW_PART_MAPPED_INNER_CHUNKS
    // chunks must end on contig breaks except for first and last
    ASSERT(get_r_start() == 0
           or (not get_rc()?
               get_c_start() == 0
               : get_c_end() == ce_bptr()->get_seq_offset() + ce_bptr()->get_len()));
    ASSERT(get_r_end() == re_bptr()->get_len()
           or (not get_rc()?
               get_c_end() == ce_bptr()->get_seq_offset() + ce_bptr()->get_len()
               : get_c_start() == 0));
#endif
    // unmappable contigs have no mutations and a single chunk
    if (_is_unmappable)
    {
        ASSERT(_mut_ptr_cont.size() == 0);
        ASSERT(_ce_bptr->chunk_cont().size() == 1);
        ASSERT(&*(_ce_bptr->chunk_cont().begin()) == this);
        ASSERT(_ce_bptr->mut_cont().size() == 0);
    }
    return true;
}

ostream& operator << (ostream& os, const Read_Chunk_Pos& pos)
{
    os << "(c_pos=" << (size_t)pos.c_pos << ",r_pos=" << (size_t)pos.r_pos
       << ",mut=" << pos.mut() << ",mut_offset=" << (size_t)pos.mut_offset << ")";
    return os;
}

ostream& operator << (ostream& os, const Read_Chunk& rhs)
{
    os << indent::tab << "(Read_Chunk &=" << (void*)&rhs
       << indent::inc << indent::nl << "re_bptr=" << (void*)rhs._re_bptr
       << ",ce_bptr=" << (void*)rhs._ce_bptr
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

} // namespace MAC
