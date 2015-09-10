//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Read_Chunk_Pos.hpp"

#include "Read_Entry.hpp"
#include "Contig_Entry.hpp"


namespace MAC
{

/// Check if position is past the last mutation in the read chunk.
bool Read_Chunk_Pos::past_last_mut() const
{
    ASSERT(rc_cptr);
    return mca_cit == rc_cptr->mut_ptr_cont().end();
}

/// Check if position is past first mutation in the read chunk.
bool Read_Chunk_Pos::past_first_mut() const
{
    ASSERT(rc_cptr);
    return mca_cit != rc_cptr->mut_ptr_cont().begin();
}

void Read_Chunk_Pos::check() const
{
#ifndef DISABLE_ASSERTS
    ASSERT(rc_cptr != NULL);
    ASSERT(rc_cptr->get_c_start() <= c_pos and c_pos <= rc_cptr->get_c_end());
    ASSERT(rc_cptr->get_r_start() <= r_pos and r_pos <= rc_cptr->get_r_end());
    // if offset == 0; c_pos <= mut_start (or c_end if no more mutations)
    ASSERT(not (mut_offset == 0)
           or c_pos <= (not past_last_mut()?
                        mut().rf_start()
                        : rc_cptr->get_c_end()));
    // if offset != 0; c_pos == mut_start + min(mut_offset, mut_len)
    ASSERT(not (mut_offset != 0)
           or (not past_last_mut()
               and c_pos == mut().rf_start() + min(mut_offset, mut().rf_len())));
#endif
}

Size_Type Read_Chunk_Pos::get_match_len(bool forward) const
{
    check();
    if (mut_offset != 0)
    {
        return 0;
    }
    if (forward)
    {
        return (not past_last_mut()? mut().rf_start() : rc_cptr->get_c_end()) - c_pos;
    }
    else
    {
        return c_pos - (past_first_mut()? prev_mut().rf_end() : rc_cptr->get_c_start());
    }
}

void Read_Chunk_Pos::increment(Size_Type brk, bool on_contig)
{
    check();
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
            and c_pos == mut().rf_start() + min(mut_offset, mut().rf_len()))
    {
        // at the start or inside a mutation
        Size_Type c_leftover = (mut().rf_len() > mut_offset?
                                mut().rf_len() - mut_offset
                                : 0);
        Size_Type r_leftover = (mut().seq_len() > mut_offset?
                                mut().seq_len() - mut_offset
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
        if (mut_offset >= mut().rf_len() and mut_offset >= mut().seq_len())
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
    check();
    ASSERT(*this != rc_cptr->get_start_pos());
    if (brk == 0)
    {
        on_contig = true;
    }
    ASSERT(on_contig? (brk == 0 || brk < c_pos) : (not rc_cptr->get_rc()? brk < r_pos : r_pos < brk));
    if (mut_offset != 0
            or (mut_offset == 0 and past_first_mut()
                and c_pos == prev_mut().rf_end()))
    {
        // at the end or inside a mutation
        if (mut_offset == 0)
        {
            --mca_cit;
            mut_offset = max(mut().rf_len(), mut().seq_len());
        }
        ASSERT(mut_offset > 0);
        Size_Type c_leftover = min(mut().rf_len(), mut_offset);
        Size_Type r_leftover = min(mut().seq_len(), mut_offset);
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

Read_Chunk_Pos& Read_Chunk_Pos::jump_to_brk(Size_Type brk, bool on_contig)
{
    check();
    ASSERT(not on_contig or (rc_cptr->get_c_start() <= brk and brk <= rc_cptr->get_c_end()));
    ASSERT(on_contig or (rc_cptr->get_r_start() <= brk and brk <= rc_cptr->get_r_end()));
    bool forward = (on_contig?
                    c_pos <= brk
                    : (not rc_cptr->get_rc()) == (r_pos <= brk));
    while ((on_contig and c_pos != brk) or (not on_contig and r_pos != brk))
    {
        advance(forward, brk, on_contig);
    }
    return *this;
}

bool Read_Chunk_Pos::advance_til_mut(const Mutation& mut, bool forward)
{
    Size_Type mut_first;
    Size_Type mut_second;
    if (forward == not rc_cptr->get_rc())
    {
        mut_first = mut.rf_start();
        mut_second = mut.rf_end();
        ASSERT(r_pos <= mut_first);
    }
    else
    {
        mut_first = mut.rf_end();
        mut_second = mut.rf_start();
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

ptree Read_Chunk_Pos::to_ptree() const
{
    return ptree().put("c_pos", c_pos)
                  .put("r_pos", r_pos)
                  .put("mut_offset", mut_offset)
                  .put("mut", mut().to_ptree());
}

} // namespace MAC
