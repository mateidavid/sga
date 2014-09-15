//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Read_Chunk.hpp"

#include "Mutation.hpp"
#include "Read_Entry.hpp"
#include "Contig_Entry.hpp"
#include "../Util/Util.h"
#include "logger.hpp"

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
    ASSERT(on_contig? (brk == 0 || brk < c_pos) : (not rc_cptr->get_rc()? brk < r_pos : r_pos < brk));
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

Read_Chunk_Pos& Read_Chunk_Pos::jump_to_brk(Size_Type brk, bool on_contig)
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
    return *this;
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

Read_Chunk::Read_Chunk(Read_Entry_BPtr re_bptr, Contig_Entry_BPtr ce_bptr)
    : _re_bptr(re_bptr),
      _ce_bptr(ce_bptr),
      _r_start(0),
      _r_len(re_bptr->get_len()),
      _c_start(0),
      _c_len(ce_bptr->get_len()),
      _rc(false),
      _is_unmappable(false)
{
    ASSERT(_r_len == _c_len);
}

Size_Type Read_Chunk::get_read_len() const
{
    ASSERT(re_bptr());
    return re_bptr()->get_len();
}

Read_Chunk_BPtr
Read_Chunk::make_chunk_from_cigar(const Cigar& cigar, Seq_Type&& rf, const Seq_Type& qr)
{
    ASSERT(rf.size() == cigar.get_rf_len()
           or rf.size() >= cigar.get_rf_start() + cigar.get_rf_len());
    ASSERT(qr.size() == cigar.get_qr_len()
           or qr.size() >= cigar.get_qr_start() + cigar.get_qr_len());

    // create objects with default constructor
    //shared_ptr< Read_Chunk > chunk_sptr(new Read_Chunk());
    //shared_ptr< Contig_Entry > ce_sptr(new Contig_Entry(rf_ptr, (rf_ptr->size() == cigar.get_rf_len()? cigar.get_rf_start() : 0)));
    Read_Chunk_BPtr rc_bptr = Read_Chunk_Fact::new_elem();
    Contig_Entry_BPtr ce_bptr =
        Contig_Entry_Fact::new_elem(std::move(rf),
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

Read_Chunk_BPtr
Read_Chunk::make_relative_chunk(Read_Chunk_CBPtr rc1_cbptr,
                                Read_Chunk_CBPtr rc2_cbptr,
                                const Cigar& cigar)
{
    Read_Chunk_BPtr new_rc_bptr;
    new_rc_bptr = make_chunk_from_cigar(cigar, Seq_Type(rc1_cbptr->get_seq()), rc2_cbptr->get_seq());
    // fix Read_Entry pointer
    new_rc_bptr->re_bptr() = rc2_cbptr->re_bptr();
    return new_rc_bptr;
}

std::tuple< Read_Chunk_BPtr, Read_Chunk_BPtr >
Read_Chunk::split(Read_Chunk_BPtr rc_bptr, Size_Type c_brk, Mutation_CBPtr mut_left_cbptr)
{
    logger("Read_Chunk", debug1) << ptree("split")
        .put("rc", rc_bptr->to_ptree())
        .put("c_brk", c_brk)
        .put("mut_left_ptr", mut_left_cbptr.to_ptree());

    ASSERT(rc_bptr->is_unlinked());
    ASSERT(not mut_left_cbptr or (mut_left_cbptr->get_start() == c_brk and mut_left_cbptr->is_ins()));
    Read_Chunk_BPtr left_rc_bptr = nullptr;
    Read_Chunk_BPtr right_rc_bptr = nullptr;

    // chunk stays intact on the lhs of the cut if:
    if (// endpoint is before c_brk
        rc_bptr->get_c_end() < c_brk
        // or at c_brk, and:
        or (rc_bptr->get_c_end() == c_brk
            and (// there are no mutations
                 rc_bptr->mut_ptr_cont().empty()
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
                      rc_bptr->mut_ptr_cont().empty()
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
                           and not rc_bptr->mut_ptr_cont().empty()
                           and pos.mca_cit == ++(rc_bptr->mut_ptr_cont().begin())
                           and rc_bptr->mut_ptr_cont().begin()->mut_cbptr()->is_del()
                           and rc_bptr->mut_ptr_cont().begin()->mut_cbptr()->get_start() == rc_bptr->get_c_start()
                           and rc_bptr->mut_ptr_cont().begin()->mut_cbptr()->get_end() == c_brk));
                if (rc_bptr->get_c_start() < c_brk)
                {
                    // remove initial deletion and adjust contig coordinates
                    Mutation_Chunk_Adapter_BPtr mca_bptr = &*rc_bptr->mut_ptr_cont().begin();
                    Mutation_BPtr mut_bptr = mca_bptr->mut_cbptr().unconst();
                    rc_bptr->mut_ptr_cont().erase(mca_bptr);
                    mut_bptr->chunk_ptr_cont().erase(mca_bptr);
                    Mutation_Chunk_Adapter_Fact::del_elem(mca_bptr);
                    rc_bptr->_c_start = c_brk;
                    rc_bptr->_c_len -= mut_bptr->get_len();
                }
                return std::make_tuple(left_rc_bptr, rc_bptr);
            }
            else
            {
                // chunk goes on the lhs
                ASSERT((c_brk == rc_bptr->get_c_end()
                        and pos.mca_cit == rc_bptr->mut_ptr_cont().end())
                       or (c_brk < rc_bptr->get_c_end()
                           and not rc_bptr->mut_ptr_cont().empty()
                           and pos.mca_cit == --(rc_bptr->mut_ptr_cont().end())
                           and pos.mca_cit->mut_cbptr()->is_del()
                           and pos.mca_cit->mut_cbptr()->get_end() == rc_bptr->get_c_end()));
                if (c_brk < rc_bptr->get_c_end())
                {
                    // remove final deletion and adjust contig coordinates
                    Mutation_Chunk_Adapter_BPtr mca_bptr = &*rc_bptr->mut_ptr_cont().rbegin();
                    Mutation_BPtr mut_bptr = mca_bptr->mut_cbptr().unconst();
                    rc_bptr->mut_ptr_cont().erase(mca_bptr);
                    mut_bptr->chunk_ptr_cont().erase(mca_bptr);
                    Mutation_Chunk_Adapter_Fact::del_elem(mca_bptr);
                    rc_bptr->_c_len -= mut_bptr->get_len();
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
            right_rc_bptr->_c_start = c_brk;
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

Read_Chunk_BPtr
Read_Chunk::collapse_mapping(Read_Chunk_BPtr c1rc1_cbptr, Read_Chunk_BPtr rc1rc2_cbptr,
                             Mutation_Cont& mut_cont)
{
    // the entire part of rc1 to which rc2 is mapped must be in turn mapped to c1
    ASSERT(c1rc1_cbptr->get_r_start() <= rc1rc2_cbptr->get_c_start()
           and rc1rc2_cbptr->get_c_end() <= c1rc1_cbptr->get_r_end());
    Read_Chunk_BPtr c1rc2_bptr = Read_Chunk_Fact::new_elem();
    c1rc2_bptr->_r_start = rc1rc2_cbptr->get_r_start();
    c1rc2_bptr->_r_len = rc1rc2_cbptr->get_r_len();
    c1rc2_bptr->_rc = (c1rc1_cbptr->get_rc() != rc1rc2_cbptr->get_rc());
    c1rc2_bptr->re_bptr() = rc1rc2_cbptr->re_bptr();
    c1rc2_bptr->ce_bptr() = c1rc1_cbptr->ce_bptr();
    c1rc2_bptr->_is_unmappable = false;

    // traverse c1 left-to-right; traverse rc1 left-to-right if not _rc, r-to-l ow;
    Pos pos = c1rc1_cbptr->get_start_pos();

    /// iterator used to traverse the rc1 mutations in c1 order;
    /// if rc1 in the same orientation as c1, it points to next rc1 mutation;
    /// if rc1 in opposite orientation as c1, it points to last mutation considered;
    Mutation_Ptr_Cont::iterator r_mut_it = (not c1rc1_cbptr->get_rc()?
                                            rc1rc2_cbptr->mut_ptr_cont().begin()
                                            : rc1rc2_cbptr->mut_ptr_cont().end());
    /// empty rc1 mutation at the endpoint mapped to c1 start
    Mutation_CBPtr fake_start_r_mut_cbptr = Mutation_Fact::new_elem(
        not c1rc1_cbptr->get_rc()? rc1rc2_cbptr->get_c_start() : rc1rc2_cbptr->get_c_end(), 0);
    /// empty rc1 mutation at the endpoint mapped to c1 end
    Mutation_CBPtr fake_end_r_mut_cbptr = Mutation_Fact::new_elem(
        not c1rc1_cbptr->get_rc()? rc1rc2_cbptr->get_c_end() : rc1rc2_cbptr->get_c_start(), 0);
    /// bool used to trigger special code the first time in the while loop
    bool past_start = false;

    // mutation accumulator
    Mutation_BPtr crt_mut_bptr = Mutation_Fact::new_elem();
    // currently unused; could track mutation decomposition
    vector< Mutation_CBPtr > c_muts;
    vector< Mutation_CBPtr > r_muts;
    while (true)
    {
        /// next rc1 mutation to look for
        Mutation_CBPtr next_r_mut_cbptr;
        if (not past_start)
        {
            next_r_mut_cbptr = fake_start_r_mut_cbptr;
        }
        else
        {
            if (not c1rc1_cbptr->get_rc())
            {
                next_r_mut_cbptr = (r_mut_it != rc1rc2_cbptr->mut_ptr_cont().end()?
                                    r_mut_it->mut_cbptr()
                                    : fake_end_r_mut_cbptr);
            }
            else // c1rc1_cbptr->get_rc() == true
            {
                auto tmp = r_mut_it;
                next_r_mut_cbptr = (r_mut_it != rc1rc2_cbptr->mut_ptr_cont().begin()?
                                    (--tmp)->mut_cbptr()
                                    : fake_end_r_mut_cbptr);
            }
        }
        /// next c1 position
        Pos pos_next = pos;
        /// true iff current stretch represents current read mutation
        bool got_r_mut = pos_next.advance_til_mut(*next_r_mut_cbptr);

        if (not past_start)
        {
            ASSERT(next_r_mut_cbptr == fake_start_r_mut_cbptr);
            // if c1rc1_cbptr->get_r_start < rc1rc2_cbptr->get_c_start
            // first iteration produces a chunk of rc1 to which rc2 is not mapped
            if (got_r_mut)
            {
                ASSERT(pos_next == pos);
                c1rc2_bptr->_c_start = pos.c_pos;
                past_start = true;
                if (pos == c1rc1_cbptr->get_start_pos())
                {
                    // rc2 starts on contig break; must incorporate initial deletions
                    Pos tmp_pos = pos_next;
                    while (tmp_pos != c1rc1_cbptr->get_end_pos())
                    {
                        tmp_pos.increment();
                        if (tmp_pos.r_pos != pos_next.r_pos)
                        {
                            break;
                        }
                        // we found a deletion at c1 start
                        // initial deletions are always whole (offset==0 in pos_next and tmp_pos)
                        ASSERT(tmp_pos.r_pos == pos_next.r_pos);
                        ASSERT(tmp_pos.mut_offset == 0);
                        ASSERT(pos_next.mut_offset == 0);
                        ASSERT(tmp_pos.mca_cit == std::next(pos_next.mca_cit));
                        Mutation_CBPtr c_mut_cbptr = pos_next.mca_cit->mut_cbptr();
                        crt_mut_bptr->extend(c_mut_cbptr);
                        c_muts.push_back(c_mut_cbptr);
                        pos_next = tmp_pos;
                    }
                }
                else
                {
                    // if rc2 doesn't start on contig break, skip initial deletions
                    Pos tmp_pos = pos_next;
                    tmp_pos.advance_past_del();
                    // FIX: notorious situation:
                    // do not skip deletions if we can reach the contig end on deletions only,
                    // and read is aligned on contig end
                    if (tmp_pos != c1rc1_cbptr->get_end_pos())
                    {
                        pos_next = tmp_pos;
                        c1rc2_bptr->_c_start = pos_next.c_pos;
                    }
                    else
                    {
                        // rc2 contains insertions only, and is mapped to the contig end
                        // in this case, we leave the deletion, as it is needed
                        // to make the c1rc2_bptr aligned on the contig end
                        ASSERT(rc1rc2_cbptr->get_c_start() == rc1rc2_cbptr->get_c_end());
                    }
                } // else (rc2 doesn't start on c1 break)
            } // if (got_r_mut)
        } // if (not past_start)
        else if ((not got_r_mut and pos.get_match_len() > 0)
            or (got_r_mut and next_r_mut_cbptr == fake_end_r_mut_cbptr))
        {
            // passed a match strech, or we hit the end of rc2
            ASSERT(not got_r_mut or pos_next == pos);
            if (got_r_mut and pos_next.r_pos == c1rc1_cbptr->get_end_pos().r_pos)
            {
                // rc2 ends on contig end; must incorporate remaining deletions
                Pos tmp_pos = pos_next;
                while (tmp_pos != c1rc1_cbptr->get_end_pos())
                {
                    tmp_pos.increment();
                    ASSERT(tmp_pos.r_pos == pos_next.r_pos);
                    ASSERT(tmp_pos.mut_offset == 0);
                    ASSERT(tmp_pos.mca_cit == std::next(pos_next.mca_cit));
                    // final deletions are whole, or we have used a read mutation
                    ASSERT(pos_next.mut_offset == 0 or r_muts.size() > 0);
                    Mutation_CBPtr c_mut_cbptr = pos_next.mca_cit->mut_cbptr();
                    //Mutation_CBPtr tmp_mut_cbptr = Mutation_Fact::new_elem(pos_next.c_pos, tmp_pos.c_pos - pos_next.c_pos, 0);
                    // extend with deletion
                    crt_mut_bptr->extend(pos_next.c_pos, tmp_pos.c_pos - pos_next.c_pos, Seq_Type());
                    c_muts.push_back(c_mut_cbptr);
                    pos_next = tmp_pos;
                }
            }
            // the stretch pos->pos_next is a match, or we are at the end;
            // consolidate any outstanding mutations
            crt_mut_bptr->simplify(c1rc1_cbptr->ce_bptr()->substr(crt_mut_bptr->get_start(), crt_mut_bptr->get_len()));
            if (not crt_mut_bptr->is_empty())
            {
                // use an equivalent Mutation if one exists; else add this one
                crt_mut_bptr = mut_cont.find_equiv_or_add(crt_mut_bptr).unconst();
                // add MCA
                Mutation_Chunk_Adapter_BPtr mca_bptr = Mutation_Chunk_Adapter_Fact::new_elem(crt_mut_bptr, c1rc2_bptr);
                crt_mut_bptr->chunk_ptr_cont().insert(mca_bptr);
                c1rc2_bptr->mut_ptr_cont().push_back(mca_bptr);
            }
            crt_mut_bptr = Mutation_Fact::new_elem();
            c_muts.clear();
            r_muts.clear();

            if (got_r_mut)
            {
                // exit point of main while loop
                c1rc2_bptr->_c_len = pos_next.c_pos - c1rc2_bptr->_c_start;
                break;
            }
        } // else if (match strech or fake_end mutation)
        else
        {
            // the stretch pos->post_next is not a match and we are not at the end;
            // add mutation slice to m
            if (got_r_mut)
            {
                // this was the entire read mutation next_r_mut_cbptr
                ASSERT(c1rc1_cbptr->get_rc() or (pos.r_pos == next_r_mut_cbptr->get_start()
                                                 and pos_next.r_pos == next_r_mut_cbptr->get_end()));
                ASSERT(not c1rc1_cbptr->get_rc() or (pos.r_pos == next_r_mut_cbptr->get_end()
                                                     and pos_next.r_pos == next_r_mut_cbptr->get_start()));
                crt_mut_bptr->extend(pos.c_pos, pos_next.c_pos - pos.c_pos,
                                     (not c1rc1_cbptr->get_rc()?
                                      next_r_mut_cbptr->get_seq()
                                      : reverseComplement(next_r_mut_cbptr->get_seq())));
                r_muts.push_back(next_r_mut_cbptr);
                // advance read mutation iterator
                if (not c1rc1_cbptr->get_rc())
                {
                    ++r_mut_it;
                }
                else
                {
                    --r_mut_it;
                }
            }
            else
            {
                // this was a (possibly sliced) contig mutation
                ASSERT((pos_next.mca_cit == pos.mca_cit and pos_next.mut_offset > pos.mut_offset)
                       or (pos_next.mca_cit == std::next(pos.mca_cit) and pos_next.mut_offset == 0));
                //const Mutation& c_mut = *_mut_ptr_cont[pos.mut_idx];
                Mutation_CBPtr c_mut_cbptr = pos.mca_cit->mut_cbptr();
                crt_mut_bptr->extend(pos.c_pos, pos_next.c_pos - pos.c_pos,
                                     c_mut_cbptr->get_seq().substr(
                                         min(pos.mut_offset, c_mut_cbptr->get_seq_len()),
                                         (pos_next.mut_offset == 0? string::npos : pos_next.mut_offset - pos.mut_offset)));
                c_muts.push_back(c_mut_cbptr);
            }
        }
        pos = pos_next;
    }
    // clean up temporary mutations
    Mutation_Fact::del_elem(fake_start_r_mut_cbptr);
    Mutation_Fact::del_elem(fake_end_r_mut_cbptr);
    Mutation_Fact::del_elem(crt_mut_bptr);

    // if rc2 is mapped to endpoints of rc1, rc2 will be mapped to the same extent of c1;
    // in particular, (in this case) no extremal deletions of c1 will be omitted
    ASSERT(not (rc1rc2_cbptr->get_c_start() == c1rc1_cbptr->get_r_start())
           or (not c1rc1_cbptr->get_rc()?
               c1rc2_bptr->get_c_start() == c1rc1_cbptr->get_c_start()
               : c1rc2_bptr->get_c_end() == c1rc1_cbptr->get_c_end()));
    ASSERT(not (rc1rc2_cbptr->get_c_end() == c1rc1_cbptr->get_r_end())
           or (not c1rc1_cbptr->get_rc()?
               c1rc2_bptr->get_c_end() == c1rc1_cbptr->get_c_end()
               : c1rc2_bptr->get_c_start() == c1rc1_cbptr->get_c_start()));
    return c1rc2_bptr;
}

Read_Chunk_BPtr Read_Chunk::invert_mapping(Read_Chunk_CBPtr rc_cbptr)
{
    // create new chunk and contig entry
    Read_Chunk_BPtr new_rc_bptr = Read_Chunk_Fact::new_elem();
    Contig_Entry_BPtr new_ce_bptr = Contig_Entry_Fact::new_elem(
        string(rc_cbptr->get_seq()), rc_cbptr->get_r_start());
    // basic settings
    new_rc_bptr->_r_start = rc_cbptr->get_c_start();
    new_rc_bptr->_r_len = rc_cbptr->get_c_len();
    new_rc_bptr->_c_start = rc_cbptr->get_r_start();
    new_rc_bptr->_c_len = rc_cbptr->get_r_len();
    new_rc_bptr->_rc = rc_cbptr->get_rc();
    new_rc_bptr->ce_bptr() = new_ce_bptr;
    // iterate over contig sequence and create reversed mutations
    Pos pos = rc_cbptr->get_start_pos();
    while (true)
    {
        ASSERT(pos != rc_cbptr->get_end_pos() or pos.get_match_len() == 0);
        // ignore matched stretches
        while (pos.get_match_len() > 0)
        {
            pos.increment();
        }
        if (pos == rc_cbptr->get_end_pos())
        {
            break;
        }
        // with no breakpoints, we don't slice mutations
        ASSERT(pos.mut_offset == 0);
        // the next stretch is a mutation
        Pos pos_next = pos;
        pos_next.increment();
        ASSERT(pos_next.mut_offset == 0);
        ASSERT(pos_next.mca_cit == std::next(pos.mca_cit));
        // create reversed Mutation
        /*
        Mutation rev_mut((not _rc? pos.r_pos : pos_next.r_pos),
                         (not _rc? pos_next.r_pos - pos.r_pos : pos.r_pos - pos_next.r_pos),
                         (not _rc? _ce_ptr->substr(pos.c_pos, pos_next.c_pos - pos.c_pos)
                          : reverseComplement(_ce_ptr->substr(pos.c_pos, pos_next.c_pos - pos.c_pos))));
        */
        string mut_c_seq = rc_cbptr->ce_bptr()->substr(pos.c_pos, pos_next.c_pos - pos.c_pos);
        Mutation_BPtr new_mut_bptr = (
            not rc_cbptr->get_rc()?
            Mutation_Fact::new_elem(pos.r_pos, pos_next.r_pos - pos.r_pos, mut_c_seq)
            : Mutation_Fact::new_elem(pos_next.r_pos, pos.r_pos - pos_next.r_pos, reverseComplement(mut_c_seq)));
        // save it in Mutation container
        new_ce_bptr->mut_cont().insert(new_mut_bptr);
        // add MCA
        Mutation_Chunk_Adapter_BPtr mca_bptr = Mutation_Chunk_Adapter_Fact::new_elem(new_mut_bptr, new_rc_bptr);
        new_mut_bptr->chunk_ptr_cont().insert(mca_bptr);
        new_rc_bptr->mut_ptr_cont().insert_before(
            not rc_cbptr->get_rc()? new_rc_bptr->mut_ptr_cont().end() : new_rc_bptr->mut_ptr_cont().begin(),
            mca_bptr);
        pos = pos_next;
    }
    // add chunk to contig
    new_ce_bptr->chunk_cont().insert(new_rc_bptr);
    return new_rc_bptr;
}

void Read_Chunk::reverse()
{
    _rc = not _rc;
    _c_start = _ce_bptr->get_len() - (_c_start + _c_len);
}

void Read_Chunk::cat_c_right(Read_Chunk_BPtr rc_bptr, Read_Chunk_BPtr rc_next_bptr,
                             Mutation_Cont& mut_cont)
{
    ASSERT(rc_next_bptr);
    ASSERT(rc_next_bptr->re_bptr() == rc_bptr->re_bptr());
    ASSERT(rc_next_bptr->ce_bptr() == rc_bptr->ce_bptr());
    ASSERT(rc_next_bptr->get_rc() == rc_bptr->get_rc());
    ASSERT(rc_bptr->get_rc() or rc_next_bptr->get_r_start() == rc_bptr->get_r_end());
    ASSERT(not rc_bptr->get_rc() or rc_next_bptr->get_r_end() == rc_bptr->get_r_start());
    ASSERT(rc_next_bptr->get_c_start() == rc_bptr->get_c_end());
    ASSERT(rc_bptr->is_unlinked());
    ASSERT(rc_next_bptr->is_unlinked());

    // fix coordinates
    if (rc_bptr->get_rc())
    {
        rc_bptr->_r_start = rc_next_bptr->_r_start;
    }
    rc_bptr->_c_len += rc_next_bptr->_c_len;
    rc_bptr->_r_len += rc_next_bptr->_r_len;
    // if there are touching mutations across the break, merge them
    if (not rc_bptr->mut_ptr_cont().empty()
        and not rc_next_bptr->mut_ptr_cont().empty()
        and rc_bptr->mut_ptr_cont().rbegin()->mut_cbptr()->get_end()
            == rc_next_bptr->mut_ptr_cont().begin()->mut_cbptr()->get_start())
    {
        // create new merged Mutation
        Mutation_BPtr left_mut_bptr = rc_bptr->mut_ptr_cont().rbegin()->mut_cbptr().unconst();
        Mutation_BPtr right_mut_bptr = rc_next_bptr->mut_ptr_cont().begin()->mut_cbptr().unconst();
        Mutation_BPtr new_mut_bptr = Mutation_Fact::new_elem(left_mut_bptr->get_start(),
                                                             left_mut_bptr->get_len(),
                                                             left_mut_bptr->get_seq());
        new_mut_bptr->extend(right_mut_bptr);
        Mutation_BPtr equiv_mut_bptr = mut_cont.find(new_mut_bptr, false).unconst();
        if (equiv_mut_bptr)
        {
            // an equivalent merged Mutation exists; use it
            Mutation_Fact::del_elem(new_mut_bptr);
            new_mut_bptr = equiv_mut_bptr;
        }
        else
        {
            // no equivalent Mutation exists; insert this one
            mut_cont.insert(new_mut_bptr);
        }
        // remove MCAs for left&right mutations
        Mutation_Chunk_Adapter_BPtr left_mca_bptr = &*rc_bptr->mut_ptr_cont().rbegin();
        Mutation_Chunk_Adapter_BPtr right_mca_bptr = &*rc_next_bptr->mut_ptr_cont().begin();
        left_mut_bptr->chunk_ptr_cont().erase(left_mca_bptr);
        right_mut_bptr->chunk_ptr_cont().erase(right_mca_bptr);
        rc_bptr->mut_ptr_cont().erase(left_mca_bptr);
        rc_next_bptr->mut_ptr_cont().erase(right_mca_bptr);
        // if the old Mutations are no longer used, erase&deallocate them
        if (left_mut_bptr->chunk_ptr_cont().empty())
        {
            mut_cont.erase(left_mut_bptr);
            Mutation_Fact::del_elem(left_mut_bptr);
        }
        if (right_mut_bptr->chunk_ptr_cont().empty())
        {
            mut_cont.erase(right_mut_bptr);
            Mutation_Fact::del_elem(right_mut_bptr);
        }
        // insert MCA for new Mutation
        Mutation_Chunk_Adapter_BPtr new_mca_bptr = Mutation_Chunk_Adapter_Fact::new_elem(new_mut_bptr, rc_bptr);
        rc_bptr->mut_ptr_cont().insert_before(rc_bptr->mut_ptr_cont().end(), new_mca_bptr);
        new_mut_bptr->chunk_ptr_cont().insert(new_mca_bptr);
    } // if [exist touching mutations]
    // grab rhs mutations
    rc_bptr->mut_ptr_cont().splice_right(rc_next_bptr->mut_ptr_cont(), rc_bptr);
    // deallocate rhs chunk
    Read_Chunk_Fact::del_elem(rc_next_bptr);
}

tuple< Size_Type, Size_Type >
Read_Chunk::mapped_range(Size_Type rg_start, Size_Type rg_end, bool on_contig,
                         bool rg_start_maximal, bool rg_end_maximal) const
{
    ASSERT(rg_start <= rg_end);
    ASSERT(not on_contig or (get_c_start() <= rg_start and rg_end <= get_c_end()));
    ASSERT(on_contig or (get_r_start() <= rg_start and rg_end <= get_r_end()));
    // compute range endpoints in contig order
    Size_Type rg_endpoint_c_left;
    Size_Type rg_endpoint_c_right;
    bool rg_endpoint_c_left_maximal;
    bool rg_endpoint_c_right_maximal;
    if (on_contig or not get_rc())
    {
        rg_endpoint_c_left = rg_start;
        rg_endpoint_c_right = rg_end;
        rg_endpoint_c_left_maximal = rg_start_maximal;
        rg_endpoint_c_right_maximal = rg_end_maximal;
    }
    else
    {
        rg_endpoint_c_left = rg_end;
        rg_endpoint_c_right = rg_start;
        rg_endpoint_c_left_maximal = rg_end_maximal;
        rg_endpoint_c_right_maximal = rg_start_maximal;
    }
    // compute map start position
    auto start_pos = get_start_pos().jump_to_brk(rg_endpoint_c_left, on_contig);
    if (not rg_endpoint_c_left_maximal and start_pos.get_match_len() == 0)
    {
        // incorporate mutation if range position doesn't change
        auto next_pos = start_pos;
        next_pos.increment();
        if (on_contig? next_pos.c_pos == start_pos.c_pos : next_pos.r_pos == start_pos.r_pos)
        {
            start_pos = next_pos;
        }
    }
    // compute map end position
    auto end_pos = get_end_pos().jump_to_brk(rg_endpoint_c_right, on_contig);
    if (not rg_endpoint_c_right_maximal and end_pos.get_match_len() == 0)
    {
        // incorporate mutation if range position doesn't change
        auto next_pos = end_pos;
        next_pos.decrement();
        if (on_contig? next_pos.c_pos == end_pos.c_pos : next_pos.r_pos == end_pos.r_pos)
        {
            end_pos = next_pos;
        }
    }
    Size_Type ret_rg_start;
    Size_Type ret_rg_end;
    if (on_contig)
    {
        if (not get_rc())
        {
            ret_rg_start = start_pos.r_pos;
            ret_rg_end = end_pos.r_pos;
        }
        else
        {
            ret_rg_start = end_pos.r_pos;
            ret_rg_end = start_pos.r_pos;
        }
    }
    else // not on_contig
    {
        ret_rg_start = start_pos.c_pos;
        ret_rg_end = end_pos.c_pos;
    }
    ASSERT(ret_rg_start <= ret_rg_end or (rg_start == rg_end and not rg_start_maximal and not rg_end_maximal));
    return std::make_tuple(ret_rg_start, ret_rg_end);
}

void Read_Chunk::make_unmappable(Read_Chunk_BPtr rc_bptr)
{
    Contig_Entry_BPtr ce_new_bptr = Contig_Entry_Fact::new_elem(rc_bptr->get_seq());
    ASSERT(ce_new_bptr->get_len() == rc_bptr->_r_len);
    rc_bptr->mut_ptr_cont().clear_and_dispose();
    rc_bptr->ce_bptr() = ce_new_bptr;
    rc_bptr->_c_start = 0;
    rc_bptr->_c_len = rc_bptr->_r_len;
    rc_bptr->_rc = false;
    rc_bptr->_is_unmappable = true;
    ce_new_bptr->chunk_cont().insert(rc_bptr);
}

Seq_Type Read_Chunk::get_seq() const
{
    Seq_Type res;
    // start with start_pos iff not rc
    Pos pos = (not _rc? get_start_pos() : get_end_pos());
    while (pos != (not _rc? get_end_pos() : get_start_pos()))
    {
        Pos next_pos = pos;
        next_pos.advance(not _rc);
        ASSERT(next_pos.mut_offset == 0);
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
                res += pos.mut().get_seq().substr(0, next_pos.r_pos - pos.r_pos);
            }
            else
            {
                res += reverseComplement(next_pos.mut().get_seq().substr(0, pos.r_pos - next_pos.r_pos));
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
        ASSERT((&*(_ce_bptr->chunk_cont().begin())).raw() == this);
        ASSERT(_ce_bptr->mut_cont().empty());
    }
    return true;
}

boost::property_tree::ptree Read_Chunk_Pos::to_ptree() const
{
    return ptree().put("c_pos", c_pos)
                  .put("r_pos", r_pos)
                  .put("mut_offset", mut_offset)
                  .put("mut", mut().to_ptree());
}

boost::property_tree::ptree Read_Chunk::to_ptree() const
{
    return ptree().put("re_bptr", re_bptr().to_ptree())
                  .put("ce_bptr", ce_bptr().to_ptree())
                  .put("is_unmappable", is_unmappable())
                  .put("r_start", get_r_start())
                  .put("r_len", get_r_len())
                  .put("c_start", get_c_start())
                  .put("c_len", get_c_len())
                  .put("rc", get_rc())
                  .put("mut_ptr_cont", cont_to_ptree(mut_ptr_cont()))
                  .put("re_parent", _re_parent.to_ptree())
                  .put("re_l_child", _re_l_child.to_ptree())
                  .put("re_r_child", _re_r_child.to_ptree())
                  .put("ce_parent", _ce_parent.to_ptree())
                  .put("ce_l_child", _ce_l_child.to_ptree())
                  .put("ce_r_child", _ce_r_child.to_ptree());
}

} // namespace MAC
