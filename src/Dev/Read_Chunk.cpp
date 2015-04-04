//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Read_Chunk.hpp"

#include "Read_Entry.hpp"
#include "Contig_Entry.hpp"


namespace MAC
{

Read_Chunk::Read_Chunk(Read_Entry_BPtr re_bptr, Contig_Entry_BPtr ce_bptr)
    : _r_start(0),
      _r_len(re_bptr->len()),
      _c_start(0),
      _c_len(ce_bptr->len()),
      _re_bptr(re_bptr),
      _ce_bptr(ce_bptr),
      _mut_ptr_cont(),
      _bits(0)
{
    ASSERT(_r_len == _c_len);
}

Size_Type Read_Chunk::get_read_len() const
{
    ASSERT(re_bptr());
    return re_bptr()->len();
}

Read_Chunk_BPtr
Read_Chunk::make_chunk_from_cigar(const Cigar& cigar, Seq_Type&& rf, const Seq_Proxy_Type& qr)
{
    ASSERT(rf.size() == cigar.rf_len()
           or rf.size() >= cigar.rf_start() + cigar.rf_len());
    ASSERT(qr.size() == cigar.qr_len()
           or qr.size() >= cigar.qr_start() + cigar.qr_len());

    Read_Chunk_BPtr rc_bptr = Read_Chunk_Fact::new_elem();
    Contig_Entry_BPtr ce_bptr =
        Contig_Entry_Fact::new_elem(move(rf),
                                    (rf.size() == cigar.rf_len()? cigar.rf_start() : 0));

    // fix lengths and rc flags
    rc_bptr->_c_start = cigar.rf_start();
    rc_bptr->_c_len = cigar.rf_len();
    rc_bptr->_r_start = cigar.qr_start();
    rc_bptr->_r_len = cigar.qr_len();
    rc_bptr->_set_rc(cigar.reversed());

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
Read_Chunk::make_chunk_from_cigar(const Cigar& cigar, const Seq_Proxy_Type& qr, Contig_Entry_BPtr ce_bptr)
{
    ASSERT(ce_bptr->len() >= cigar.rf_len());
    ASSERT(qr.size() == cigar.qr_len()
           or qr.size() >= cigar.qr_start() + cigar.qr_len());

    Read_Chunk_BPtr rc_bptr = Read_Chunk_Fact::new_elem();

    // fix lengths and rc flags
    rc_bptr->_c_start = cigar.rf_start();
    rc_bptr->_c_len = cigar.rf_len();
    rc_bptr->_r_start = cigar.qr_start();
    rc_bptr->_r_len = cigar.qr_len();
    rc_bptr->_set_rc(cigar.reversed());

    // fix cross-pointers
    rc_bptr->_ce_bptr = ce_bptr;
    ce_bptr->chunk_cont().insert(rc_bptr);

    // construct mutations
    Mutation_Cont mut_cont = Mutation_Cont(cigar, qr);
    // set mut_ptr_cont (this adds mca-s)
    rc_bptr->_mut_ptr_cont = Mutation_Ptr_Cont(mut_cont, rc_bptr);
    // merge mut_cont with mutation container in ce_bptr
    ce_bptr->mut_cont().merge(mut_cont);

    return rc_bptr;
}

Read_Chunk_BPtr
Read_Chunk::make_relative_chunk(Read_Chunk_CBPtr rc1_cbptr, Read_Chunk_CBPtr rc2_cbptr,
                                const Cigar& cigar)
{
    Read_Chunk_BPtr new_rc_bptr;
    new_rc_bptr = make_chunk_from_cigar(cigar, rc1_cbptr->get_seq(), rc2_cbptr->get_seq());
    // fix Read_Entry pointer
    new_rc_bptr->re_bptr() = rc2_cbptr->re_bptr();
    return new_rc_bptr;
}

pair< Read_Chunk_BPtr, Read_Chunk_BPtr >
Read_Chunk::split(Read_Chunk_BPtr rc_bptr, Size_Type c_brk, Mutation_CBPtr mut_left_cbptr, bool strict)
{
    LOG("Read_Chunk", debug1) << ptree("split")
        .put("rc", rc_bptr->to_ptree())
        .put("c_brk", c_brk)
        .put("mut_left_ptr", mut_left_cbptr.to_ptree())
        .put("strict", strict);

    ASSERT(not mut_left_cbptr
           or (mut_left_cbptr->rf_start() == c_brk and mut_left_cbptr->is_ins()));
    ASSERT(rc_bptr->is_unlinked());

/*
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
                or rc_bptr->mut_ptr_cont().rbegin()->mut_cbptr()->rf_start() < c_brk
                // or last insertion at c_brk is selected to stay on lhs
                or rc_bptr->mut_ptr_cont().rbegin()->mut_cbptr() == mut_left_cbptr)))
    {
        return make_pair(rc_bptr, right_rc_bptr);
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
        return make_pair(left_rc_bptr, rc_bptr);
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
                           and rc_bptr->mut_ptr_cont().begin()->mut_cbptr()->rf_start() == rc_bptr->get_c_start()
                           and rc_bptr->mut_ptr_cont().begin()->mut_cbptr()->rf_end() == c_brk));
                if (rc_bptr->get_c_start() < c_brk)
                {
                    // remove initial deletion and adjust contig coordinates
                    Mutation_Chunk_Adapter_BPtr mca_bptr = &*rc_bptr->mut_ptr_cont().begin();
                    Mutation_BPtr mut_bptr = mca_bptr->mut_cbptr().unconst();
                    rc_bptr->mut_ptr_cont().erase(mca_bptr);
                    mut_bptr->chunk_ptr_cont().erase(mca_bptr);
                    Mutation_Chunk_Adapter_Fact::del_elem(mca_bptr);
                    rc_bptr->_c_start = c_brk;
                    rc_bptr->_c_len -= mut_bptr->rf_len();
                }
                return make_pair(left_rc_bptr, rc_bptr);
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
                           and pos.mca_cit->mut_cbptr()->rf_end() == rc_bptr->get_c_end()));
                if (c_brk < rc_bptr->get_c_end())
                {
                    // remove final deletion and adjust contig coordinates
                    Mutation_Chunk_Adapter_BPtr mca_bptr = &*rc_bptr->mut_ptr_cont().rbegin();
                    Mutation_BPtr mut_bptr = mca_bptr->mut_cbptr().unconst();
                    rc_bptr->mut_ptr_cont().erase(mca_bptr);
                    mut_bptr->chunk_ptr_cont().erase(mca_bptr);
                    Mutation_Chunk_Adapter_Fact::del_elem(mca_bptr);
                    rc_bptr->_c_len -= mut_bptr->rf_len();
                }
                return make_pair(rc_bptr, right_rc_bptr);
            }
        }
        else
        {
            // chunk gets cut
            ASSERT(rc_bptr->get_r_start() < pos.r_pos and pos.r_pos < rc_bptr->get_r_end());
            right_rc_bptr = Read_Chunk_Fact::new_elem();
            // fix inner fields of both chunks
            right_rc_bptr->_set_rc(rc_bptr->get_rc());
            right_rc_bptr->_re_bptr = rc_bptr->_re_bptr;
            right_rc_bptr->_ce_bptr = rc_bptr->_ce_bptr;
            right_rc_bptr->_c_start = c_brk;
            right_rc_bptr->_c_len = rc_bptr->_c_start + rc_bptr->_c_len - c_brk;
            rc_bptr->_c_len -= right_rc_bptr->_c_len;
            if (not rc_bptr->get_rc())
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
            return make_pair(rc_bptr, right_rc_bptr);
        }
    }
*/

    if (rc_bptr->get_c_end() < c_brk)
    {
        return make_pair(rc_bptr, nullptr);
    }
    if (c_brk < rc_bptr->get_c_start())
    {
        return make_pair(nullptr, rc_bptr);
    }
    ASSERT(rc_bptr->get_c_start() <= c_brk and c_brk <= rc_bptr->get_c_end());
    // scan mutation_ptr container to find:
    // - (iterator to) the first mutation which goes on the rhs of the split
    // - the difference in read_len - contig_len on the lhs
    ptrdiff_t delta_len = 0;
    auto it = rc_bptr->mut_ptr_cont().begin();
    while (it != rc_bptr->mut_ptr_cont().end()
           and (it->mut_cbptr()->rf_start() < c_brk
                or it->mut_cbptr() == mut_left_cbptr))
    {
        delta_len += it->mut_cbptr()->seq_len();
        delta_len -= it->mut_cbptr()->rf_len();
        ++it;
    }
    // compute lengths of the lhs and rhs
    Size_Type lhs_c_len = c_brk - rc_bptr->get_c_start();
    Size_Type rhs_c_len = rc_bptr->get_c_end() - c_brk;
    ASSERT(lhs_c_len + rhs_c_len == rc_bptr->get_c_len());
    ASSERT(static_cast< ptrdiff_t >(lhs_c_len) + delta_len >= 0);
    Size_Type lhs_r_len = static_cast< ptrdiff_t >(lhs_c_len) + delta_len;
    ASSERT(lhs_r_len <= rc_bptr->get_r_len());
    Size_Type rhs_r_len = rc_bptr->get_r_len() - lhs_r_len;

    Read_Chunk_BPtr left_rc_bptr = nullptr;
    Read_Chunk_BPtr right_rc_bptr = nullptr;
    if (lhs_r_len > 0 or strict)
    {
        // fix lhs
        left_rc_bptr = rc_bptr;
        left_rc_bptr->_c_len = lhs_c_len;
        left_rc_bptr->_r_start += (not rc_bptr->get_rc()? 0 : rhs_r_len);
        left_rc_bptr->_r_len = lhs_r_len;
        if (rhs_r_len > 0 or strict)
        {
            // create new rhs
            right_rc_bptr = Read_Chunk_Fact::new_elem(
                (not rc_bptr->get_rc()? left_rc_bptr->get_r_end() : left_rc_bptr->get_r_start() - rhs_r_len), rhs_r_len,
                c_brk, rhs_c_len,
                rc_bptr->get_rc());
            right_rc_bptr->re_bptr() = rc_bptr->re_bptr();
            right_rc_bptr->ce_bptr() = rc_bptr->ce_bptr();
            // transfer mutations [it, end) to rhs
            right_rc_bptr->mut_ptr_cont().splice(left_rc_bptr->mut_ptr_cont(), right_rc_bptr, it);
        }
        else if (rhs_c_len > 0)
        {
            // rhs consists of a single deletion
            ASSERT(it == prev(rc_bptr->mut_ptr_cont().end()));
            rc_bptr->mut_ptr_cont().erase_and_dispose(&*it);
        }
    }
    else // lhs is empty: do not create 2 chunks, original chunk stays on the right
    {
        // fix rhs
        right_rc_bptr = rc_bptr;
        right_rc_bptr->_c_start = c_brk;
        right_rc_bptr->_c_len = rhs_c_len;
        right_rc_bptr->_r_start += (not rc_bptr->get_rc()? lhs_r_len : 0);
        right_rc_bptr->_r_len = rhs_r_len;
        if (lhs_c_len > 0)
        {
            // lhs consists of a single deletion
            ASSERT(it == next(rc_bptr->mut_ptr_cont().begin()));
            rc_bptr->mut_ptr_cont().erase_and_dispose(&*prev(it));
        }
    }
    return make_pair(left_rc_bptr, right_rc_bptr);
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
    c1rc2_bptr->_set_rc(c1rc1_cbptr->get_rc() != rc1rc2_cbptr->get_rc());
    c1rc2_bptr->re_bptr() = rc1rc2_cbptr->re_bptr();
    c1rc2_bptr->ce_bptr() = c1rc1_cbptr->ce_bptr();

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
                        ASSERT(tmp_pos.mca_cit == next(pos_next.mca_cit));
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
                    ASSERT(tmp_pos.mca_cit == next(pos_next.mca_cit));
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
            crt_mut_bptr->simplify(c1rc1_cbptr->ce_bptr()->substr(crt_mut_bptr->rf_start(), crt_mut_bptr->rf_len()));
            if (not crt_mut_bptr->is_empty())
            {
                // use an equivalent Mutation if one exists; else add this one
                crt_mut_bptr = mut_cont.find_equiv_or_add(crt_mut_bptr).unconst();
                // add MCA
                Mutation_Chunk_Adapter_BPtr mca_bptr = Mutation_Chunk_Adapter_Fact::new_elem(crt_mut_bptr, c1rc2_bptr);
                crt_mut_bptr->chunk_ptr_cont().insert(mca_bptr);
                c1rc2_bptr->mut_ptr_cont().push_back(mca_bptr);
                crt_mut_bptr = Mutation_Fact::new_elem();
            }
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
                ASSERT(c1rc1_cbptr->get_rc() or (pos.r_pos == next_r_mut_cbptr->rf_start()
                                                 and pos_next.r_pos == next_r_mut_cbptr->rf_end()));
                ASSERT(not c1rc1_cbptr->get_rc() or (pos.r_pos == next_r_mut_cbptr->rf_end()
                                                     and pos_next.r_pos == next_r_mut_cbptr->rf_start()));
                crt_mut_bptr->extend(pos.c_pos, pos_next.c_pos - pos.c_pos,
                                     //(not c1rc1_cbptr->get_rc()?
                                     //  next_r_mut_cbptr->get_seq()
                                     //  : next_r_mut_cbptr->get_seq().revcomp()));
                                     next_r_mut_cbptr->seq().revcomp(c1rc1_cbptr->get_rc()));
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
                       or (pos_next.mca_cit == next(pos.mca_cit) and pos_next.mut_offset == 0));
                //const Mutation& c_mut = *_mut_ptr_cont[pos.mut_idx];
                Mutation_CBPtr c_mut_cbptr = pos.mca_cit->mut_cbptr();
                crt_mut_bptr->extend(pos.c_pos, pos_next.c_pos - pos.c_pos,
                                     c_mut_cbptr->seq().substr(
                                                  min(pos.mut_offset, c_mut_cbptr->seq_len()),
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
        rc_cbptr->get_seq(), rc_cbptr->get_r_start());
    // basic settings
    new_rc_bptr->_r_start = rc_cbptr->get_c_start();
    new_rc_bptr->_r_len = rc_cbptr->get_c_len();
    new_rc_bptr->_c_start = rc_cbptr->get_r_start();
    new_rc_bptr->_c_len = rc_cbptr->get_r_len();
    new_rc_bptr->_set_rc(rc_cbptr->get_rc());
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
        ASSERT(pos_next.mca_cit == next(pos.mca_cit));
        // create reversed Mutation
        /*
        Mutation rev_mut((not _rc? pos.r_pos : pos_next.r_pos),
                         (not _rc? pos_next.r_pos - pos.r_pos : pos.r_pos - pos_next.r_pos),
                         (not _rc? _ce_ptr->substr(pos.c_pos, pos_next.c_pos - pos.c_pos)
                          : reverseComplement(_ce_ptr->substr(pos.c_pos, pos_next.c_pos - pos.c_pos))));
        */
        Seq_Proxy_Type mut_c_seq = rc_cbptr->ce_bptr()->substr(pos.c_pos, pos_next.c_pos - pos.c_pos);
        Mutation_BPtr new_mut_bptr = (
            not rc_cbptr->get_rc()?
            Mutation_Fact::new_elem(pos.r_pos, pos_next.r_pos - pos.r_pos, Seq_Type(mut_c_seq))
            : Mutation_Fact::new_elem(pos_next.r_pos, pos.r_pos - pos_next.r_pos, Seq_Type(mut_c_seq.revcomp())));
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
    _set_rc(not get_rc());
    _c_start = _ce_bptr->len() - (_c_start + _c_len);
}

void Read_Chunk::cat_c_right(Read_Chunk_BPtr rc_bptr, Read_Chunk_BPtr rc_next_bptr,
                             Mutation_Cont& mut_cont)
{
    LOG("Read_Chunk", debug1) << ptree("cat_c_right")
        .put("rc_ptr", rc_bptr.to_int())
        .put("rc_next_ptr", rc_next_bptr.to_int());

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
        and rc_bptr->mut_ptr_cont().rbegin()->mut_cbptr()->rf_end()
            == rc_next_bptr->mut_ptr_cont().begin()->mut_cbptr()->rf_start())
    {
        // create new merged Mutation
        Mutation_BPtr left_mut_bptr = rc_bptr->mut_ptr_cont().rbegin()->mut_cbptr().unconst();
        Mutation_BPtr right_mut_bptr = rc_next_bptr->mut_ptr_cont().begin()->mut_cbptr().unconst();
        Mutation_BPtr new_mut_bptr = Mutation_Fact::new_elem(left_mut_bptr->rf_start(),
                                                             left_mut_bptr->rf_len(),
                                                             Seq_Type(left_mut_bptr->seq()));
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
        // and deallocate them!
        Mutation_Chunk_Adapter_Fact::del_elem(left_mca_bptr);
        Mutation_Chunk_Adapter_Fact::del_elem(right_mca_bptr);
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
    rc_bptr->mut_ptr_cont().splice(rc_next_bptr->mut_ptr_cont(), rc_bptr);
    // deallocate rhs chunk
    Read_Chunk_Fact::del_elem(rc_next_bptr);
}

Range_Type Read_Chunk::mapped_range(Range_Type orig_rg, bool on_contig,
                                    bool rg_start_maximal, bool rg_end_maximal) const
{
    ASSERT(orig_rg.start() <= orig_rg.end());
    Range_Type rg(max(orig_rg.start(), on_contig? get_c_start() : get_r_start()),
                  min(orig_rg.end(), on_contig? get_c_end() : get_r_end()));
    ASSERT(rg.start() <= rg.end());
    // compute range endpoints in contig order
    Size_Type rg_endpoint_c_left;
    Size_Type rg_endpoint_c_right;
    bool rg_endpoint_c_left_maximal;
    bool rg_endpoint_c_right_maximal;
    if (on_contig or not get_rc())
    {
        rg_endpoint_c_left = rg.start();
        rg_endpoint_c_right = rg.end();
        rg_endpoint_c_left_maximal = rg_start_maximal;
        rg_endpoint_c_right_maximal = rg_end_maximal;
    }
    else
    {
        rg_endpoint_c_left = rg.end();
        rg_endpoint_c_right = rg.start();
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
    ASSERT(ret_rg_start <= ret_rg_end or (rg.start() == rg.end() and not rg_start_maximal and not rg_end_maximal));
    return Range_Type(ret_rg_start, ret_rg_end);
}

void Read_Chunk::make_unmappable(Read_Chunk_BPtr rc_bptr)
{
    Contig_Entry_BPtr ce_new_bptr = Contig_Entry_Fact::new_elem(rc_bptr->get_seq());
    ASSERT(ce_new_bptr->len() == rc_bptr->_r_len);
    rc_bptr->mut_ptr_cont().clear_and_dispose();
    rc_bptr->ce_bptr() = ce_new_bptr;
    rc_bptr->_c_start = 0;
    rc_bptr->_c_len = rc_bptr->_r_len;
    rc_bptr->_set_rc(false);
    rc_bptr->_set_is_unbreakable(true);
    ce_new_bptr->chunk_cont().insert(rc_bptr);
    ce_new_bptr->set_unmappable();
}

Seq_Type Read_Chunk::get_seq() const
{
    Seq_Type res;
    // start with start_pos iff not rc
    Pos pos = (not get_rc()? get_start_pos() : get_end_pos());
    while (pos != (not get_rc()? get_end_pos() : get_start_pos()))
    {
        Pos next_pos = pos;
        next_pos.advance(not get_rc());
        ASSERT(next_pos.r_pos >= pos.r_pos);
        ASSERT(next_pos.mut_offset == 0);
        Size_Type match_len = pos.get_match_len(not get_rc());
        if (match_len > 0)
        {
            // match stretch follows
            /*
            string tmp = _ce_bptr->substr((not get_rc()? pos.c_pos : next_pos.c_pos) - _ce_bptr->seq_offset(), match_len);
            res += (not get_rc()? tmp : reverseComplement(tmp));
            */
            res += _ce_bptr->substr((not get_rc()? pos.c_pos : next_pos.c_pos) - _ce_bptr->seq_offset(), match_len).revcomp(get_rc());
        }
        else if (pos.r_pos != next_pos.r_pos)
        {
            // mutation follows
            if (not get_rc())
            {
                ASSERT(not pos.past_last_mut());
                res += pos.mut().seq().substr(0, next_pos.r_pos - pos.r_pos);
            }
            else
            {
                res += next_pos.mut().seq().substr(0, next_pos.r_pos - pos.r_pos).revcomp();
            }
        }
        pos = next_pos;
    }
    return res;
}

void Read_Chunk::check() const
{
#ifndef BOOST_DISABLE_ASSERTS
    // check integrity of mutation pointer container
    mut_ptr_cont().check();
    // no empty chunks
    ASSERT(get_r_len() > 0);
    // contigs coordinates
    ASSERT(get_c_start() <= get_c_end());
    ASSERT(get_c_start() <= ce_bptr()->seq_offset() + ce_bptr()->len());
    ASSERT(get_c_end() <= ce_bptr()->seq_offset() + ce_bptr()->len());
    // mapped length
    Size_Type c_len = get_c_end() - get_c_start();
    Size_Type r_len = get_r_end() - get_r_start();
    long long delta = 0;
    Mutation_CBPtr last_mut_cbptr = nullptr;
    for (auto mca_cbptr : mut_ptr_cont() | referenced)
    {
        Mutation_CBPtr mut_cbptr = mca_cbptr->mut_cbptr();
        // no empty mutations
        ASSERT(not mut_cbptr->is_empty());
        // mutations must be in contig order
        ASSERT(not last_mut_cbptr or last_mut_cbptr->rf_end() <= mut_cbptr->rf_start());
#ifndef ALLOW_CONSECUTIVE_MUTATIONS
        ASSERT(not last_mut_cbptr or last_mut_cbptr->rf_end() < mut_cbptr->rf_start());
#endif
        delta += (long long)mut_cbptr->seq_len() - (long long)mut_cbptr->rf_len();
        last_mut_cbptr = mut_cbptr;
    }
    ASSERT((long long)c_len + delta == (long long)r_len);
#ifndef ALLOW_PART_MAPPED_INNER_CHUNKS
    // chunks must end on contig breaks except for first and last
    ASSERT(get_r_start() == 0
           or (not get_rc()?
               get_c_start() == 0
               : get_c_end() == ce_bptr()->seq_offset() + ce_bptr()->len()));
    ASSERT(get_r_end() == re_bptr()->len()
           or (not get_rc()?
               get_c_end() == ce_bptr()->seq_offset() + ce_bptr()->len()
               : get_c_start() == 0));
#endif
    // unmappable contigs have no mutations and a single chunk
    if (is_unbreakable())
    {
        ASSERT(ce_bptr()->is_unmappable() or ce_bptr()->is_lowcomplex());
        ASSERT(mut_ptr_cont().empty());
        ASSERT(ce_bptr()->chunk_cont().single_node());
        ASSERT((&*(ce_bptr()->chunk_cont().begin())).raw() == this);
        ASSERT(ce_bptr()->mut_cont().empty());
    }
#endif
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
                  .put("is_unbreakable", is_unbreakable())
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

string Read_Chunk::to_string(Read_Chunk_CBPtr rc_cbptr, bool r_dir, bool forward)
{
    bool r_forward;
    bool c_forward;
    if (r_dir)
    {
        r_forward = forward;
        c_forward = (r_forward != rc_cbptr->get_rc());
    }
    else
    {
        c_forward = forward;
        r_forward = (c_forward != rc_cbptr->get_rc());
    }

    Size_Type re_len = (rc_cbptr->re_bptr()?
                        rc_cbptr->re_bptr()->len() : rc_cbptr->get_r_len());
    Size_Type ce_len = rc_cbptr->ce_bptr()->len();
    ostringstream oss;
    oss << Read_Entry_Fact::size();
    size_t re_pad = oss.str().size();
    oss.str("");
    oss << Contig_Entry_Fact::size();
    size_t ce_pad = oss.str().size();
    oss.str("");
    oss << Read_Chunk_Fact::size();
    size_t rc_pad = oss.str().size();
    oss.str("");

    auto print_subinterval = [] (ostream& os,
                                 Size_Type subint_start, Size_Type subint_end,
                                 Size_Type int_start, Size_Type int_end, bool forward) {
        if (forward)
        {
            os << (int_start == subint_start? "[" : " ")
               << "["
               << setw(5) << right << subint_start
               << ","
               << setw(5) << right << subint_end
               << ")"
               << (int_end == subint_end? ")" : " ");
        }
        else // not forward
        {
            os << (int_end == subint_end? "(" : " ")
               << "("
               << setw(5) << right << subint_end
               << ","
               << setw(5) << right << subint_start
               << "]"
               << (int_start == subint_start? "]" : " ");
        }
    };

    oss << "rc " << setw(rc_pad) << left << rc_cbptr.to_int()
        << " re " << setw(re_pad) << left << rc_cbptr->re_bptr().to_int() << " ";
    print_subinterval(oss, rc_cbptr->get_r_start(), rc_cbptr->get_r_end(), 0, re_len, r_forward);
    oss << " " << setw(5) << left
        << (rc_cbptr->ce_bptr()->is_unmappable()? "unmap" : (not rc_cbptr->get_rc()? "fwd" : "rev"))
        << " ce " << setw(ce_pad) << left << rc_cbptr->ce_bptr().to_int() << " ";
    print_subinterval(oss, rc_cbptr->get_c_start(), rc_cbptr->get_c_end(), 0, ce_len, c_forward);
    return oss.str();
}


} // namespace MAC
