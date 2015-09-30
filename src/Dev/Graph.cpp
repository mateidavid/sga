#include "Graph.hpp"

#include <unordered_set>
#include <iomanip>
#include <cmath>
#include <boost/range/adaptor/filtered.hpp>

#include "Unmapper.hpp"
#include "overlapper.h"
//#include "Anchor_Support.hpp"
#include "zstr.hpp"


namespace MAC
{

Read_Chunk_CBPtr Graph::search_read_chunk(
    Contig_Entry_CBPtr ce_cbptr, Read_Entry_CBPtr re_cbptr,
    Size_Type pos, bool on_contig)
{
    if (on_contig)
    {
        auto iint_res = ce_cbptr->chunk_cont().iintersect(pos, pos);
        for (auto rc_cbptr : iint_res | referenced)
        {
            if (rc_cbptr->re_bptr() == re_cbptr)
            {
                return rc_cbptr;
            }
        }
    }
    else // not on_contig
    {
        for (auto rc_cbptr = re_cbptr->chunk_cont().get_chunk_with_pos(pos);
             rc_cbptr and rc_cbptr->get_r_start() <= pos;
             rc_cbptr = re_cbptr->chunk_cont().get_sibling(rc_cbptr, true, true))
        {
            if (rc_cbptr->ce_bptr() == ce_cbptr)
            {
                return rc_cbptr;
            }
        }
    }
    return nullptr;
}

Read_Chunk_CBPtr Graph::search_read_chunk_exact(
    Contig_Entry_CBPtr ce_cbptr, Read_Entry_CBPtr re_cbptr,
    Size_Type start_pos, Size_Type stop_pos, bool on_contig)
{
    if (on_contig)
    {
        auto iint_res = ce_cbptr->chunk_cont().iintersect(start_pos, stop_pos);
        for (auto rc_cbptr : iint_res | referenced)
        {
            if (rc_cbptr->re_bptr() == re_cbptr
                and rc_cbptr->get_c_start() == start_pos
                and rc_cbptr->get_c_end() == stop_pos)
            {
                return rc_cbptr;
            }
        }
    }
    else // not on_contig
    {
        for (auto rc_cbptr = re_cbptr->chunk_cont().get_chunk_with_pos(start_pos);
             rc_cbptr and rc_cbptr->get_r_start() <= start_pos;
             rc_cbptr = re_cbptr->chunk_cont().get_sibling(rc_cbptr, true, true))
        {
            if (rc_cbptr->ce_bptr() == ce_cbptr
                and rc_cbptr->get_r_start() == start_pos
                and rc_cbptr->get_r_end() == stop_pos)
            {
                return rc_cbptr;
            }
        }
    }
    return nullptr;
}

void Graph::add_read(string&& name, Seq_Type&& seq)
{
    LOG("graph", debug) << ptree("begin").put("name", name);
    // create RE & CE
    Read_Entry_BPtr re_bptr = Read_Entry_Fact::new_elem(move(name), seq.size());
    Contig_Entry_BPtr ce_bptr = Contig_Entry_Fact::new_elem(move(seq));
    // create initial chunk
    Read_Chunk_BPtr rc_bptr = Read_Chunk_Fact::new_elem(re_bptr, ce_bptr);
    // add chunk to RE and CE containers
    re_bptr->chunk_cont().insert(rc_bptr);
    ce_bptr->chunk_cont().insert(rc_bptr);
    // add RE and CE to graph containers
    re_cont().insert(re_bptr);
    ce_cont().insert(ce_bptr);
    // check new entry
    check(set< Read_Entry_CBPtr >({ re_bptr }));
}

void Graph::remove_read(Read_Entry_BPtr re_bptr)
{
    re_bptr->chunk_cont().clear_and_dispose(
        [&] (Read_Chunk_BPtr rc_bptr) {
            auto ce_bptr = rc_bptr->ce_bptr();
            ce_bptr->chunk_cont().erase(rc_bptr);
            Read_Chunk::dispose(rc_bptr);
            ce_bptr->mut_cont().drop_unused();
            if (ce_bptr->chunk_cont().empty())
            {
                ASSERT(ce_bptr->mut_cont().empty());
                ce_cont().erase(ce_bptr);
                Contig_Entry_Fact::del_elem(ce_bptr);
            }
        });
    re_cont().erase(re_bptr);
    Read_Entry_Fact::del_elem(re_bptr);
}

bool Graph::cut_contig_entry(Contig_Entry_BPtr ce_bptr, Size_Type c_brk, Mutation_CBPtr mut_left_cbptr)
{
    auto ce_new_bptr = ce_bptr->cut(c_brk, mut_left_cbptr);
    if (ce_new_bptr)
    {
        ce_cont().insert(ce_new_bptr);
    }
    return ce_new_bptr;
} // cut_contig_entry

bool Graph::cut_read_chunk(Read_Chunk_BPtr rc_bptr, Size_Type r_brk)
{
    LOG("graph", debug1) << ptree("begin")
        .put("rc_ptr", rc_bptr.to_ptree())
        .put("r_brk", r_brk);

    ASSERT(rc_bptr->get_r_len() > 0);
    ASSERT(rc_bptr->get_r_start() <= r_brk and r_brk <= rc_bptr->get_r_end());
    if (rc_bptr->is_unbreakable())
    {
        // never cut unbreakable chunks
        return false;
    }
    if (r_brk == rc_bptr->get_r_start() or r_brk == rc_bptr->get_r_end())
    {
        // cut is at the edge; it must be forced
        if ((r_brk == rc_bptr->get_r_start()) == (not rc_bptr->get_rc()))
        {
            // cut if mapping doesn't start on contig start
            if (rc_bptr->get_c_start() != 0)
            {
                return cut_contig_entry(rc_bptr->ce_bptr(), rc_bptr->get_c_start(), nullptr);
            }
            else
            {
                return false;
            }
        }
        else
        {
            // cut if mapping doesn't end on contig end
            if (rc_bptr->get_c_end() != rc_bptr->ce_bptr()->len())
            {
                // if chunk ends in insertion, keep it on lhs, along with the rest of the chunk
                Mutation_CBPtr mut_left_cbptr = nullptr;
                if (not rc_bptr->mut_ptr_cont().empty()
                    and rc_bptr->mut_ptr_cont().rbegin()->mut_cbptr()->rf_start() == rc_bptr->get_c_end())
                {
                    mut_left_cbptr = rc_bptr->mut_ptr_cont().rbegin()->mut_cbptr();
                    ASSERT(mut_left_cbptr->is_ins());
                }
                return cut_contig_entry(rc_bptr->ce_bptr(), rc_bptr->get_c_end(), mut_left_cbptr);
            }
            else
            {
                return false;
            }
        }
    }
    else
    {
        // cut is inside the read chunk; it must not be forced
        ASSERT(rc_bptr->get_r_start() < r_brk and r_brk < rc_bptr->get_r_end());
        bool cut_made = false;
        Read_Chunk::Pos pos = rc_bptr->get_start_pos();
        pos.jump_to_brk(r_brk, false);

        if (pos.mut_offset > 0)
        {
            // in the middle of a mutation; must cut it
            ASSERT(not pos.past_last_mut());
            Mutation_BPtr mut_bptr = pos.mca_cit->mut_cbptr().unconst();
            Contig_Entry_BPtr ce_bptr = rc_bptr->ce_bptr();
            ce_bptr->cut_mutation(mut_bptr,
                                  min(pos.mut_offset, mut_bptr->rf_len()),
                                  min(pos.mut_offset, mut_bptr->seq_len()));
            cut_made = true;
            // redo the jump
            pos = rc_bptr->get_start_pos();
            pos.jump_to_brk(r_brk, false);
        }
        // now we are certain the breakpoint no longer falls inside a mutation (insertion/mnp)
        ASSERT(pos.r_pos == r_brk);
        ASSERT(pos.mut_offset == 0);

        // check if an insertion at c_pos has to remain on the left of the cut
        Mutation_CBPtr mut_left_cbptr = nullptr;
        if (pos.past_first_mut()
            and pos.prev_mut().rf_start() == pos.c_pos)
        {
            auto tmp = pos.mca_cit;
            mut_left_cbptr = (--tmp)->mut_cbptr();
            ASSERT(mut_left_cbptr->is_ins());
        }

        // cut contig at given c_offset
        // all insertions at c_pos go to the right except for at most one, passed as argument
        return cut_contig_entry(rc_bptr->ce_bptr(), pos.c_pos, mut_left_cbptr) or cut_made;
    }
} // cut_read_chunk

bool Graph::cut_read_entry(Read_Entry_BPtr re_bptr, Size_Type r_brk)
{
    LOG("graph", debug1) << ptree("begin")
        .put("re_ptr", re_bptr.to_ptree())
        .put("r_brk", r_brk);

    if (r_brk <= re_bptr->start() or r_brk >= re_bptr->end())
    {
        // cut is on the edge of the read; it must be forced
        Read_Chunk_BPtr rc_bptr = (r_brk <= re_bptr->start()?
                                   &*re_bptr->chunk_cont().begin()
                                   : &*re_bptr->chunk_cont().rbegin());
        return cut_read_chunk(rc_bptr, r_brk <= re_bptr->start()? re_bptr->start() : re_bptr->end());
    }
    else
    {
        // cut is inside the read, it must not be forced
        Read_Chunk_BPtr rc_bptr = re_bptr->chunk_cont().get_chunk_with_pos(r_brk).unconst();
        // the chunk to be cut must exist
        ASSERT(rc_bptr->get_r_start() <= r_brk and r_brk < rc_bptr->get_r_end());
        return cut_read_chunk(rc_bptr, r_brk);
    }
}

void Graph::merge_chunk_contigs(Read_Chunk_BPtr c1rc1_chunk_bptr, Read_Chunk_BPtr c2rc2_chunk_bptr, Cigar& rc1rc2_cigar)
{
    rc1rc2_cigar.check(c1rc1_chunk_bptr->get_seq(), c2rc2_chunk_bptr->get_seq());
    ASSERT(c1rc1_chunk_bptr->get_c_start() == 0 and c1rc1_chunk_bptr->get_c_end() == c1rc1_chunk_bptr->ce_bptr()->len());
    ASSERT(c2rc2_chunk_bptr->get_c_start() == 0 and c2rc2_chunk_bptr->get_c_end() == c2rc2_chunk_bptr->ce_bptr()->len());
    // do not do anything if the chunks are already mapped to the same contig
    // NOTE: with this, we are ignoring alternate mappings
    Contig_Entry_BPtr c1_ce_bptr = c1rc1_chunk_bptr->ce_bptr();
    Contig_Entry_BPtr c2_ce_bptr = c2rc2_chunk_bptr->ce_bptr();
    if (c1_ce_bptr == c2_ce_bptr)
    {
        return;
    }

    LOG("graph", debug1) << ptree("begin")
        .put("c1rc1_chunk_ptr", c1rc1_chunk_bptr.to_int())
        .put("c1", (*c1_ce_bptr).to_ptree())
        .put("c2rc2_chunk_ptr", c2rc2_chunk_bptr.to_int())
        .put("c2", (*c2_ce_bptr).to_ptree())
        .put("rc1rc2_cigar", rc1rc2_cigar.to_ptree());

    // construct read chunk for the rc1->rc2 mapping
    Read_Chunk_BPtr rc1rc2_chunk_bptr;
    rc1rc2_chunk_bptr = Read_Chunk::make_relative_chunk(c1rc1_chunk_bptr, c2rc2_chunk_bptr,
                                                        rc1rc2_cigar);

    // construct c1<-rc2 mapping directly from c1<-rc1 & rc1<-rc2
    Read_Chunk_BPtr c1rc2_chunk_bptr = Read_Chunk::collapse_mapping(c1rc1_chunk_bptr, rc1rc2_chunk_bptr,
                                                                    c1_ce_bptr->mut_cont());

    // next, construct rc2<-c2 mapping by reversing c2<-rc2
    Read_Chunk_BPtr rc2c2_chunk_bptr = Read_Chunk::invert_mapping(c2rc2_chunk_bptr);

    // construct c1<-c2 mapping from c1<-rc2 and rc2<-c2
    Read_Chunk_BPtr c1c2_chunk_bptr = Read_Chunk::collapse_mapping(c1rc2_chunk_bptr, rc2c2_chunk_bptr,
                                                                   c1_ce_bptr->mut_cont());

    // for the remaining chunks mapped to c2, remap them to c1 through c1->c2 mapping
    for (auto c2rc_bptr : c2_ce_bptr->chunk_cont() | referenced)
    {
        Read_Chunk_BPtr c1rc_bptr;
        if (c2rc_bptr == c2rc2_chunk_bptr)
        {
            c1rc_bptr = c1rc2_chunk_bptr;
        }
        else
        {
            c1rc_bptr = Read_Chunk::collapse_mapping(c1c2_chunk_bptr, c2rc_bptr, c1_ce_bptr->mut_cont());
        }
        // unlink old chunk from RE cont, and link new one
        Read_Entry_BPtr re_bptr = c2rc_bptr->re_bptr();
        re_bptr->chunk_cont().erase(c2rc_bptr);
        re_bptr->chunk_cont().insert(c1rc_bptr);
        c1_ce_bptr->chunk_cont().insert(c1rc_bptr);
    }
    // all read chunks mapped to c2 are translated to c1

    // next, deallocate temporary structures
    // rc1rc2_chunk_bptr:
    Contig_Entry_BPtr rc1_ce_bptr = rc1rc2_chunk_bptr->ce_bptr();
    ASSERT(size_one(rc1_ce_bptr->chunk_cont()));
    rc1_ce_bptr->chunk_cont().clear_and_dispose();
    rc1_ce_bptr->mut_cont().clear_and_dispose();
    Contig_Entry_Fact::del_elem(rc1_ce_bptr);

    // rc2c2_chunk_bptr
    Contig_Entry_BPtr rc2_ce_bptr = rc2c2_chunk_bptr->ce_bptr();
    ASSERT(size_one(rc2_ce_bptr->chunk_cont()));
    rc2_ce_bptr->chunk_cont().clear_and_dispose();
    rc2_ce_bptr->mut_cont().clear_and_dispose();
    Contig_Entry_Fact::del_elem(rc2_ce_bptr);

    // c1c2_chunk_bptr
    c1c2_chunk_bptr->mut_ptr_cont().clear_and_dispose();
    Read_Chunk_Fact::del_elem(c1c2_chunk_bptr);
    c1_ce_bptr->mut_cont().drop_unused();

    // deallocate everything in c2
    ce_cont().erase(c2_ce_bptr);
    c2_ce_bptr->chunk_cont().clear_and_dispose();
    c2_ce_bptr->mut_cont().clear_and_dispose();
    Contig_Entry_Fact::del_elem(c2_ce_bptr);

    //cerr << "after merging read chunks:\n" << *this;
    check(set< Contig_Entry_CBPtr >( { c1_ce_bptr }));
}

vector< tuple< Size_Type, Size_Type, Cigar > >
Graph::chunker(Read_Entry_BPtr re1_bptr, Read_Entry_BPtr re2_bptr, Cigar& cigar)
{
    vector< tuple< Size_Type, Size_Type, Cigar > > rc_mapping;

    auto r1_start = re1_bptr->start();
    auto r1_end = re1_bptr->end();
    auto r2_start = re2_bptr->start();
    auto r2_end = re2_bptr->end();

    // we repeatedly cut the read entries of either read
    // until their chunks match in the way described by the cigar string
    // keep track of read chunk mapping, and cigar transformation between them
    bool done;
    while (true)
    {
        if (trim_tuc_step()
            and (r1_start != re1_bptr->start()
                 or r1_end != re1_bptr->end()
                 or r2_start != re2_bptr->start()
                 or r2_end != re2_bptr->end()))
        {
            // apply trim if any occurred during unmapping operations
            r1_start = re1_bptr->start();
            r1_end = re1_bptr->end();
            r2_start = re2_bptr->start();
            r2_end = re2_bptr->end();
            cigar.trim(r1_start, r1_end, r2_start, r2_end);
            cigar.drop_terminal_indels();
            cut_read_entry(re1_bptr, cigar.rf_start());
            cut_read_entry(re1_bptr, cigar.rf_end());
            cut_read_entry(re2_bptr, cigar.qr_start());
            cut_read_entry(re2_bptr, cigar.qr_end());
        }

        Size_Type r1_match_start = cigar.rf_start();
        Size_Type r1_match_len = cigar.rf_len();
        static_cast< void >(r1_match_len);
        Size_Type r2_match_start = cigar.qr_start();
        Size_Type r2_match_len = cigar.qr_len();
        bool r2_rc = cigar.reversed();

        // after every graph modification, restart at the beginning
        done = true;
        rc_mapping.clear();

        /// re1 part matched so far: [r1_match_start, r1_pos)
        Size_Type r1_pos = r1_match_start;
        /// re2 part matched so far: [r2_match_start, r2_pos) or [r2_pos, r2_end) if rc
        Size_Type r2_pos = (not r2_rc? r2_match_start : r2_match_start + r2_match_len);
        /// cigar ops matched so far: [0, op_start)
        size_t op_start = 0;
        while (op_start < cigar.n_ops())
        {
            /// next chunk in re1
            Read_Chunk_BPtr rc1_bptr = re1_bptr->chunk_cont().get_chunk_with_pos(r1_pos).unconst();
            ASSERT(rc1_bptr);
            /// re1 position relative to rc1_r_start
            Size_Type rc1_offset = r1_pos - rc1_bptr->get_r_start();
            /// re1 length left in rc1 after partial match
            Size_Type rc1_remaining_len = rc1_bptr->get_r_len() - rc1_offset;
            /// next chunk in re2
            Read_Chunk_BPtr rc2_bptr = re2_bptr->chunk_cont().get_chunk_with_pos(not r2_rc? r2_pos : r2_pos - 1).unconst();
            ASSERT(rc2_bptr);
            /// re2 position relative to rc2_r_start or rc2_r_end if rc
            Size_Type rc2_offset = (not r2_rc
                                    ? r2_pos - rc2_bptr->get_r_start()
                                    : rc2_bptr->get_r_end() - r2_pos);
            /// re2 length left in rc2 after partial match
            Size_Type rc2_remaining_len = rc2_bptr->get_r_len() - rc2_offset;

            // invariant: we matched read 1 chunks before rc1
            // to read 2 chunks before/after rc2
            // using cigar ops before op_start
            ASSERT(r1_pos < r1_match_start + r1_match_len);
            ASSERT(r2_rc or r2_pos < r2_match_start + r2_match_len);
            ASSERT(not r2_rc or r2_pos > r2_match_start);
            ASSERT(rc1_offset == 0 or rc1_bptr->is_unbreakable());
            ASSERT(rc2_offset == 0 or rc2_bptr->is_unbreakable());
            ASSERT(rc1_bptr);
            ASSERT(rc1_bptr->get_r_start() + rc1_offset == r1_pos);
            ASSERT(rc1_remaining_len > 0);
            ASSERT(rc2_bptr);
            ASSERT(r2_rc or rc2_bptr->get_r_start() + rc2_offset == r2_pos);
            ASSERT(not r2_rc or rc2_bptr->get_r_end() - rc2_offset == r2_pos);
            ASSERT(rc2_remaining_len > 0);
            ASSERT(r1_pos == cigar.op_rf_pos(op_start));
            ASSERT(r2_pos == cigar.op_qr_pos(op_start));

            // advance past cigar ops until either chunk ends
            size_t op_cnt = 1;
            size_t op_end = op_start + op_cnt;
            while (op_end < cigar.n_ops()
                    and (cigar.op_rf_len(op_start, op_cnt) < rc1_remaining_len
                         or (cigar.op_rf_len(op_start, op_cnt) == rc1_remaining_len
                             and cigar.op_is_insertion(op_end)))
                    and (cigar.op_qr_len(op_start, op_cnt) < rc2_remaining_len
                         or (cigar.op_qr_len(op_start, op_cnt) == rc2_remaining_len
                             and cigar.op_is_deletion(op_end))))
            {
                ++op_cnt;
                ++op_end;
            }
            Size_Type op_rf_sub_len = cigar.op_rf_len(op_start, op_cnt);
            Size_Type op_qr_sub_len = cigar.op_qr_len(op_start, op_cnt);
            // stop conditions: can be derived by inspecting loop
            ASSERT(// past all ops
                   op_end == cigar.n_ops()
                   // or current chunks perfectly matched by cigar
                   or (op_rf_sub_len == rc1_remaining_len
                       and op_qr_sub_len == rc2_remaining_len)
                   // or cigar op extends past rc1
                   or op_rf_sub_len > rc1_remaining_len
                   // or cigar op extends past rc2
                   or op_qr_sub_len > rc2_remaining_len
                   // or rc1 mapped fully, but rc2 has non-insertion op at the end
                   or (op_rf_sub_len == rc1_remaining_len
                       and op_qr_sub_len < rc2_remaining_len
                       and not cigar.op_is_insertion(op_end))
                   // or rc2 mapped fully, but rc1 has non-deletion op at the end
                   or (op_rf_sub_len < rc1_remaining_len
                       and op_qr_sub_len == rc2_remaining_len
                       and not cigar.op_is_deletion(op_end)));
            // we previously cut the chunks at the ends of the mapping
            // cuts don't succeed only if those chunks are unmappable
            ASSERT(op_end < cigar.n_ops()
                   or op_rf_sub_len >= rc1_remaining_len
                   or rc1_bptr->is_unbreakable());
            ASSERT(op_end < cigar.n_ops()
                   or op_qr_sub_len >= rc2_remaining_len
                   or rc2_bptr->is_unbreakable());
            if (op_rf_sub_len > rc1_remaining_len)
            {
                // the first inequality trivially holds with <=
                // but it can be shown it holds in fact with <
                ASSERT(cigar.op_rf_pos(op_end - 1) < rc1_bptr->get_r_end()
                       and rc1_bptr->get_r_end() < cigar.op_rf_pos(op_end));
            }
            if (op_qr_sub_len > rc2_remaining_len)
            {
                // as with rf, the inequalities involving op_end-1 tivially hold with <=
                // but it can be shown they hold with <
                ASSERT(r2_rc or (cigar.op_qr_pos(op_end - 1) < rc2_bptr->get_r_end()
                       and rc2_bptr->get_r_end() < cigar.op_qr_pos(op_end)));
                ASSERT(not r2_rc or (cigar.op_qr_pos(op_end) < rc2_bptr->get_r_start()
                       and rc2_bptr->get_r_start() < cigar.op_qr_pos(op_end - 1)));
            }

            // check if either chunk ended during the last cigar op
            if (op_rf_sub_len > rc1_remaining_len
                or op_qr_sub_len > rc2_remaining_len)
            {
                // find out which of the 2 ended earlier, cut the cigar op at that position
                Size_Type r1_break_len = 0;
                if (op_rf_sub_len > rc1_remaining_len
                    and not cigar.op_is_insertion(op_end - 1))
                {
                    r1_break_len = cigar.get_rf_op_prefix_len(op_end - 1, rc1_bptr->get_r_end());
                    ASSERT(0 < r1_break_len and r1_break_len < cigar.op_rf_len(op_end - 1));
                }

                Size_Type r2_break_len = 0;
                if (op_qr_sub_len > rc2_remaining_len
                    and not cigar.op_is_deletion(op_end -1))
                {
                    r2_break_len = cigar.get_qr_op_prefix_len(op_end - 1, (not r2_rc?
                                                                           rc2_bptr->get_r_end()
                                                                           : rc2_bptr->get_r_start()));
                    ASSERT(0 < r2_break_len and r2_break_len < cigar.op_qr_len(op_end - 1));
                }

                ASSERT(r1_break_len > 0 or r2_break_len > 0);
                cigar.cut_op(op_end - 1, (r1_break_len == 0?
                                          r2_break_len
                                          : (r2_break_len == 0?
                                             r1_break_len
                                             : min(r1_break_len, r2_break_len))));
                // recompute sub-lengths after the cut
                op_rf_sub_len = cigar.op_rf_len(op_start, op_cnt);
                op_qr_sub_len = cigar.op_qr_len(op_start, op_cnt);
            }

            ASSERT(op_rf_sub_len <= rc1_remaining_len
                   and op_qr_sub_len <= rc2_remaining_len);
            // now we are sure the op ends on at least one of the read chunk boundaries
            // unless these are the last chunks and they are both unbreakable
            if (op_end == cigar.n_ops() and rc1_bptr->is_unbreakable() and rc2_bptr->is_unbreakable())
            {
                ASSERT(done);
                break;
            }
            ASSERT(op_rf_sub_len == rc1_remaining_len
                   or op_qr_sub_len == rc2_remaining_len);
            // if it doesn't end on both, we might need to cut the other
            if (op_qr_sub_len < rc2_remaining_len)
            {
                // it follows from main stop condition that the next op is not an insertion
                ASSERT(op_end == cigar.n_ops() or not cigar.op_is_insertion(op_end));
                if (op_qr_sub_len == 0)
                {
                    // no progress on rc2: rc1 is mapped entirely to a deletion
                    // move on to next chunk on read 1
                    ASSERT(op_rf_sub_len == rc1_remaining_len
                           and op_rf_sub_len > 0);
                    r1_pos = rc1_bptr->get_r_end();
                    op_start = op_end;
                    continue;
                }
                else if (rc2_bptr->is_unbreakable())
                {
                    // rc2 cannot be cut
                    if (not rc1_bptr->is_unbreakable())
                    {
                        // rc1 is mapped to an unbreakable chunk; we unmap it
                        unmap_chunk(rc1_bptr);
                        done = false;
                        break;
                    }
                    else // both rc1 and rc2 are unbreakable
                    {
                        // we advance both
                        r1_pos = rc1_bptr->get_r_end();
                        r2_pos = (not r2_rc?
                                  r2_pos + op_qr_sub_len
                                  : r2_pos - op_qr_sub_len);
                        op_start = op_end;
                        continue;
                    }
                }
                else
                {
                    cut_read_entry(re2_bptr, cigar.op_qr_pos(op_end));
                    done = false;
                    break;
                }
            }
            if (op_rf_sub_len < rc1_remaining_len)
            {
                // it follows from main stop condition that the next op is not a deletion
                ASSERT(op_end == cigar.n_ops() or not cigar.op_is_deletion(op_end));
                if (op_rf_sub_len == 0)
                {
                    // no progress on rc1: rc2 mapped entirely to an insertion
                    // move on to next chunk on read 2
                    ASSERT(op_qr_sub_len == rc2_remaining_len
                           and op_qr_sub_len > 0);
                    r2_pos = (not r2_rc? rc2_bptr->get_r_end() : rc2_bptr->get_r_start());
                    op_start = op_end;
                    continue;
                }
                else if (rc1_bptr->is_unbreakable())
                {
                    // rc1 cannot be cut
                    if (not rc2_bptr->is_unbreakable())
                    {
                        // rc2 is mapped to an unbreakable chunk; we unmap it
                        unmap_chunk(rc2_bptr);
                        done = false;
                        break;
                    }
                    else // both rc1 and rc2 are unbreakable
                    {
                        // we advance both
                        r1_pos += op_rf_sub_len;
                        r2_pos = (not r2_rc? rc2_bptr->get_r_end() : rc2_bptr->get_r_start());
                        op_start = op_end;
                        continue;
                    }
                }
                else
                {
                    cut_read_entry(re1_bptr, cigar.op_rf_pos(op_end));
                    done = false;
                    break;
                }
            }
            // reached when both rc1 and rc2 end at the current cigar op
            ASSERT(op_rf_sub_len == rc1_remaining_len
                   and op_qr_sub_len == rc2_remaining_len);
            if (rc1_bptr->is_unbreakable() and not rc2_bptr->is_unbreakable())
            {
                unmap_chunk(rc2_bptr);
                done = false;
                break;
            }
            if (rc2_bptr->is_unbreakable() and not rc1_bptr->is_unbreakable())
            {
                unmap_chunk(rc1_bptr);
                done = false;
                break;
            }
            ASSERT(rc1_bptr->is_unbreakable() == rc2_bptr->is_unbreakable());
            if (not rc1_bptr->is_unbreakable())
            {
                // add read chunk mapping
                rc_mapping.push_back(make_tuple(rc1_bptr->get_r_start(), rc2_bptr->get_r_start(), cigar.subcigar(op_start, op_cnt)));
            }
            // advance both chunks and cigar
            r1_pos = rc1_bptr->get_r_end();
            r2_pos = (not r2_rc? rc2_bptr->get_r_end() : rc2_bptr->get_r_start());
            op_start = op_end;
        } // while (op_start < cigar.n_ops())
        if (done)
        {
            break;
        }
    } // while (true)
    return rc_mapping;
} // chunker

void Graph::add_overlap(const string& r1_name, const string& r2_name,
                        Size_Type r1_start, Size_Type r1_len,
                        Size_Type r2_start, Size_Type r2_len, bool r2_rc,
                        const string& cigar_string)
{
    LOG("graph", debug) << ptree("begin")
        .put("r1_name", r1_name)
        .put("r2_name", r2_name)
        .put("r1_start", r1_start)
        .put("r1_len", r1_len)
        .put("r2_start", r2_start)
        .put("r2_len", r2_len)
        .put("r2_rc", r2_rc)
        .put("cigar", cigar_string);

    // fetch read entries
    Read_Entry_BPtr re1_bptr = re_cont().find(r1_name).unconst();
    if (not re1_bptr)
    {
        LOG("graph", warning) << ptree("missing_read").put("r_name", r1_name);
        return;
    }
    ASSERT(re1_bptr);
    Read_Entry_BPtr re2_bptr = re_cont().find(r2_name).unconst();
    if (not re2_bptr)
    {
        LOG("graph", warning) << ptree("missing_read").put("r_name", r2_name);
        return;
    }
    ASSERT(re2_bptr);
    // construct cigar object
    Cigar cigar(cigar_string, r2_rc, r1_start, r2_start);
    ASSERT(r1_len == cigar.rf_len());
    ASSERT(r2_len == cigar.qr_len());
    // fix cigar if read entries are trimmed
    cigar.trim(re1_bptr->start(), re1_bptr->end(), re2_bptr->start(), re2_bptr->end());
    if (cigar.rf_len() == 0 or cigar.qr_len()== 0)
    {
        return;
    }
    r1_start = cigar.rf_start();
    r1_len = cigar.rf_len();
    r2_start = cigar.qr_start();
    r2_len = cigar.qr_len();
    auto r1_start_offset = r1_start - re1_bptr->start();
    auto r2_start_offset = r2_start - re2_bptr->start();
    // disambiguate cigar
    Seq_Type r1_seq = re1_bptr->get_seq();
    Seq_Type r2_seq = re2_bptr->get_seq();
    LOG("graph", debug1) << ptree("before_disambiguate")
        .put("re1", r1_seq.substr(r1_start_offset, r1_len))
        .put("re2", r2_seq.substr(r2_start_offset, r2_len).revcomp(cigar.reversed()))
        .put("cigar", cigar.to_ptree());
    cigar.disambiguate(r1_seq.substr(r1_start_offset, r1_len),
                       r2_seq.substr(r2_start_offset, r2_len));
    LOG("graph", debug1) << ptree("after_disambiguate").put("cigar", cigar.to_ptree());
    // check cigar
    cigar.check(r1_seq.substr(r1_start_offset, r1_len),
                r2_seq.substr(r2_start_offset, r2_len));
    // discard indels at either end of the cigar string
    cigar.drop_terminal_indels();
    if (cigar.rf_len() == 0)
    {
        return;
    }
    // cut r1 & r2 at the ends of the match region
    // NOTE: unbreakable chunks are never cut
    cut_read_entry(re1_bptr, cigar.rf_start());
    cut_read_entry(re1_bptr, cigar.rf_end());
    cut_read_entry(re2_bptr, cigar.qr_start());
    cut_read_entry(re2_bptr, cigar.qr_end());
    // break reads into corresponding chunks
    auto rc_mapping = chunker(re1_bptr, re2_bptr, cigar);
    // reached when we have a complete rc map
    for (auto& tmp : rc_mapping)
    {
        Size_Type rc1_start;
        Size_Type rc2_start;
        Cigar rc1rc2_cigar;
        tie(rc1_start, rc2_start, rc1rc2_cigar) = move(tmp);
        Read_Chunk_BPtr rc1_bptr = re1_bptr->chunk_cont().get_chunk_with_pos(rc1_start).unconst();
        ASSERT(rc1_bptr->get_r_start() == rc1_start);
        Read_Chunk_BPtr rc2_bptr = re2_bptr->chunk_cont().get_chunk_with_pos(rc2_start).unconst();
        ASSERT(rc2_bptr->get_r_start() == rc2_start);
        merge_chunk_contigs(rc1_bptr, rc2_bptr, rc1rc2_cigar);
    }
    // find unmappable regions in the contigs recently merged
    auto region_cont = find_unmappable_regions(re1_bptr, cigar.rf_start(), cigar.rf_end());
    /*
    for (const auto& rg : region_cont)
    {
        unmap_re_region(re1_bptr, rg);
    }
    */
    if (not region_cont.empty())
    {
        unmap_re_regions(re1_bptr, move(region_cont));
    }
    // check read entries
    check(set< Read_Entry_CBPtr >({ re1_bptr, re2_bptr }));
    // catenate at the end of each step
    if (cat_at_step())
    {
        cat_read_contigs(re1_bptr);
        check(set< Read_Entry_CBPtr >({ re1_bptr, re2_bptr }));
    }
} // Graph::add_overlap

bool Graph::cat_contigs(Contig_Entry_BPtr ce_bptr, bool c_right)
{
    LOG("graph", debug1) << ptree("begin")
        .put("ce_ptr", ce_bptr.to_int())
        .put("c_right", c_right)
        .put("ce_ptr->category", ce_bptr->category());

    /** !!! GCC BUG !!!
     * The code below demonstrates a bug in gcc 4.9.1 & 4.8.1 when compiled with -O2.
     * Specifically, the code returns false in the statement immediately below,
     * even when ce_next_bptr != nullptr.
     *
     * By commenting out the 2 identical print statements below, one sees the printed value
     * change between the statements (!?) The first value is sometimes incorrect: null when
     * it should be non-null. The presence of print lines causes the if statement to be
     * correctly evaluated, hiding the bug.
     *
     * Commenting out the completely unrelated "hello" print statement also causes proper
     * evaluation of the following if statement.
     */
#ifdef ENABLE_GCC_BUG
    Contig_Entry_BPtr ce_next_cbptr;
    bool same_orientation;
    set< Read_Chunk_CBPtr > rc_cbptr_cont;
    tie(ce_next_cbptr, same_orientation, rc_cbptr_cont) = ce_bptr->can_cat_dir(c_right);
    Contig_Entry_BPtr ce_next_bptr(ce_next_cbptr.unconst());

#ifdef ENABLE_GCC_BUG_SHOW_VALUE_CHANGE
    cerr << "ce_next_bptr (1):" << ce_next_bptr.to_int() << endl;
    cerr << "ce_next_bptr (2):" << ce_next_bptr.to_int() << endl;
#elif ENABLE_GCC_BUG_FIX_WITH_HELLO
    cerr << "hello\n";
#endif
#else
    auto tmp = ce_bptr->can_cat_dir(c_right);
    Contig_Entry_BPtr ce_next_bptr = get<0>(tmp).unconst();
    bool same_orientation = get<1>(tmp);
    set< Read_Chunk_CBPtr > rc_cbptr_cont = move(get<2>(tmp));
#endif

    if (not ce_next_bptr)
    {
        return false;
    }
    LOG("graph", debug1) << ptree("mid")
        .put("ce_next_ptr", ce_next_bptr.to_int())
        .put("same_orientation", same_orientation)
        .put("ce_next_ptr->category", ce_next_bptr->category());
    if (ce_bptr->category() != ce_next_bptr->category())
    {
        // don't merge contigs with different categories
        return false;
    }

    if (not same_orientation)
    {
        // we reverse one of the contigs
        // reversal modifications are done in-place, so chunk vectors are still valid
        ce_bptr->reverse();
        ce_bptr->check();
        c_right = not c_right;
        same_orientation = true;
    }
    if (not c_right)
    {
        swap(ce_bptr, ce_next_bptr);
        c_right = true;
        auto tmp2 = ce_bptr->can_cat_dir(true);
        rc_cbptr_cont = move(get<2>(tmp2));
        ASSERT(not rc_cbptr_cont.empty());
    }
    // at this point:
    // - contigs are in the same orientation
    // - ce_bptr is the start of the merge and ce_next_bptr is the tail
    // - rc_cbptr_cont contains chunks from ce_bptr that go into ce_next_bptr

    LOG("graph", debug1) << ptree("merge")
        .put("ce_ptr", ce_bptr.to_int())
        .put("ce_next_ptr", ce_next_bptr.to_int());

    ce_cont().erase(ce_next_bptr);
    Contig_Entry::cat_c_right(ce_bptr, ce_next_bptr, rc_cbptr_cont);
    ce_bptr->check();

    return true;
} // cat_contigs

void Graph::cat_read_contigs(Read_Entry_BPtr re_bptr)
{
    LOG("graph", debug1) << ptree("begin").put("re", re_bptr.to_int());

    bool done = false;
    while (not done)
    {
        done = true;
        for (auto rc_bptr : re_bptr->chunk_cont() | referenced)
        {
            Contig_Entry_BPtr ce_bptr = rc_bptr->ce_bptr();
            if (ce_bptr->is_unmappable())
            {
                continue;
            }
            // first try to merge first contig to the left of read start
            if (rc_bptr == &*re_bptr->chunk_cont().begin())
            {
                if (cat_contigs(ce_bptr, rc_bptr->get_rc()))
                {
                    re_bptr->check();
                    done = false;
                    break;
                }
            }
            // then every contig to the right of its chunk
            if (cat_contigs(ce_bptr, not rc_bptr->get_rc()))
            {
                re_bptr->check();
                done = false;
                break;
            }
        }
    }
} // cat_read_contigs

void Graph::cat_all_read_contigs()
{
    LOG("graph", info) << ptree("begin");
    for (auto re_bptr : re_cont() | referenced)
    {
        cat_read_contigs(re_bptr);
    }
    LOG("graph", info) << ptree("end");
}

/*
Read_Chunk_BPtr Graph::trim_contig_to_chunk(Read_Chunk_BPtr rc_bptr)
{
    // we save read entry and offset, and recompute chunk after graph modifications
    Read_Entry_BPtr re_bptr = rc_bptr->re_bptr();
    Size_Type rc_start = rc_bptr->get_r_start();
    Size_Type rc_end = rc_bptr->get_r_end();
    Contig_Entry_BPtr ce_bptr = rc_bptr->ce_bptr();

    LOG("graph", debug2) << ptree("trim_contig_to_chunk")
        .put("re_name", re_bptr->name())
        .put("re_ptr", re_bptr.to_int())
        .put("start", rc_start)
        .put("end", rc_end);

    // trim contig entry at the start of the read chunk
    if (rc_bptr->get_c_start() != 0)
    {
        auto ce_mid_bptr = ce_bptr->cut(rc_bptr->get_c_start(), nullptr, false);
        if (ce_mid_bptr)
        {
            ce_cont().insert(ce_mid_bptr);
            ce_bptr = ce_mid_bptr;
        }
        rc_bptr = search_read_chunk_exact(ce_bptr, re_bptr, rc_start, rc_end, false).unconst();
        // chunk should survive the cut
        ASSERT(rc_bptr);
        ASSERT(rc_bptr->get_r_start() == rc_start);
        ASSERT(rc_bptr->get_r_end() == rc_end);
    }

    // trim contig entry at the end of the chunk
    if (rc_bptr->get_c_end() != ce_bptr->len())
    {
        Mutation_CBPtr last_ins_cbptr = nullptr;
        if (not rc_bptr->mut_ptr_cont().empty()
            and rc_bptr->mut_ptr_cont().rbegin()->mut_cbptr()->is_ins()
            and rc_bptr->mut_ptr_cont().rbegin()->mut_cbptr()->rf_start() == rc_bptr->get_c_end())
        {
            last_ins_cbptr = rc_bptr->mut_ptr_cont().rbegin()->mut_cbptr();
        }
        auto ce_rhs_bptr = ce_bptr->cut(rc_bptr->get_c_end(), last_ins_cbptr, false);
        if (ce_rhs_bptr)
        {
            ce_cont().insert(ce_rhs_bptr);
        }
        rc_bptr = search_read_chunk_exact(ce_bptr, re_bptr, rc_start, rc_end, false).unconst();
        // chunk should survive the cut
        ASSERT(rc_bptr);
        ASSERT(rc_bptr->get_r_start() == rc_start);
        ASSERT(rc_bptr->get_r_end() == rc_end);
    }
    // at this point, the chunk should span the entire contig
    ASSERT(rc_bptr->get_c_start() == 0 and rc_bptr->get_c_end() == ce_bptr->len());
    return rc_bptr;
}
*/

void Graph::trim_tuc(Read_Chunk_BPtr rc_bptr)
{
    ASSERT(rc_bptr);
    Read_Entry_BPtr re_bptr = rc_bptr->re_bptr();
    Contig_Entry_BPtr ce_bptr = rc_bptr->ce_bptr();
    ASSERT(ce_bptr->is_unmappable());
    ASSERT(rc_bptr->get_r_start() == re_bptr->start() or rc_bptr->get_r_end() == re_bptr->end());
    // adjust re
    if (rc_bptr->get_r_start() == re_bptr->start())
    {
        re_bptr->start() = rc_bptr->get_r_end();
    }
    re_bptr->len() -= rc_bptr->get_r_len();
    // erase rc
    re_bptr->chunk_cont().erase(rc_bptr);
    ce_bptr->chunk_cont().erase(rc_bptr);
    Read_Chunk_Fact::del_elem(rc_bptr);
    // erase ce
    ce_cont().erase(ce_bptr);
    Contig_Entry_Fact::del_elem(ce_bptr);
}

void Graph::trim_tucs()
{
    LOG("Graph", info) << ptree("begin");
    for (auto re_bptr : re_cont() | referenced)
    {
        if (re_bptr->len() == 0) continue;
        ASSERT(not re_bptr->chunk_cont().empty());
        for (int dir = 0; dir < 2; ++dir)
        {
            auto rc_bptr = (dir == 0
                            ? &*re_bptr->chunk_cont().begin()
                            : &*re_bptr->chunk_cont().rbegin());
            if (rc_bptr->ce_bptr()->is_unmappable())
            {
                trim_tuc(rc_bptr);
            }
        }
    }
    LOG("Graph", info) << ptree("end");
}

void Graph::unmap_chunk(Read_Chunk_BPtr rc_bptr)
{ Unmapper{*this}.unmap_chunk(rc_bptr); }
void Graph::unmap_re_regions(Unmapper::re_set_type&& unmap_re_set)
{ Unmapper{*this}.unmap_re_regions(move(unmap_re_set)); }
void Graph::unmap_re_regions(Read_Entry_BPtr re_bptr, Range_Cont&& rg_cont)
{ Unmapper{*this}.unmap_re_regions(re_bptr, move(rg_cont)); }

/*
void Graph::unmap_chunk(Read_Chunk_BPtr rc_bptr)
{
    rc_bptr = trim_contig_to_chunk(rc_bptr);
    Read_Entry_BPtr re_bptr = rc_bptr->re_bptr();
    Size_Type rc_start = rc_bptr->get_r_start();
    Size_Type rc_end = rc_bptr->get_r_end();
    Contig_Entry_BPtr ce_bptr = rc_bptr->ce_bptr();
    // chunks mapping to this contig are remapped to individual contigs, with unmappable flag set
    // we also save them in a list using Read_Entry and offsets
    vector< tuple< Read_Entry_BPtr, Size_Type, Size_Type > > l;
    ce_bptr->chunk_cont().clear_and_dispose([&] (Read_Chunk_BPtr other_rc_bptr) {
        l.push_back(make_tuple(other_rc_bptr->re_bptr(), other_rc_bptr->get_r_start(), other_rc_bptr->get_r_end()));
        Read_Chunk::make_unmappable(other_rc_bptr);
        ce_cont().insert(other_rc_bptr->ce_bptr());
    });
    // done with old Contig_Entry
    ce_bptr->mut_cont().clear_and_dispose();
    ce_cont().erase(ce_bptr);
    Contig_Entry_Fact::del_elem(ce_bptr);

    // for all read chunks made unmappable, try to merge them with nearby unmappable chunks
    set< Read_Entry_CBPtr > re_set;
    for (auto& e : l) {
        tie(re_bptr, rc_start, rc_end) = e;
        re_set.insert(re_bptr);
        Read_Chunk_BPtr other_rc_bptr = re_bptr->chunk_cont().get_chunk_with_pos(rc_start).unconst();
        ASSERT(other_rc_bptr->ce_bptr()->is_unmappable());
        ASSERT(other_rc_bptr->get_r_start() <= rc_start and rc_end <= other_rc_bptr->get_r_end());
        if (other_rc_bptr->get_r_start() != rc_start or other_rc_bptr->get_r_end() != rc_end)
        {
            // already been merged with other unmappable chunks
            continue;
        }
        //TODO: this should be done as a separate graph operation
        extend_unmapped_chunk(re_bptr, rc_start, rc_end);
    }
    check(re_set);
} // unmap_chunk
*/

/*
void Graph::unmap_re_region(Read_Entry_BPtr re_bptr, const Range_Type& re_rg)
{
    // region should be non-empty
    ASSERT(re_rg.start() < re_rg.end());
    Size_Type pos = re_rg.start();
    // cut first chunk if necessary
    Read_Chunk_BPtr rc_bptr = re_bptr->chunk_cont().get_chunk_with_pos(re_rg.start()).unconst();
    ASSERT(rc_bptr->get_r_start() <= re_rg.start() and re_rg.start() < rc_bptr->get_r_end());
    if (rc_bptr->is_unbreakable())
    {
        pos = rc_bptr->get_r_end();
    }
    else if (rc_bptr->get_r_start() < re_rg.start())
    {
        cut_read_chunk(rc_bptr, re_rg.start());
        rc_bptr = re_bptr->chunk_cont().get_chunk_with_pos(re_rg.start()).unconst();
        ASSERT(rc_bptr->get_r_start() == re_rg.start());
    }

    // unmap chunks repeatedly until the end of the range
    while (pos < re_rg.end())
    {
        rc_bptr = re_bptr->chunk_cont().get_chunk_with_pos(pos).unconst();
        if (rc_bptr->is_unbreakable())
        {
            ASSERT(rc_bptr->get_r_start() <= pos);
            pos = rc_bptr->get_r_end();
        }
        else
        {
            ASSERT(rc_bptr->get_r_start() == pos);
            // if chunk contains the region end, we cut it first
            if (re_rg.end() < rc_bptr->get_r_end())
            {
                cut_read_chunk(rc_bptr, re_rg.end());
                rc_bptr = re_bptr->chunk_cont().get_chunk_with_pos(pos).unconst();
                ASSERT(rc_bptr->get_r_start() == pos and rc_bptr->get_r_end() == re_rg.end());
            }
            pos = rc_bptr->get_r_end();
            unmap_chunk(rc_bptr);
        }
    }
} // unmap_region
*/

/*
void Graph::extend_unmapped_chunk_dir(Read_Entry_BPtr re_bptr, Size_Type pos, bool r_right)
{
    LOG("graph", debug2) << ptree("extend_unmap_chunk_dir")
        .put("re_name", re_bptr->name())
        .put("re_ptr", re_bptr.to_int())
        .put("pos", pos)
        .put("r_right", r_right);

    bool r_left = not r_right;
    Size_Type leftover_bp = (r_left? pos - re_bptr->start() : re_bptr->end() - pos);
    while (leftover_bp > 0)
    {
        // calls in previous iterations might extend the unmapped region
        // we first recompute the chunks based solely on: re_cptr & unmap_start
        Read_Chunk_BPtr rc_bptr = re_bptr->chunk_cont().get_chunk_with_pos(r_left? pos : pos - 1).unconst();
        ASSERT(rc_bptr->ce_bptr()->is_unmappable());
        ASSERT(not r_left or rc_bptr->get_r_start() <= pos);
        ASSERT(not r_right or pos <= rc_bptr->get_r_end());
        pos = (r_left? rc_bptr->get_r_start() : rc_bptr->get_r_end());
        leftover_bp = (r_left? pos - re_bptr->start() : re_bptr->end() - pos);
        if (leftover_bp == 0)
        {
            break;
        }
        Read_Chunk_BPtr next_rc_bptr = re_bptr->chunk_cont().get_sibling(rc_bptr, true, r_right).unconst();

        if (next_rc_bptr->ce_bptr()->is_unmappable())
        {
            // consecutive unmappable chunks; merge contigs
            ASSERT(not rc_bptr->get_rc());
            ASSERT(not next_rc_bptr->get_rc());
            pos = (r_left? next_rc_bptr->get_r_start() : next_rc_bptr->get_r_end());
            bool success = cat_contigs(r_left? next_rc_bptr->ce_bptr() : rc_bptr->ce_bptr(), true);
            static_cast< void >(success);
            ASSERT(success);
            continue;
        }
        else // next_rc_bptr not unmappable
        {
            // next chunk is mappable
            if (leftover_bp <= unmap_trigger_len())
            {
                unmap_chunk(next_rc_bptr);
                continue;
            }
            else // global::unmap_trigger_len < leftover_bp
            {
                // check there is a minimum of mappable bp
                // if not, unmap small intermediate chunks and merge them
                // skip small mappable chunks
                list< tuple< Size_Type, Size_Type > > skipped_chunks;
                Size_Type skipped_len = 0;
                while (not next_rc_bptr->ce_bptr()->is_unmappable()
                       and skipped_len + next_rc_bptr->get_r_len() <= unmap_trigger_len())
                {
                    skipped_chunks.push_back(make_tuple(next_rc_bptr->get_r_start(), next_rc_bptr->get_r_end()));
                    skipped_len += next_rc_bptr->get_r_len();
                    next_rc_bptr = re_bptr->chunk_cont().get_sibling(next_rc_bptr, true, r_right).unconst();
                    // since skipped_len < unmap_trigger_len < leftover_bp
                    ASSERT(next_rc_bptr);
                }
                ASSERT(skipped_len <= unmap_trigger_len());
                if (next_rc_bptr->ce_bptr()->is_unmappable())
                {
                    // unmap all skipped chunks (if they are not modified in previous iterations)
                    for (auto it = skipped_chunks.begin(); it != skipped_chunks.end(); ++it)
                    {
                        Size_Type rc_start;
                        Size_Type rc_end;
                        tie(rc_start, rc_end) = *it;
                        rc_bptr = re_bptr->chunk_cont().get_chunk_with_pos(rc_start).unconst();
                        if (rc_bptr->get_r_start() != rc_start or rc_bptr->get_r_end() != rc_end or rc_bptr->ce_bptr()->is_unmappable())
                        {
                            // already been merged with other unmappable chunks
                            continue;
                        }
                        unmap_chunk(rc_bptr);
                    }
                }
                else // next_rc_bptr not unmappable
                {
                    break;
                }
            } // else (global::unmap_trigger_len < leftover_bp)
        } // else (next_rc_bptr not unmappable)
    } // while (leftover_bp > 0)

    // if the last chunk is unmappable, trim read
    Read_Chunk_BPtr last_rc_bptr = (not r_right? &*re_bptr->chunk_cont().begin() : &*re_bptr->chunk_cont().rbegin());
    if (last_rc_bptr and last_rc_bptr->ce_bptr()->is_unmappable())
    {
        Contig_Entry_BPtr last_ce_bptr = last_rc_bptr->ce_bptr();
        // trim read
        re_bptr->chunk_cont().erase(last_rc_bptr);
        if (not r_right)
        {
            re_bptr->start() += last_rc_bptr->get_r_len();
        }
        else
        {
            re_bptr->len() -= last_rc_bptr->get_r_len();
        }
        // destroy RC & CE
        ASSERT(last_ce_bptr->mut_cont().empty());
        ASSERT(last_rc_bptr->mut_ptr_cont().empty());
        last_ce_bptr->chunk_cont().erase(last_rc_bptr);
        Read_Chunk_Fact::del_elem(last_rc_bptr);
        ce_cont().erase(last_ce_bptr);
        Contig_Entry_Fact::del_elem(last_ce_bptr);
    }

    check(set< Read_Entry_CBPtr >( { re_bptr }));
} // extend_unmapped_chunk_dir

void Graph::extend_unmapped_chunk(Read_Entry_BPtr re_bptr, Size_Type rc_start, Size_Type rc_end)
{
    extend_unmapped_chunk_dir(re_bptr, rc_start, false);
    extend_unmapped_chunk_dir(re_bptr, rc_end, true);
}
*/

Range_Cont
Graph::find_unmappable_regions(Read_Entry_CBPtr re_cbptr, Size_Type r_start, Size_Type r_end) const
{
    Range_Cont region_cont;
    Size_Type pos = r_start;
    while (pos < r_end)
    {
        Read_Chunk_CBPtr rc_cbptr = re_cbptr->chunk_cont().get_chunk_with_pos(pos);
        ASSERT(rc_cbptr);
        ASSERT(rc_cbptr->get_r_start() <= pos and pos < rc_cbptr->get_r_end());
        pos = rc_cbptr->get_r_end();
        find_unmappable_regions(rc_cbptr, region_cont);
    }
    return region_cont;
}

void
Graph::find_unmappable_regions(Read_Chunk_CBPtr orig_rc_cbptr, Range_Cont& region_cont) const
{
    // If different chunks of the same read
    // are mapped to overlapping regions of the contig,
    // we make those regions unmappable
    // By property of contig entries, 2 chunks of the same read
    // can be mapped to the same contig only if each overlaps
    // at least one contig endpoint.

    Contig_Entry_CBPtr ce_cbptr = orig_rc_cbptr->ce_bptr();
    // group all extremal read chunks by their read entry
    map< Read_Entry_CBPtr, set< Read_Chunk_CBPtr > > re_chunk_map;
    // first look at contig start
    for (auto rc_cbptr : ce_cbptr->chunk_cont().iintersect(0, 0) | referenced)
    {
        re_chunk_map[rc_cbptr->re_bptr()].insert(rc_cbptr);
    }
    // then at contig end
    for (auto rc_cbptr : ce_cbptr->chunk_cont().iintersect(ce_cbptr->len(), ce_cbptr->len()) | referenced)
    {
        re_chunk_map[rc_cbptr->re_bptr()].insert(rc_cbptr);
    }

    // read regions to unmap
    for (auto& t : re_chunk_map)
    {
        //auto& re_bptr = t.first;
        auto& rc_bptr_cont = t.second;
        if (rc_bptr_cont.size() == 1)
        {
            // ignore reads with single extremal chunk
            continue;
        }
        // consider all pairwise overlaps at this point
        for (auto it_1 = rc_bptr_cont.begin(); it_1 != rc_bptr_cont.end(); ++it_1)
        {
            Read_Chunk_CBPtr rc1_cbptr = *it_1;
            for (auto it_2 = rc_bptr_cont.begin(); it_2 != it_1; ++it_2)
            {
                Read_Chunk_CBPtr rc2_cbptr = *it_2;
                Range_Type overlap_c_rg(max(rc1_cbptr->get_c_start(), rc2_cbptr->get_c_start()),
                                        min(rc1_cbptr->get_c_end(), rc2_cbptr->get_c_end()));
                // restrict to original chunk coordinates
                overlap_c_rg.begin() = max(overlap_c_rg.begin(), orig_rc_cbptr->get_c_start());
                overlap_c_rg.end() = min(overlap_c_rg.end(), orig_rc_cbptr->get_c_end());
                if (overlap_c_rg.end() < overlap_c_rg.begin())
                {
                    // negative overlap region; ignore
                    continue;
                }
                // overlap_c_rg.start() <= overlap_c_rg.end()
                // compute read positions on original chunk
                auto rg = orig_rc_cbptr->mapped_range(overlap_c_rg, true, true, true);
                region_cont.insert(rg);
            }
        }
    }
    /*
    auto re_it = re_list.begin();
    auto pos_it = pos_list.begin();
    for (; re_it != re_list.end() and pos_it != pos_list.end(); ++re_it, ++pos_it)
    {
        const Read_Entry* re_cptr;
        Size_Type rc_start;
        Size_Type rc_end;
        re_cptr = *re_it;
        tie(rc_start, rc_end) = *pos_it;
        Read_Chunk_CPtr rc_cptr = re_cptr->get_chunk_with_pos(rc_start);
        // either the chunk has the same boundaries as when it was discovered
        // or it was marked unmappable in previous iterations
        //ASSERT((rc_cptr->get_r_start() == rc_start and rc_cptr->get_r_end() == rc_end) or rc_cptr->is_unmappable());
        if (not rc_cptr->is_unmappable())
        {
            unmap_chunk(rc_cptr);
        }
    }
    ASSERT(check(set< const Read_Entry* >(re_list.begin(), re_list.end())));
    */
}

void Graph::unmap_single_chunks()
{
    LOG("graph", info) << ptree("begin");
    for (auto re_bptr : _re_cont | referenced)
    {
        bool done = false;
        while (not done)
        {
            done = true;
            for (auto rc_bptr : re_bptr->chunk_cont() | referenced)
            {
                Contig_Entry_BPtr ce_bptr = rc_bptr->ce_bptr();
                LOG("graph", debug) << ptree("unmap_single_chunks_loop")
                    .put("ce_ptr", ce_bptr.to_int())
                    .put("is_unmappable", ce_bptr->is_unmappable())
                    .put("size_one", size_one(ce_bptr->chunk_cont()));
                if (not ce_bptr->is_unmappable() and size_one(ce_bptr->chunk_cont()))
                {
                    LOG("graph", debug) << ptree("unmap_single_chunks_unmap_start").put("re_ptr", re_bptr.to_int());
                    unmap_chunk(rc_bptr);
                    LOG("graph", debug) << ptree("unmap_single_chunks_unmap_end").put("re_ptr", re_bptr.to_int());
                    done = false;
                    break;
                }
            }
        }
    }
    LOG("graph", info) << ptree("end");
}

void Graph::unmap_read_ends()
{
    LOG("graph", info) << ptree("begin");
    for (auto re_bptr : re_cont() | referenced)
    {
        unmap_single_terminal_chunk(&*re_bptr->chunk_cont().begin(), true);
        unmap_single_terminal_chunk(&*re_bptr->chunk_cont().rbegin(), false);
    }
    LOG("graph", info) << ptree("end");
}

void Graph::unmap_single_terminal_chunk(Read_Chunk_BPtr rc_bptr, bool r_start)
{
    ASSERT(not rc_bptr->re_bptr()->chunk_cont().get_sibling(rc_bptr, true, not r_start));
    if (rc_bptr->ce_bptr()->is_unmappable())
    {
        return;
    }
    if (size_one(rc_bptr->ce_bptr()->chunk_cont()))
    {
        unmap_chunk(rc_bptr);
        return;
    }
    Read_Entry_BPtr re_bptr = rc_bptr->re_bptr();
    Contig_Entry_BPtr ce_bptr = rc_bptr->ce_bptr();
    bool c_start = (r_start != rc_bptr->get_rc());
    if (c_start)
    {
        if (// chunk does not span c_start
            rc_bptr->get_c_start() > 0
            // or it's not the first in the chunk container
            or &*ce_bptr->chunk_cont().begin() != rc_bptr
            // or a second chunk exists and starts at c_start
            or (not size_one(ce_bptr->chunk_cont()) and (++ce_bptr->chunk_cont().begin())->get_c_start() == 0))
        {
            return;
        }
        Size_Type c_brk = (++ce_bptr->chunk_cont().begin())->get_c_start();
        // there can be no mutations on rc_bptr before c_brk
        ASSERT(rc_bptr->mut_ptr_cont().empty()
               or rc_bptr->mut_ptr_cont().begin()->mut_cbptr()->rf_start() >= c_brk);
        cut_contig_entry(ce_bptr, c_brk, nullptr);
        rc_bptr = r_start? &*re_bptr->chunk_cont().begin() : &*re_bptr->chunk_cont().rbegin();
        unmap_chunk(rc_bptr);
    }
    else
    {
        if (// chunk does not span c_end
            rc_bptr->get_c_end() < ce_bptr->len())
        {
            return;
        }
        // this is trickier than the c_start==true case
        // because chunks are not in the order of their end pos
        auto rg = ce_bptr->chunk_cont().iintersect(ce_bptr->len(), ce_bptr->len());
        // at least rc_bptr must span c_end
        ASSERT(rg.begin() != rg.end());
        if (++rg.begin() != rg.end())
        {
            // more than 2 chunks span c_end
            return;
        }
        ASSERT(&*rg.begin() == rc_bptr);
        Size_Type c_brk = ce_bptr->chunk_cont().max_end(ce_bptr->len() - 1);
        // max_end smaller than contig end must exist because chunk_cont.size() >= 2 and rg.size() == 1
        ASSERT(c_brk < ce_bptr->len());
        // there can be no mutations in rc_bptr past c_brk
        ASSERT(rc_bptr->mut_ptr_cont().empty()
               or rc_bptr->mut_ptr_cont().rbegin()->mut_cbptr()->rf_end() <= c_brk);
        Mutation_CBPtr mut_left_cbptr = nullptr;
        if (not rc_bptr->mut_ptr_cont().empty()
            and rc_bptr->mut_ptr_cont().rbegin()->mut_cbptr()->is_ins()
            and rc_bptr->mut_ptr_cont().rbegin()->mut_cbptr()->rf_start() == c_brk)
        {
            mut_left_cbptr = rc_bptr->mut_ptr_cont().rbegin()->mut_cbptr();
        }
        cut_contig_entry(ce_bptr, c_brk, mut_left_cbptr);
        rc_bptr = r_start? &*re_bptr->chunk_cont().begin() : &*re_bptr->chunk_cont().rbegin();
        unmap_chunk(rc_bptr);
    }
}

void Graph::print_mutations(ostream& os, size_t min_support, Size_Type flank_len) const
{
    os << "MUT\tid.ce\tlen.ce\tpos.rf\tlen.rf\tlen.alt\tseq.rf\tseq.alt\t"
        "flank.left\tflank.right\tnum.cks.rf\tnum.cks.alt\t"
        "uniq.rf\tuniq.alt\tcn.rf\tcn.alt\n";
    for (auto ce_bptr : ce_cont() | referenced |
             ba::filtered([] (Contig_Entry_CBPtr ce_cbptr) { return ce_cbptr->is_normal(); }))
    {
        for (auto mut_bptr : ce_bptr->mut_cont() | referenced)
        {
            size_t num_cks_qr = mut_bptr->chunk_ptr_cont().nonconst_size();
            size_t num_cks_total = boost::distance(ce_bptr->chunk_cont().iintersect(mut_bptr->rf_start(), mut_bptr->rf_end()));
            ASSERT(num_cks_qr <= num_cks_total);
            size_t num_cks_rf = num_cks_total - num_cks_qr;
            if (min(num_cks_rf, num_cks_qr) < min_support)
            {
                continue;
            }
            os << "MUT\t" << ce_bptr.to_int() << "\t" << ce_bptr->len() << "\t"
               << mut_bptr->rf_start() << "\t" << mut_bptr->rf_len() << "\t" << mut_bptr->seq_len() << "\t"
               << ce_bptr->seq().substr(mut_bptr->rf_start(), mut_bptr->rf_len()) << "\t"
               << mut_bptr->seq() << "\t"
               << ce_bptr->seq().substr(mut_bptr->rf_start() >= flank_len? mut_bptr->rf_start() - flank_len : 0,
                                        mut_bptr->rf_start() >= flank_len? flank_len : mut_bptr->rf_start()) << "\t"
               << ce_bptr->seq().substr(mut_bptr->rf_end(),
                                        min(ce_bptr->len() - mut_bptr->rf_end(), flank_len)) << "\t"
               << num_cks_rf << "\t" << num_cks_qr << "\t"
               << (int)mut_bptr->uniq(0) << "\t" << (int)mut_bptr->uniq(1) << "\t"
               << (int)mut_bptr->copy_num(0) << "\t" << (int)mut_bptr->copy_num(1) << endl;
        }
    }
}

void Graph::print_unmappable_contigs(ostream& os) const
{
    for (auto ce_cbptr : ce_cont() | referenced)
    {
        if (ce_cbptr->is_unmappable())
        {
            continue;
        }
        for (int dir = 0; dir < 2; ++dir)
        {
            bool c_right = (dir == 1);
            auto cks = ce_cbptr->out_chunks_dir(c_right, 3);
            for (auto& t : cks)
            {
                Contig_Entry_CBPtr ce_next_cbptr;
                bool same_orientation;
                set< Read_Chunk_CBPtr > chunk_v;
                tie(ce_next_cbptr, same_orientation) = move(t.first);
                chunk_v = move(t.second);
                multiset< tuple< string, string > > seq_v;
                if (ce_next_cbptr.to_int() < ce_cbptr.to_int())
                {
                    continue;
                }
                for (auto& rc_cbptr : chunk_v)
                {
                    Read_Chunk_CBPtr rc_next_cbptr =
                        rc_cbptr->re_bptr()->chunk_cont().get_sibling(rc_cbptr, false, c_right);
                    ASSERT(rc_next_cbptr);
                    if (rc_next_cbptr->ce_bptr()->is_unmappable())
                    {
                        seq_v.insert(make_tuple(
                                         rc_next_cbptr->get_seq().revcomp(rc_cbptr->get_rc()),
                                         rc_cbptr->re_bptr()->name()));
                    }
                }
                if (not seq_v.empty())
                {
                    os << ">\t" << ce_cbptr.to_int() << "\t" << dir << "\t"
                       << ce_next_cbptr.to_int() << "\t" << static_cast< int >(same_orientation) << "\n";
                    for (const auto& t : seq_v)
                    {
                        os << get<1>(t) << ": " << get<0>(t) << "\n";
                    }
                }
            }
        }
    }
}

void Graph::resolve_unmappable_fully_mapped(
    const map< Seq_Type, size_t >& seq_cnt_map,
    vector< Seq_Type >& bseq_v,
    map< Seq_Type, tuple< size_t, string > >& seq_bseq_map)
{
    ASSERT(bseq_v.empty());
    ASSERT(seq_bseq_map.empty());
    /*
    while (seq_bseq_map.size() < seq_cnt_map.size())
    {
        // find most popular and longest sequences not yet mapped
        string most_pop_candidate_seq;
        size_t most_pop_candidate_seq_count = 0;
        string longest_candidate_seq;
        for (const auto& t : seq_cnt_map)
        {
            const string& seq = t.first;
            const size_t& seq_count = t.second;
            if (seq_bseq_map.count(seq) == 0)
            {
                // seq to yet mapped to a base sequence
                if (seq_count > most_pop_candidate_seq_count)
                {
                    most_pop_candidate_seq = seq;
                    most_pop_candidate_seq_count = seq_count;
                }
                if (seq.size() > longest_candidate_seq.size())
                {
                    longest_candidate_seq = seq;
                }
            }
        }
        // HEURISTIC for new base sequence:
        // - use most popular unmapped sequence if it occurs a minimum number of times
        // - if not, use longest unmapped sequence
        static const size_t min_count_new_base_seq = 2;
        if (most_pop_candidate_seq_count >= min_count_new_base_seq)
        {
            LOG("graph", "debug1") << ptree("new_base_seq").put("seq", most_pop_candidate_seq);
            bseq_v.push_back(most_pop_candidate_seq);
        }
        else
        {
            LOG("graph", "debug1") << ptree("new_base_seq").put("seq", longest_candidate_seq);
            bseq_v.push_back(longest_candidate_seq);
        }
        // found new candidate sequence
        // remap all sequences to all candidate sequences
        seq_bseq_map.clear();
        for (const auto& t : seq_cnt_map)
        {
            const string& seq = t.first;
            size_t bseq_idx = bseq_v.size();
            string cigar_string;
            long max_score = numeric_limits< long >::min();
            for (size_t i = 0; i < bseq_v.size(); ++i)
            {
                // map seq to bseq
                OverlapperParams params = affine_default_params;
                params.type = ALT_GLOBAL;
                params.use_m_ops = false;
                auto output = Overlapper::computeAlignmentAffine(bseq_v[i], seq, params);
                if (output.score > max_score)
                {
                    // HEURISTIC for allowing this mapping:
                    // - (edit distance < 5) OR (score >= 70% * max_score)
                    if (output.edit_distance < 5
                        or 100 * output.score >= 70 * params.match_score * static_cast< long >(min(seq.size(), bseq_v[i].size())))
                    {
                        bseq_idx = i;
                        cigar_string = output.cigar;
                        max_score = output.score;
                    }
                }
            }
            if (bseq_idx < bseq_v.size())
            {
                seq_bseq_map[seq] = make_tuple(bseq_idx, cigar_string);
            }
        } // for (seq : seq_cnt_map)
    } // while not all mapped
    */

    // move sequences into a vector, ignore counts
    vector< Seq_Type > seq_v;
    for (const auto& t : seq_cnt_map)
    {
        LOG("graph", debug1) << ptree("begin")
            .put("rseq", t.first);
        seq_v.emplace_back(t.first);
    }

    // for every sequence, count 3-mers and save them in a map
    vector< map< Seq_Type, size_t > > kmer_map_v(seq_v.size());
    for (size_t k = 0; k < seq_v.size(); ++k)
    {
        const Seq_Type& seq = seq_v[k];
        map< Seq_Type, size_t >& m = kmer_map_v[k];
        for (size_t i = 0; i + 3 < seq.size(); ++i)
        {
            Seq_Type kmer = Seq_Type(seq.substr(i, 3));
            if (m.count(kmer) == 0)
            {
                m[kmer] = 1;
            }
            else
            {
                ++m[kmer];
            }
        }
    }

    // count all pairwise distances
    auto get_dist = [] (const map< Seq_Type, size_t >& km1, const map< Seq_Type, size_t >& km2)
    {
        // we compute the distance as follows:
        //   sum_{kmer k} |count(s1,k) - count(s2,k)| / max(count(s1,k), count(s2,k)) / sum(denom)
        double res = 0.0;
        double denom = 0.0;
        for (const auto& k : km1 | ba::map_keys)
        {
            double cnt1 = static_cast< double >(km1.at(k));
            double cnt2 = 0.0;
            if (km2.count(k) > 0)
            {
                cnt2 = km2.at(k);
            }
            LOG("graph", debug2) << ptree("contrib")
                .put("kmer", k)
                .put("cnt1", cnt1)
                .put("cnt2", cnt2);
            res += abs(cnt1 - cnt2) / max(cnt1, cnt2);
            denom += max(cnt1, cnt2);
        }
        for (const auto& k : km2 | ba::map_keys)
        {
            if (km1.count(k) == 0)
            {
                res += 1.0;
            }
            denom += static_cast< double >(km2.at(k));
        }
        return res / (denom > 0.0? denom : 1.0);
    };
    vector< vector< double > > dist(seq_v.size(), vector< double >(seq_v.size()));
    for (size_t i = 0; i < seq_v.size(); ++i)
    {
        dist[i][i] = 0;
        for (size_t j = i + 1; j < seq_v.size(); ++j)
        {
            LOG("graph", debug1) << ptree("get_dist")
                .put("seq1", seq_v[i])
                .put("seq2", seq_v[j]);
            dist[i][j] = get_dist(kmer_map_v[i], kmer_map_v[j]);
            dist[j][i] = dist[i][j];
            LOG("graph", debug1) << ptree("get_dist_res")
                .put("res", dist[i][j]);
        }
    }

    // pick base sequences:
    //   greedily choose the one closer than thres_dist to most others
    //   repeat as long as there exists a sequence st:
    //   - more than thres_dist to every chosen one AND
    //   - max(length(seq),length(bseq)) > thres_len
    static const double thres_dist = .4;
    static const size_t thres_len = 10;
    OverlapperParams params = affine_default_params;
    params.type = ALT_GLOBAL;
    params.use_m_ops = false;
    while (seq_bseq_map.size() < seq_v.size())
    {
        size_t idx = seq_v.size();
        size_t max_cluster_size = 0;
        for (size_t i = 0; i < seq_v.size(); ++i)
        {
            if (seq_bseq_map.count(seq_v[i]) > 0)
            {
                // already mapped
                continue;
            }
            size_t cluster_size = 1;
            for (size_t j = i + 1; j < seq_v.size(); ++j)
            {
                if (seq_bseq_map.count(seq_v[j]) > 0)
                {
                    continue;
                }
                // neither i nor j are mapped
                if (dist[i][j] < thres_dist)
                {
                    ++cluster_size;
                }
            }
            if (cluster_size > max_cluster_size)
            {
                max_cluster_size = cluster_size;
                idx = i;
            }
        }
        ASSERT(idx < seq_v.size());
        // pick sequence at idx as a new base sequence
        LOG("graph", debug1) << ptree("new_base_seq").put("seq", seq_v[idx]);
        bseq_v.emplace_back(seq_v[idx]);
        const Seq_Type& bseq = bseq_v.back();
        const size_t bseq_idx = bseq_v.size() - 1;
        // map sequences not already mapped to the new base sequence
        for (size_t i = 0; i < seq_v.size(); ++i)
        {
            if (seq_bseq_map.count(seq_v[i]) > 0
                or (dist[idx][i] >= thres_dist and max(seq_v[i].size(), bseq.size()) > thres_len))
            {
                // already mapped or too far
                continue;
            }
            const string& seq = seq_v[i];
            auto output = Overlapper::computeAlignmentAffine(bseq, seq, params);
            seq_bseq_map[seq] = make_tuple(bseq_idx, output.cigar);
        }
    } // while (seq_bseq_map.size() < seq_v.size())
} // Graph::get_base_sequences

void Graph::resolve_unmappable_inner_region(
    Contig_Entry_CBPtr ce_cbptr, bool c_right,
    Contig_Entry_CBPtr ce_next_cbptr, bool same_orientation)
{
    LOG("graph", debug) << ptree("begin")
        .put("ce_ptr", ce_cbptr.to_int())
        .put("c_right", c_right)
        .put("ce_next_ptr", ce_next_cbptr.to_int())
        .put("same_orientation", same_orientation);

    auto oc_map = ce_cbptr->out_chunks_dir(c_right, 3);
    auto t = make_pair(ce_next_cbptr, same_orientation);
    ASSERT(oc_map.count(t) > 0);
    const auto& out_chunks_v = oc_map.at(t);
    // collect unmappable chunks, along with their orientation, in unmappable_chunks_v
    vector< tuple< Read_Chunk_BPtr, bool > > unmappable_chunks_v;
    // collect sequences along with occurrence count in seq_map
    map< Seq_Type, size_t > seq_cnt_map;
    for (const auto& rc_cbptr : out_chunks_v)
    {
        bool r_right = (c_right != rc_cbptr->get_rc());
        Read_Chunk_CBPtr rc_next_cbptr = rc_cbptr->re_bptr()->chunk_cont().get_sibling(rc_cbptr, true, r_right);
        ASSERT(rc_next_cbptr);
        if (not rc_next_cbptr->ce_bptr()->is_unmappable())
        {
            continue;
        }
        unmappable_chunks_v.push_back(make_tuple(rc_next_cbptr.unconst(), r_right));
        Seq_Type seq = rc_next_cbptr->get_seq().revcomp(not r_right);
        if (seq_cnt_map.count(seq) > 0)
        {
            ++seq_cnt_map.at(seq);
        }
        else
        {
            seq_cnt_map[seq] = 1;
        }
        // check next next exists and is mapped to ce_next
        Read_Chunk_CBPtr rc_next_next_cbptr = rc_next_cbptr->re_bptr()->chunk_cont().get_sibling(rc_next_cbptr, true, r_right);
        static_cast< void >(rc_next_next_cbptr); // silence unused warnings when assertions disabled
        ASSERT(rc_next_next_cbptr);
        ASSERT(not rc_next_next_cbptr->ce_bptr()->is_unmappable());
        ASSERT(rc_next_next_cbptr->ce_bptr() == ce_next_cbptr);
        ASSERT((rc_next_next_cbptr->get_rc() == rc_cbptr->get_rc()) == same_orientation);
    }
    // construct set of base sequences
    map< Seq_Type, tuple< size_t, string > > seq_bseq_map;
    vector< Seq_Type > bseq_v;
    resolve_unmappable_fully_mapped(seq_cnt_map, bseq_v, seq_bseq_map);
    // create new contig entries
    vector< Contig_Entry_BPtr > new_ce_bptr_v;
    for (auto& bseq : bseq_v)
    {
        Contig_Entry_BPtr new_ce_bptr = Contig_Entry_Fact::new_elem(move(bseq));
        new_ce_bptr_v.push_back(new_ce_bptr);
        ce_cont().insert(new_ce_bptr);
    }
    for (size_t i = 0; i < unmappable_chunks_v.size(); ++i)
    {
        // create new read chunk mapped to new contig entry
        Read_Chunk_CBPtr old_rc_bptr;
        bool r_right;
        tie(old_rc_bptr, r_right) = unmappable_chunks_v[i];
        // seq as considered during base seq selection
        Seq_Type seq = old_rc_bptr->get_seq().revcomp(not r_right);
        size_t bseq_idx;
        string cigar_string;
        tie(bseq_idx, cigar_string) = seq_bseq_map.at(seq);
        // cigar object that takes into account orientation and offset
        Cigar cigar(cigar_string, not r_right, 0, old_rc_bptr->get_r_start());
        // since the cigar object knows the orientation, we pass the non-reversed chunk sequence
        Read_Chunk_BPtr new_rc_bptr = Read_Chunk::make_chunk_from_cigar(cigar, old_rc_bptr->get_seq(), new_ce_bptr_v[bseq_idx]);
        // set lowcomplexity flag if the contig entry contains more than 1 mutation
        if (new_rc_bptr->ce_bptr()->mut_cont().nonconst_size() > 1)
        {
            new_rc_bptr->ce_bptr()->set_lowcomplex();
        }
        Read_Entry_BPtr re_bptr = old_rc_bptr->re_bptr();
        new_rc_bptr->re_bptr() = re_bptr;
        // unlink & deallocate old unmappable read chunk
        Contig_Entry_BPtr old_ce_bptr = old_rc_bptr->ce_bptr();
        re_bptr->chunk_cont().erase(old_rc_bptr);
        old_ce_bptr->chunk_cont().erase(old_rc_bptr);
        Read_Chunk_Fact::del_elem(old_rc_bptr);
        // unlink & deallocate old unmappable contig entry
        ce_cont().erase(old_ce_bptr);
        Contig_Entry_Fact::del_elem(old_ce_bptr);
        // link new chunk into its read entry container
        re_bptr->chunk_cont().insert(new_rc_bptr);
    }
    check(set< Contig_Entry_CBPtr >(new_ce_bptr_v.begin(), new_ce_bptr_v.end()));
}

void Graph::resolve_unmappable_partially_mapped(
    const set< tuple< string, bool > >& seq_set,
    const vector< string >& bseq_v,
    map< tuple< string, bool >, tuple< size_t, string > >& seq_bseq_map)
{
    ASSERT(seq_bseq_map.empty());
    ASSERT(not bseq_v.empty());
    for (const auto& t : seq_set)
    {
        const string& seq = get<0>(t);
        const bool bseq_right = get<1>(t);
        size_t bseq_idx = bseq_v.size();
        string cigar_string;
        long max_score = numeric_limits< long >::min();
        OverlapperParams params = affine_default_params;
        params.type = ALT_CUSTOM;
        params.use_m_ops = false;
        params.gap_s2_start = false;
        params.gap_s2_end = false;
        for (size_t i = 0; i < bseq_v.size(); ++i)
        {
            params.gap_s1_start = bseq_right;
            params.gap_s1_end = not bseq_right;
            const auto output = Overlapper::computeAlignmentAffine(bseq_v[i], seq, params);
            if (output.score > max_score)
            {
                // HEURISTIC for allowing this mapping:
                // - (edit distance < 5) OR (score >= 70% * max_score)
                const long perfect_score = params.match_score * static_cast< long >(min(seq.size(), bseq_v[i].size()));
                if (output.edit_distance < 10
                    or 100 * output.score >= 60 * perfect_score)
                {
                    bseq_idx = i;
                    cigar_string = output.cigar;
                    max_score = output.score;
                }
            }
        }
        if (bseq_idx < bseq_v.size())
        {
            seq_bseq_map[t] = make_tuple(bseq_idx, cigar_string);
        }
    }
}

void Graph::resolve_unmappable_terminal_region(Contig_Entry_CBPtr ce_cbptr, bool c_right)
{
    LOG("graph", debug) << ptree("begin")
        .put("ce_ptr", ce_cbptr.to_int())
        .put("c_right", c_right);

    auto out_chunks_map = ce_cbptr->out_chunks_dir(c_right, 4);
    // this bucket contains the last chunks which are followed by terminal unmappable siblings
    const auto rc_last_bucket = make_pair(Contig_Entry_CBPtr(), false);
    const auto& rc_last_cbptr_v = out_chunks_map.at(rc_last_bucket);
    ASSERT(not rc_last_cbptr_v.empty());
    // group base sequences of the neighbouring contigs
    vector< tuple< Contig_Entry_CBPtr, bool > > ce_next_v;
    vector< string > ce_next_seq_v;
    for (const auto& t : out_chunks_map)
    {
        Contig_Entry_CBPtr ce_next_cbptr;
        bool same_orientation;
        tie(ce_next_cbptr, same_orientation) = t.first;
        if (not ce_next_cbptr or ce_next_cbptr->seq().empty())
        {
            continue;
        }
        ce_next_v.push_back(t.first);
        //ce_next_seq_v.push_back(same_orientation? ce_next_cbptr->seq() : reverseComplement(ce_next_cbptr->seq()));
        ce_next_seq_v.push_back(string(ce_next_cbptr->seq().revcomp(not same_orientation)));
        LOG("graph", debug1) << ptree("bseq")
            .put("bseq", ce_next_seq_v.back());
    }
    if (ce_next_v.empty())
    {
        // no mappable contigs beyond this one
        return;
    }
    // group unmappable chunk seqeunces
    set< tuple< string, bool > > rc_unmap_seq_set;
    vector< tuple< Read_Chunk_CBPtr, bool, decltype(rc_unmap_seq_set)::const_iterator > > rc_unmap_v;
    for (const auto& rc_last_cbptr : rc_last_cbptr_v)
    {
        bool r_right = (c_right != rc_last_cbptr->get_rc());
        Read_Chunk_CBPtr rc_unmap_bptr = rc_last_cbptr->re_bptr()->chunk_cont().get_sibling(rc_last_cbptr, true, r_right);
        ASSERT(rc_unmap_bptr);
        ASSERT(rc_unmap_bptr->ce_bptr()->is_unmappable());
        auto t = rc_unmap_seq_set.insert(
            make_tuple(
                rc_unmap_bptr->get_seq().revcomp(rc_last_cbptr->get_rc()),
                not c_right));
        rc_unmap_v.push_back(make_tuple(rc_unmap_bptr, rc_last_cbptr->get_rc(), t.first));
        LOG("graph", debug1) << ptree("rseq")
            .put("rseq", get<0>(*get<2>(rc_unmap_v.back())));
    }
    map< tuple< string, bool >, tuple< size_t, string > > seq_bseq_map;
    resolve_unmappable_partially_mapped(rc_unmap_seq_set, ce_next_seq_v, seq_bseq_map);
    for (size_t i = 0; i < rc_unmap_v.size(); ++i)
    {
        Read_Chunk_CBPtr rc_cbptr;
        bool r_strand;
        decltype(rc_unmap_seq_set)::const_iterator set_cit;
        tie(rc_cbptr, r_strand, set_cit) = rc_unmap_v[i];
        const auto& t = *set_cit;
        if (seq_bseq_map.count(t) == 0)
        {
            LOG("graph", debug) << ptree("drop")
                .put("dropping_rseq", get<0>(t));
            //Read_Chunk_BPtr rc_bptr = rc_cbptr.unconst();
            //rc_bptr->ce_bptr()->set_ignored(); //TODO
            continue;
        }
        size_t bseq_idx;
        string cigar_string;
        tie(bseq_idx, cigar_string) = seq_bseq_map[t];
        Contig_Entry_CBPtr ce_cbptr;
        bool same_orientation;
        tie(ce_cbptr, same_orientation) = ce_next_v[bseq_idx];
        bool c_strand = not same_orientation;
        Cigar cigar = Cigar(cigar_string, c_strand != r_strand, 0, rc_cbptr->get_r_start(), c_strand);
        // shift cigar if it would be aligned to the end of the contig
        if ((not c_strand and not c_right) or (c_strand and c_right))
        {
            cigar.rf_start() += (ce_cbptr->len() - cigar.rf_len());
        }
        LOG("graph", debug1) << ptree("new_mapping")
            .put("rseq", get<0>(t))
            .put("bseq", ce_next_seq_v[bseq_idx])
            .put("cigar", cigar_string);
        Contig_Entry_BPtr ce_bptr = ce_cbptr.unconst();
        Read_Chunk_BPtr rc_bptr = rc_cbptr.unconst();
        // since the cigar object knows the orientation, we pass the non-reversed chunk sequence
        Read_Chunk_BPtr new_rc_bptr = Read_Chunk::make_chunk_from_cigar(cigar, rc_bptr->get_seq(), ce_bptr);
        Read_Entry_BPtr re_bptr = rc_bptr->re_bptr();
        new_rc_bptr->re_bptr() = re_bptr;
        // unlink & deallocate old unmappable read chunk
        Contig_Entry_BPtr old_ce_bptr = rc_bptr->ce_bptr();
        re_bptr->chunk_cont().erase(rc_bptr);
        old_ce_bptr->chunk_cont().erase(rc_bptr);
        Read_Chunk_Fact::del_elem(rc_bptr);
        // unlink & deallocate old unmappable contig entry
        ce_cont().erase(old_ce_bptr);
        Contig_Entry_Fact::del_elem(old_ce_bptr);
        // link new chunk into its read entry container
        re_bptr->chunk_cont().insert(new_rc_bptr);
    }
}

void Graph::resolve_unmappable_regions()
{
    LOG("graph", info) << ptree("begin");
    // step 1: find reads with unmappable contigs in between mappable ones
    for (auto re_bptr : re_cont() | referenced)
    {
        for (auto rc_bref_it = re_bptr->chunk_cont().begin(); rc_bref_it != re_bptr->chunk_cont().end(); )
        {
            Read_Chunk_CBPtr rc_cbptr = &*(rc_bref_it++); // loop it incremented here, before modifications to current elem
            if (not rc_cbptr->ce_bptr()->is_unmappable())
            {
                continue;
            }
            Read_Chunk_CBPtr rc_prev_cbptr = re_bptr->chunk_cont().get_sibling(rc_cbptr, true, false);
            Read_Chunk_CBPtr rc_next_cbptr = re_bptr->chunk_cont().get_sibling(rc_cbptr, true, true);
            if (not rc_prev_cbptr or not rc_next_cbptr)
            {
                continue;
            }
            ASSERT(not rc_prev_cbptr->ce_bptr()->is_unmappable());
            ASSERT(not rc_next_cbptr->ce_bptr()->is_unmappable());
            resolve_unmappable_inner_region(rc_prev_cbptr->ce_bptr(), not rc_prev_cbptr->get_rc(),
                                            rc_next_cbptr->ce_bptr(), rc_next_cbptr->get_rc() == rc_prev_cbptr->get_rc());
        }
    }
    check_all();

    // step 2: find reads with terminal unmappable contigs
    const vector< uint64_t > visit_endpoint_mask = { 1, 2 };
    for (auto ce_bptr : ce_cont() | referenced)
    {
        bitmask::reset(ce_bptr->tag(), visit_endpoint_mask[0] | visit_endpoint_mask[1]);
    }
    for (auto re_bptr : re_cont() | referenced)
    {
        for (int r_right = 0; r_right < 2; ++r_right)
        {
            Read_Chunk_CBPtr rc_last_cbptr = (not r_right? &*re_bptr->chunk_cont().begin() : &*re_bptr->chunk_cont().rbegin());
            if (not rc_last_cbptr->ce_bptr()->is_unmappable())
            {
                // last chunk in this dir is mappable
                continue;
            }
            Read_Chunk_CBPtr rc_cbptr = re_bptr->chunk_cont().get_sibling(rc_last_cbptr, true, not r_right);
            if (not rc_cbptr)
            {
                // entire read is unmappable
                continue;
            }
            Contig_Entry_BPtr ce_bptr = rc_cbptr->ce_bptr().unconst();
            bool c_right = (r_right != rc_cbptr->get_rc());
            if (not bitmask::any(ce_bptr->tag(), visit_endpoint_mask[c_right]))
            {
                bitmask::set(ce_bptr->tag(), visit_endpoint_mask[c_right]);
                resolve_unmappable_terminal_region(ce_bptr, c_right);
            }
        }
    }
    check_all();

    // step 3: cat all read contigs
    cat_all_read_contigs();
    check_all();
    LOG("graph", info) << ptree("end");
}

pair< Read_Entry_CBPtr, Read_Entry_CBPtr >
Graph::split_read(Read_Chunk_CBPtr rc_cbptr)
{
    Read_Entry_CBPtr re_cbptr = rc_cbptr->re_bptr();
    re_cont().erase(re_cbptr);
    auto p = Read_Entry::split(rc_cbptr);
    ASSERT(not p.first->chunk_cont().empty());
    ASSERT(not p.second->chunk_cont().empty());
    ASSERT(&*p.first->chunk_cont().rbegin() == rc_cbptr);
    if (trim_tuc_step())
    {
        if (rc_cbptr->ce_bptr()->is_unmappable())
        {
            unmap_chunk(rc_cbptr.unconst());
        }
        Read_Chunk_CBPtr next_rc_cbptr = &*p.second->chunk_cont().begin();
        ASSERT(next_rc_cbptr);
        if (next_rc_cbptr->ce_bptr()->is_unmappable())
        {
            unmap_chunk(next_rc_cbptr.unconst());
        }
    }
    re_cont().insert(p.first);
    re_cont().insert(p.second);
    check(set< Read_Entry_CBPtr >{ p.first, p.second });
    return p;
} // Graph::split_read

pair< Read_Entry_CBPtr, Read_Entry_CBPtr >
Graph::split_read(Read_Entry_CBPtr re_cbptr, const Allele_Anchor& l_anchor)
{
    LOG("Graph", debug) << ptree("begin")
        .put("re_bptr", re_cbptr.to_int())
        .put("re_name", re_cbptr->name())
        .put("l_anchor", l_anchor);
    ASSERT(re_cbptr);
    auto re_bptr = re_cbptr.unconst();
    auto ce_bptr = l_anchor.ce_cbptr().unconst();
    ASSERT(not (l_anchor.is_endpoint() and l_anchor.c_right()));
    auto r_anchor = l_anchor.get_sibling(false);
    ASSERT(l_anchor.ce_cbptr() == r_anchor.ce_cbptr());
    auto l_out_pos = l_anchor.is_endpoint()? 0 : l_anchor.mut_cbptr()->rf_start();
    (void)l_out_pos;
    auto l_in_pos = l_anchor.is_endpoint()? 0 : l_anchor.mut_cbptr()->rf_end();
    auto r_in_pos = r_anchor.is_endpoint()? ce_bptr->len() : r_anchor.mut_cbptr()->rf_start();
    auto r_out_pos = r_anchor.is_endpoint()? ce_bptr->len() : r_anchor.mut_cbptr()->rf_end();
    (void)r_out_pos;
    ASSERT(l_in_pos < r_in_pos);
    //auto match_len = r_in_pos - l_in_pos;
    auto rc_bptr = ce_bptr->chunk_cont().search_read(re_bptr).unconst();
    ASSERT(rc_bptr);
    ASSERT(rc_bptr->get_c_start() <= l_out_pos);
    ASSERT(r_out_pos <= rc_bptr->get_c_end());
    // create new Read_Entry object that will hold the tail
    re_cont().erase(re_bptr);
    auto re_p = Read_Entry::split(rc_bptr, l_in_pos, r_in_pos);
    re_cont().insert(re_p.first);
    re_cont().insert(re_p.second);
    check(set< Read_Entry_CBPtr >{ re_p.first, re_p.second });
    LOG("Graph", debug) << ptree("end")
        .put("head_re", re_p.first.to_int())
        .put("tail_re", re_p.second.to_int());
    return re_p;
} // Graph::split_read

void Graph::clear_and_dispose()
{
    ce_cont().clear_and_dispose();
    re_cont().clear_and_dispose();
}

void Graph::check_all() const
{
#ifndef DISABLE_ASSERTS
    LOG("graph", info) << ptree("begin");
    size_t chunks_count_1 = 0;
    size_t chunks_count_2 = 0;
    // check read entry objects
    for (auto re_cbptr : re_cont() | referenced)
    {
        re_cbptr->check();
        chunks_count_1 += re_cbptr->chunk_cont().size();
    }
    // check contig entry objects
    for (auto ce_cbptr : ce_cont() | referenced)
    {
        ce_cbptr->check();
        chunks_count_2 += ce_cbptr->chunk_cont().size();
    }
    ASSERT(chunks_count_1 == chunks_count_2);
    LOG("graph", info) << ptree("end");
#endif
}

void Graph::check(const set< Read_Entry_CBPtr >& re_set, const set< Contig_Entry_CBPtr >& ce_set) const
{
    static_cast< void >(re_set);
    static_cast< void >(ce_set);
#ifndef DISABLE_ASSERTS
    // compute contig entries referenced by the selected read entries
    set< Contig_Entry_CBPtr > ce_extra_set;
    for (const auto& re_cbptr : re_set)
    {
        for (auto rc_cbptr : re_cbptr->chunk_cont() | referenced)
        {
            ce_extra_set.insert(rc_cbptr->ce_bptr());
        }
    }

    // compute read entries referenced by selected contig entries
    set< Read_Entry_CBPtr > re_extra_set;
    for (const auto& ce_cbptr : ce_set)
    {
        for (auto rc_cbptr : ce_cbptr->chunk_cont() | referenced)
        {
            re_extra_set.insert(rc_cbptr->re_bptr());
        }
    }

    // check read entry objects
    for (const auto& re_cbptr : re_set)
    {
        re_cbptr->check();
    }
    for (const auto& re_cbptr : re_extra_set)
    {
        if (re_set.count(re_cbptr) > 0)
        {
            continue;
        }
        re_cbptr->check();
    }

    // check contig entry objects
    for (const auto& ce_cbptr : ce_set)
    {
        ce_cbptr->check();
    }
    for (const auto& ce_cbptr : ce_extra_set)
    {
        if (ce_set.count(ce_cbptr) > 0)
        {
            continue;
        }
        ce_cbptr->check();
    }
#endif
}

void Graph::check_leaks() const
{
#ifndef DISABLE_ASSERTS
    // check Read_Entry factory
    //   in addition to Read_Entry objects in the graph, there is 1
    //   in the unique Read_Entry_Cont header node
    ASSERT(Read_Entry_Fact::used() == re_cont().size() + 1);
    // check Contig_Entry factory
    //   in addition to Contig_Entry objects in the graph, there is 1
    //   in the unique Contig_Entry_Cont header node
    ASSERT(Contig_Entry_Fact::used() == ce_cont().size() + 1);
    // check Read_Chunk factory
    //   in addition to Read_Chunk objects in the graph, there is 1
    //   in each Read_Chunk_Cont header node;
    //   and there is 1 Read_Chunk_Cont in each Read_Entry & each Contig_Entry object
    size_t num_chunks = 0;
    size_t num_muts = 0;
    size_t num_mca = 0;
    for (auto ce_cbptr : ce_cont() | referenced)
    {
        num_chunks += ce_cbptr->chunk_cont().nonconst_size();
        num_muts += ce_cbptr->mut_cont().nonconst_size();
        for (auto rc_cbptr : ce_cbptr->chunk_cont() | referenced)
        {
            num_mca += rc_cbptr->mut_ptr_cont().nonconst_size();
        }
    }
    ASSERT(Read_Chunk_Fact::used()
           == num_chunks + Read_Entry_Fact::used() + Contig_Entry_Fact::used());
    // check Mutation factory:
    //   in addition to Mutation objects in the graph
    //   there is 1 in
    //   - every Mutation_Cont header node: 1 Mutation_Cont in every Contig_Entry object
    ASSERT(Mutation_Fact::used()
           == num_muts + Contig_Entry_Fact::used());
    // check MCA factory:
    //   in addition to MCA objects in the graph,
    //   there is 1 MCA object in:
    //   - every Mutation_Ptr_Cont header: 1 in every Read_Chunk object
    //   - every Read_Chunk_Ptr_Cont header: 1 in every Mutation object
    ASSERT(Mutation_Chunk_Adapter_Fact::used()
           == num_mca + Read_Chunk_Fact::used() + Mutation_Fact::used());
    /*
    if (Mutation_Chunk_Adapter_Fact::used()
       != num_mca + Read_Chunk_Fact::used() + Mutation_Fact::used())
    {
        set< Mutation_Chunk_Adapter_Fact::index_type > s;

        s.insert(_re_cont.header_ptr()->chunk_cont().header_ptr()->mut_ptr_cont().get_root_node().to_int());
        for (auto re_cbptr : re_cont() | referenced)
        {
            s.insert(re_cbptr->chunk_cont().header_ptr()->mut_ptr_cont().get_root_node().to_int());
        }

        s.insert(_ce_cont.get_root_node()->mut_cont().header_ptr()->chunk_ptr_cont().get_root_node().to_int());
        s.insert(_ce_cont.get_root_node()->chunk_cont().header_ptr()->mut_ptr_cont().get_root_node().to_int());
        for (auto ce_cbptr : ce_cont() | referenced)
        {
            s.insert(ce_cbptr->mut_cont().header_ptr()->chunk_ptr_cont().get_root_node().to_int());
            s.insert(ce_cbptr->chunk_cont().header_ptr()->mut_ptr_cont().get_root_node().to_int());
            for (auto rc_cbptr : ce_cbptr->chunk_cont() | referenced)
            {
                s.insert(rc_cbptr->mut_ptr_cont().get_root_node().to_int());
                for (auto mca_cbptr : rc_cbptr->mut_ptr_cont() | referenced)
                {
                    s.insert(mca_cbptr.to_int());
                }
            }
            for (auto mut_cbptr : ce_cbptr->mut_cont() | referenced)
            {
                s.insert(mut_cbptr->chunk_ptr_cont().get_root_node().to_int());
            }
        }
        ASSERT(s.size() == num_mca + Read_Chunk_Fact::used() + Mutation_Fact::used());
        // next, add mca-s which were freed
        for (auto crt = _mca_fact._next_free_idn; crt; crt = _mca_fact.ns_wrapper_at(crt)._next_free_idn)
        {
            s.insert(crt.to_int());
        }
        for (Mutation_Chunk_Adapter_Fact::index_type i = 0; i < _mca_fact._cont.size(); ++i)
        {
            if (s.count(i) == 0)
            {
                cerr << "mca_bptr=" << i << " leaked\n";
                abort();
            }
        }
    }
    */
#endif
}

Graph::supercontig_list Graph::get_supercontigs(int unmappable_policy, size_t ignore_threshold) const
{
    supercontig_list res;
    unordered_set< Contig_Entry_CBPtr > visit_t(ce_cont().size());
    for (auto ce_cbptr : ce_cont() | referenced)
    {
        if (not ce_cbptr->is_normal()
            //or bitmask::any(ce_bptr->tag(), visit_mask))
            or visit_t.count(ce_cbptr))
        {
            continue;
        }
        // create new sc list and place current Contig Entry in it
        supercontig_type crt_sc;
        //bitmask::set(ce_bptr->tag(), visit_mask);
        visit_t.insert(ce_cbptr);
        crt_sc.emplace_back(ce_cbptr, true);
        // subfunction to find supercontig endpoint
        // return true iff a supercontig cycle is detected
        auto find_supercontig_endpoint = [&] (Contig_Entry_CBPtr ce_cbptr, bool c_right) -> bool
        {
            Contig_Entry_CBPtr crt_ce_cbptr = ce_cbptr;
            bool crt_c_right = c_right;
            while (true)
            {
                auto oc = crt_ce_cbptr->out_chunks_dir(crt_c_right, unmappable_policy, ignore_threshold);
                if (oc.size() != 1)
                {
                    // supercontig ends here; out-degree >1 for this ce
                    //bitmask::set(crt_ce_bptr->tag(), Contig_Entry::supercontig_endpoint_mask(crt_c_right));
                    return false;
                }
                Contig_Entry_CBPtr next_ce_cbptr = oc.begin()->first.first;
                bool same_orientation = oc.begin()->first.second;
                bool next_c_right = same_orientation? crt_c_right : not crt_c_right;
                auto back_oc = next_ce_cbptr->out_chunks_dir(not next_c_right, unmappable_policy, ignore_threshold);
                if (back_oc.size() != 1)
                {
                    // supercontig ends here; out-degree >1 for next ce
                    //bitmask::set(crt_ce_bptr->tag(), Contig_Entry::supercontig_endpoint_mask(crt_c_right));
                    return false;
                }
                ASSERT(back_oc.begin()->first == make_pair(crt_ce_cbptr, same_orientation));
                if (next_ce_cbptr == ce_cbptr)
                {
                    // supercontig is a cycle!
                    LOG("graph", info) << ptree("cycle").put("ce_ptr", ce_cbptr.to_int());
                    // must be detected during the first call, with c_right == true
                    ASSERT(c_right);
                    //bitmask::set(ce_bptr->tag(), Contig_Entry::supercontig_endpoint_mask(not c_right));
                    //bitmask::set(crt_ce_bptr->tag(), Contig_Entry::supercontig_endpoint_mask(crt_c_right));
                    return true;
                }
                crt_ce_cbptr = next_ce_cbptr;
                crt_c_right = next_c_right;
                //bitmask::set(crt_ce_bptr->tag(), visit_mask);
                visit_t.insert(crt_ce_cbptr);
                crt_sc.insert(c_right? crt_sc.end() : crt_sc.begin(), make_pair(next_ce_cbptr, next_c_right == c_right));
            }
        };
        // traverse past right endpoint, looking for supercontig endpoint
        if (not find_supercontig_endpoint(ce_cbptr, true))
        {
            // lambda returns false iff no cycle is detected, and we need to search for the other sc endpoint
            find_supercontig_endpoint(ce_cbptr, false);
        }
        res.emplace_back(move(crt_sc));
    }
    return res;
}

Size_Type Graph::skip_supercontig_bulges(supercontig_list & l) const
{
    Size_Type res = 0;
    bool done;
    do
    {
        done = true;
        for (auto it = l.begin(); it != l.end(); ++it)
        {
            pair< pair< Contig_Entry_CBPtr, bool >, pair< Contig_Entry_CBPtr, bool > > bulge_ends(*(it->begin()), *(it->rbegin()));
            // check out-degree is 1 in both directions
            auto oc_left = bulge_ends.first.first->out_chunks_dir(not bulge_ends.first.second, 3, 1);
            if (oc_left.size() != 1) continue;
            auto oc_right = bulge_ends.second.first->out_chunks_dir(bulge_ends.second.second, 3, 1);
            if (oc_right.size() != 1) continue;
            Contig_Entry_CBPtr ce_left_cbptr;
            bool same_orientation_left;
            tie(ce_left_cbptr, same_orientation_left) = oc_left.begin()->first;
            Contig_Entry_CBPtr ce_right_cbptr;
            bool same_orientation_right;
            tie(ce_right_cbptr, same_orientation_right) = oc_right.begin()->first;
            // check out-degree is 2 just outside of the bulge
            auto oc_from_left = ce_left_cbptr->out_chunks_dir(bulge_ends.first.second == same_orientation_left, 3, 1);
            if (oc_from_left.size() != 2) continue;
            auto oc_from_right = ce_right_cbptr->out_chunks_dir(not bulge_ends.second.second == same_orientation_right, 3, 1);
            if (oc_from_right.size() != 2) continue;
            // check left&right of the bulge are connected
            if (oc_from_left.count(
                    make_pair(ce_right_cbptr,
                              (bulge_ends.first.second == same_orientation_left)
                              == (bulge_ends.second.second == same_orientation_right))) == 0)
            {
                continue;
            }
            // found a bulge
            // compute size
            done = false;
            Size_Type sc_len = 0;
            for (auto p : *it)
            {
                sc_len += p.first->len();
            }
            res += sc_len;
            LOG("graph", info) << ptree("bulge")
                .put("left_ce_ptr", bulge_ends.first.first.to_int())
                .put("left_c_right", bulge_ends.first.second)
                .put("right_ce_ptr", bulge_ends.second.first.to_int())
                .put("right_c_right", bulge_ends.second.second)
                .put("len", sc_len);
            // remove bulge
            l.erase(it);
            // find supercontig ending at ce_left_cbptr
            supercontig_list::iterator it_l;
            bool it_l_orientation = false;
            for (it_l = l.begin(); it_l != l.end(); ++it_l)
            {
                if ((it_l->begin()->first == ce_left_cbptr
                     and it_l->begin()->second == (not bulge_ends.first.second == same_orientation_left)))
                {
                    it_l_orientation = false;
                    break;
                }
                if ((it_l->rbegin()->first == ce_left_cbptr
                     and it_l->rbegin()->second == (bulge_ends.first.second == same_orientation_left)))
                {
                    it_l_orientation = true;
                    break;
                }
            }
            ASSERT(it_l != l.end());
            // find supercontig ending at ce_right_cbptr
            supercontig_list::iterator it_r;
            bool it_r_orientation = false;
            for (it_r = l.begin(); it_r != l.end(); ++it_r)
            {
                if ((it_r->begin()->first == ce_right_cbptr
                     and it_r->begin()->second == (bulge_ends.second.second == same_orientation_right)))
                {
                    it_r_orientation = true;
                    break;
                }
                if ((it_r->rbegin()->first == ce_right_cbptr
                     and it_r->rbegin()->second == (not bulge_ends.second.second == same_orientation_right)))
                {
                    it_r_orientation = false;
                    break;
                }
            }
            ASSERT(it_r != l.end());
            // if the orientations are wrong, reverse it_r
            if (it_l_orientation != it_r_orientation)
            {
                supercontig_type sc;
                for (const auto& p : *it_r)
                {
                    sc.push_front(make_pair(p.first, not p.second));
                }
                swap(*it_r, sc);
                it_r_orientation = not it_r_orientation;
            }
            // merge l&r supercontigs
            ASSERT(it_l_orientation == it_r_orientation);
            if (it_l_orientation)
            {
                it_l->splice(it_l->end(), *it_r);
                l.erase(it_r);
            }
            else
            {
                it_r->splice(it_r->end(), *it_l);
                l.erase(it_l);
            }
            break;
        }
    }
    while (not done);
    return res;
}

void Graph::print_basic_stats(ostream& os) const
{
    LOG("graph", info) << ptree("begin");
    unsigned long num_re = 0;
    unsigned long bp_re = 0;
    unsigned long num_ce = 0;
    unsigned long num_ce_normal = 0;
    unsigned long bp_ce_normal = 0;
    unsigned long mut_ce_normal = 0;
    unsigned long snp_ce_normal = 0;
    unsigned long sc_branch = 0;
    for (const auto re_cbptr : re_cont() | referenced)
    {
        ++num_re;
        bp_re += re_cbptr->len();
    }
    for (const auto ce_cbptr : ce_cont() | referenced)
    {
        ++num_ce;
        if (not ce_cbptr->is_normal()) continue;
        ++num_ce_normal;
        bp_ce_normal += ce_cbptr->len();
        for (const auto mut_cbptr : ce_cbptr->mut_cont() | referenced)
        {
            ++mut_ce_normal;
            if (not mut_cbptr->is_snp()) continue;
            ++snp_ce_normal;
        }
        auto cont = ce_cbptr->out_chunks_dir(false, 3, 1);
        sc_branch += cont.size() > 1;
        cont = ce_cbptr->out_chunks_dir(true, 3, 1);
        sc_branch += cont.size() > 1;
    }

    os << "TOT\tnum_re\tbp_re\tnum_ce\tnum_ce_n\tbp_ce_n\tmut_ce_n\tsnp_ce_n\tsc_branch" << endl
       << "TOT\t"
       << num_re << "\t"
       << bp_re << "\t"
       << num_ce << "\t"
       << num_ce_normal << "\t"
       << bp_ce_normal << "\t"
       << mut_ce_normal << "\t"
       << snp_ce_normal << "\t"
       << sc_branch << endl;
}

void Graph::print_detailed_counts(ostream& os) const
{
    print_basic_stats(os);
    LOG("graph", info) << ptree("begin");
    // First read stats
    os << "RE\tname\tlen\tnum.chunks\tchunk.lens\tcontigs\n";
    for (const auto re_cbptr : re_cont() | referenced)
    {
        os << "RE\t"
           << re_cbptr->name() << '\t'
           << re_cbptr->end() - re_cbptr->start() << '\t'
           << re_cbptr->chunk_cont().size() << '\t';
        for (const auto rc_cbptr : re_cbptr->chunk_cont() | referenced)
        {
            if (rc_cbptr != &*re_cbptr->chunk_cont().begin())
            {
                os << ',';
            }
            os << (rc_cbptr->ce_bptr()->is_unmappable()? "*" : "") << rc_cbptr->get_r_len();
        }
        os << '\t';
        for (const auto rc_cbptr : re_cbptr->chunk_cont() | referenced)
        {
            if (rc_cbptr != &*re_cbptr->chunk_cont().begin())
            {
                os << ',';
            }
            os << rc_cbptr->ce_bptr().to_int();
        }
        os << '\n';
    }
    // next, contig stats
    os << "CE\tid\tlen\tunmappable\tlowcomplex"
       << "\tnum.chunks\tnum.muts\tbp.chunks\tnum.muts.chunks"
       << "\tnum.snp\tnum.ins\tnum.del\tnum.mnp\tbp.muts"
       << "\tce.deg.left\tre.deg.left\tedges.left"
       << "\tce.deg.right\tre.deg.right\tedges.right"
       << "\tce.deg.left.skip\tre.deg.left.skip\tedges.left.skip"
       << "\tce.deg.right.skip\tre.deg.right.skip\tedges.right.skip"
       << "\tallele.support\n";
    for (const auto ce_cbptr : ce_cont() | referenced)
    {
        os << "CE\t"
           << ce_cbptr.to_int() << '\t'
           << ce_cbptr->len() << '\t'
           << static_cast< int >(ce_cbptr->is_unmappable()) << '\t'
           << static_cast< int >(ce_cbptr->is_lowcomplex()) << '\t'
           << ce_cbptr->chunk_cont().size() << '\t'
           << ce_cbptr->mut_cont().size() << '\t';
        size_t num_bp_chunks = 0;
        size_t num_muts_chunks = 0;
        for (const auto rc_cbptr : ce_cbptr->chunk_cont() | referenced)
        {
            num_bp_chunks += rc_cbptr->get_r_len();
            num_muts_chunks += rc_cbptr->mut_ptr_cont().size();
        }
        os << num_bp_chunks << '\t'
           << num_muts_chunks << '\t';
        size_t n_snp = 0;
        size_t n_ins = 0;
        size_t n_del = 0;
        size_t n_mnp = 0;
        size_t total_mut_bp = 0;
        for (const auto mut_cbptr : ce_cbptr->mut_cont() | referenced)
        {
            if (mut_cbptr->is_snp())
            {
                ++n_snp;
            }
            else if (mut_cbptr->is_ins())
            {
                ++n_ins;
            }
            else if (mut_cbptr->is_del())
            {
                ++n_del;
            }
            else
            {
                ++n_mnp;
            }
            total_mut_bp += mut_cbptr->rf_len() + mut_cbptr->seq_len();
        }
        os << n_snp << '\t'
           << n_ins << '\t'
           << n_del << '\t'
           << n_mnp << '\t'
           << total_mut_bp << '\t';

        auto print_neighbours_cont_stats = [] (
            std::ostream& os,
            Contig_Entry_CBPtr ce_cbptr,
            bool c_right,
            bool print_unmap_range,
            const Contig_Entry::out_chunks_dir_type& cont)
        {
            // contig out-degree
            os << cont.size() << '\t';
            // read out-degree
            size_t s = 0;
            for (const auto& p : cont)
            {
                s += p.second.size();
            }
            os << s << '\t';
            // breakdown of inter-contig edges
            if (cont.size() == 0)
            {
                os << '.';
            }
            else
            {
                bool first = true;
                for (const auto& p : cont)
                {
                    if (not first)
                    {
                        os << ',';
                    }
                    first = false;
                    Contig_Entry_CBPtr tmp_ce_cbptr;
                    bool tmp_flip;
                    tie(tmp_ce_cbptr, tmp_flip) = p.first;
                    os << '(' << tmp_ce_cbptr.to_int()
                       << ',' << int(tmp_flip)
                       << ',' << p.second.size();
                    if (print_unmap_range)
                    {
                        Size_Type min_skipped_len;
                        Size_Type max_skipped_len;
                        tie(min_skipped_len, max_skipped_len) = ce_cbptr->unmappable_neighbour_range(c_right, p.second);
                        os << ',' << min_skipped_len
                           << ',' << max_skipped_len;
                    }
                    os << ')';
                }
            }
        }; // print_neighbours_cont_stats

        {
            print_neighbours_cont_stats(os, ce_cbptr, false, false, ce_cbptr->out_chunks_dir(false, 1));
            os << '\t';
            print_neighbours_cont_stats(os, ce_cbptr, true, false, ce_cbptr->out_chunks_dir(true, 1));
            os << '\t';
        }
        if (ce_cbptr->is_unmappable())
        {
            os << ".\t.\t.\t.\t.\t.\t.";
        }
        else
        {
            print_neighbours_cont_stats(os, ce_cbptr, false, true, ce_cbptr->out_chunks_dir(false, 3, 1));
            os << '\t';
            print_neighbours_cont_stats(os, ce_cbptr, true, true, ce_cbptr->out_chunks_dir(true, 3, 1));
            auto ac = Allele_Anchor::connect(
                Allele_Anchor(ce_cbptr, false).read_support(2),
                Allele_Anchor(ce_cbptr, true).read_support(2));
            os << '\t';
            bool first = true;
            for (const auto& p : ac)
            {
                if (not first)
                {
                    os << ",";
                }
                first = false;
                auto l_allele = p.first.first;
                auto r_allele = p.first.second;
                os << "(" << l_allele.ce_next_cbptr().to_int() << "," << l_allele.same_orientation() << ","
                   << r_allele.ce_next_cbptr().to_int() << "," << r_allele.same_orientation() << ","
                   << p.second.size() << ")";
            }
        }
        os << '\n';
    } //for (ce_cbptr : ce_cont()
    print_supercontig_stats(os);
    print_mutations(os);
}

void Graph::print_supercontig_stats(ostream& os) const
{
    LOG("graph", info) << ptree("begin");
    auto sc_list = get_supercontigs(3, 1);
    auto print_list = [&os] (const supercontig_list & l, const string & tag)
    {
        os << tag <<"\tbp.contigs\tnum.contigs\tcontig.lens\tcontigs\tdeg.left\tdeg.right" << endl;
        for (auto& sc: l)
        {
            size_t len = 0;
            for (auto& t : sc)
            {
                len += t.first->len();
            }
            os << tag << "\t" << len << "\t" << sc.size() << "\t";
            bool first = true;
            for (auto& t : sc)
            {
                if (not first)
                {
                    os << ",";
                }
                first = false;
                os << t.first->len();
            }
            os << "\t";
            first = true;
            for (auto& t : sc)
            {
                if (not first)
                {
                    os << ",";
                }
                first = false;
                os << "(" << t.first.to_int() << "," << t.second << ")";
            }
            os << "\t";
            auto oc = sc.begin()->first->out_chunks_dir(not sc.begin()->second, 3, 1);
            os << oc.size() << "\t";
            oc = sc.rbegin()->first->out_chunks_dir(sc.rbegin()->second, 3, 1);
            os << oc.size() << endl;
        }
    };
    print_list(sc_list, "SC");
    auto bulge_len = skip_supercontig_bulges(sc_list);
    print_list(sc_list, "SSC");
    os << "BLGL\tlen" << endl
       << "BLGL\t" << bulge_len << endl;
}

void Graph::test_mutation_allele_swapping()
{
    for (auto re_bptr : re_cont() | referenced)
    {
        auto pos = re_bptr->start();
        while (pos < re_bptr->end())
        {
            auto rc_bptr = re_bptr->chunk_cont().get_chunk_with_pos(pos);
            auto ce_bptr = rc_bptr->ce_bptr();
            ASSERT(rc_bptr);
            auto rc_init_pos = not rc_bptr->get_rc()? rc_bptr->get_start_pos() : rc_bptr->get_end_pos();
            auto rc_next_pos = rc_init_pos.jump_to_brk(pos, false);
            ASSERT(rc_next_pos.r_pos == pos);
            auto match_len = rc_next_pos.get_match_len(not rc_bptr->get_rc());
            if (match_len > 0) rc_next_pos.advance(not rc_bptr->get_rc());
            ASSERT(rc_next_pos.get_match_len(not rc_bptr->get_rc()) == 0);
            ASSERT(rc_next_pos.r_pos >= pos);
            pos = rc_next_pos.r_pos;
            auto mca_cit = rc_next_pos.mca_cit;
            if (mca_cit == (not rc_bptr->get_rc()
                            ? rc_bptr->mut_ptr_cont().end()
                            : rc_bptr->mut_ptr_cont().begin()))
            {
                ASSERT(pos == rc_bptr->get_r_end());
                continue;
            }
            if (rc_bptr->get_rc()) --mca_cit;
            Contig_Entry::swap_mutation_alleles(ce_bptr, mca_cit->mut_cbptr());
            check(set< Contig_Entry_CBPtr >({ ce_bptr }));
        }
    }
    check_all();
}

ptree Graph::to_ptree() const
{
    return ptree().put("re_cont", re_cont()).put("ce_cont", ce_cont());
}

ptree Graph::factory_stats() const
{
    return ptree()
        .put("Factory<Read_Entry>", _re_fact.stats())
        .put("Factory<Contig_Entry>", _ce_fact.stats())
        .put("Factory<Read_Chunk>", _rc_fact.stats())
        .put("Factory<Mutation>", _mut_fact.stats())
        .put("Factory<Mutation_Chunk_Adapter>", _mca_fact.stats());
}

void Graph::save(ostream& os) const
{
    LOG("io", info) << ptree("begin");
    ASSERT(_re_cont.header_ptr().to_int() == 0);
    ASSERT(_ce_cont.get_root_node().to_int() == 0);
    // dump the 5 factories
    _mut_fact.save(os);
    _mca_fact.save(os);
    _rc_fact.save(os);
    _re_fact.save(os);
    _ce_fact.save(os);
    // dump re_cont and ce_cont
    os.write(reinterpret_cast< const char* >(&_re_cont), sizeof(_re_cont));
    os.write(reinterpret_cast< const char* >(&_ce_cont), sizeof(_ce_cont));
    // save relevant strings in huge a list, resetting entry contents
    // list order:
    // - for every read entry in _re_cont:
    //   - read name
    //   - start_seq
    //   - end_seq
    // - for every contig entry in _ce_cont:
    //   - base sequence
    //   - for every mutation in its mut_cont
    //     - alternate sequence
    size_t n_strings = 0;
    size_t n_bytes = 0;
    for (const auto re_cbptr : re_cont() | referenced)
    {
        re_cbptr->save_strings(os, n_strings, n_bytes);
    }
    for (const auto ce_cbptr : ce_cont() | referenced)
    {
        ce_cbptr->save_strings(os, n_strings, n_bytes);
        for (const auto mut_cbptr : ce_cbptr->mut_cont() | referenced)
        {
            mut_cbptr->save_strings(os, n_strings, n_bytes);
        }
    }
    LOG("io", info) << ptree("end").put("n_strings", n_strings).put("n_bytes", n_bytes);
}

void Graph::load(istream& is)
{
    LOG("io", info) << ptree("begin");
    ASSERT(_re_cont.header_ptr().to_int() == 0);
    ASSERT(_ce_cont.get_root_node().to_int() == 0);
    ASSERT(_re_cont.empty());
    ASSERT(_ce_cont.empty());
    // load 5 factories
    _mut_fact.load(is);
    _mca_fact.load(is);
    _rc_fact.load(is);
    _re_fact.load(is);
    _ce_fact.load(is);
    // load re_cont and ce_cont
    is.read(reinterpret_cast< char* >(&_re_cont), sizeof(_re_cont));
    is.read(reinterpret_cast< char* >(&_ce_cont), sizeof(_ce_cont));
    // construct in-place empty strings in container headers (ow, destruction causes SEGV)
    _re_cont.header_ptr()->init_strings();
    _ce_cont.get_root_node()->init_strings();
    _ce_cont.get_root_node()->mut_cont().header_ptr()->init_strings();
    // load strings
    size_t n_strings = 0;
    size_t n_bytes = 0;
    for (auto re_bptr : re_cont() | referenced)
    {
        re_bptr->init_strings();
        re_bptr->load_strings(is, n_strings, n_bytes);
    }
    for (auto ce_bptr : ce_cont() | referenced)
    {
        ce_bptr->init_strings();
        ce_bptr->mut_cont().header_ptr()->init_strings();
        ce_bptr->load_strings(is, n_strings, n_bytes);
        for (auto mut_bptr : ce_bptr->mut_cont() | referenced)
        {
            mut_bptr->init_strings();
            mut_bptr->load_strings(is, n_strings, n_bytes);
        }
    }
    LOG("io", info) << ptree("end").put("n_strings", n_strings).put("n_bytes", n_bytes);
    check_all();
}

void Graph::get_terminal_reads(ostream& os) const
{
    LOG("graph", info) << ptree("begin");
    for (auto ce_cbptr : ce_cont() | referenced)
    {
        for (int dir = 0; dir < 2; ++dir)
        {
            bool c_right = (dir == 1);
            auto tmp = ce_cbptr->out_chunks_dir(c_right, 3);
            if (tmp.empty())
            {
                // scontig ends in direction dir
                Read_Chunk_CBPtr rc_cbptr = (not c_right?
                                             &*ce_cbptr->chunk_cont().begin()
                                             : &*ce_cbptr->chunk_cont().iintersect(ce_cbptr->len(), ce_cbptr->len()).begin());
                Read_Entry_CBPtr re_cbptr = rc_cbptr->re_bptr();
                os << ">" << re_cbptr->name() << ":" << re_cbptr->start() << "-" << re_cbptr->end()
                   << " " << (c_right == rc_cbptr->get_rc()? "1" : "0") << endl
                   << re_cbptr->get_seq().revcomp(c_right == rc_cbptr->get_rc()) << endl;
            }
        }
    }
    LOG("graph", info) << ptree("end");
}

void Graph::export_gfa(ostream& os, bool show_mutations) const
{
    os << "H\tVN:Z:1.0" << endl;
    for (const auto ce_cbptr : ce_cont() | referenced)
    {
        if (not ce_cbptr->is_normal()) continue;
        if (not show_mutations)
        {
            os << "S\t" << ce_cbptr.to_int() << "\t" << ce_cbptr->seq() << endl;
        }
        else // split contigs at mutations
        {
            unsigned i = 0;
            auto mut_it = ce_cbptr->mut_cont().begin();
            auto prev_mut_it = mut_it;
            while (true)
            {
                // first print base sequence before mutation with index i
                Size_Type base_start = (i > 0? prev_mut_it->rf_end() : 0);
                Size_Type base_end = (mut_it != ce_cbptr->mut_cont().end()? mut_it->rf_start() : ce_cbptr->len());
                ASSERT(base_end > base_start);
                os << "S\t" << ce_cbptr.to_int() << ":B:" << i
                   << "\t" << ce_cbptr->substr(base_start, base_end - base_start) << endl;
                // if there was a previous mutation, links base to both alleles
                if (i > 0)
                {
                    os << "L\t"
                       << ce_cbptr.to_int() << ":M:" << i - 1 << ":0\t+\t"
                       << ce_cbptr.to_int() << ":B:" << i << "\t+\t0M" << endl;
                    os << "L\t"
                       << ce_cbptr.to_int() << ":M:" << i - 1 << ":1\t+\t"
                       << ce_cbptr.to_int() << ":B:" << i << "\t+\t0M" << endl;
                }
                if (mut_it == ce_cbptr->mut_cont().end()) break;
                // print mutation base sequence
                os << "S\t" << ce_cbptr.to_int() << ":M:" << i << ":0"
                   << "\t" << ce_cbptr->substr(mut_it->rf_start(), mut_it->rf_len()) << endl;
                os << "L\t"
                       << ce_cbptr.to_int() << ":B:" << i << "\t+\t"
                       << ce_cbptr.to_int() << ":M:" << i << ":0\t+\t0M" << endl;
                // print mutation alt sequence
                os << "S\t" << ce_cbptr.to_int() << ":M:" << i << ":1"
                   << "\t" << mut_it->seq() << endl;
                os << "L\t"
                       << ce_cbptr.to_int() << ":B:" << i << "\t+\t"
                       << ce_cbptr.to_int() << ":M:" << i << ":1\t+\t0M" << endl;
                ++i;
                prev_mut_it = mut_it;
                ++mut_it;
            }
        }
        for (int c_right = 0; c_right < 2; ++c_right)
        {
            auto oc = ce_cbptr->out_chunks_dir(c_right, 3, 1);
            for (const auto& p : oc | ba::map_keys)
            {
                if (p.first.to_int() < ce_cbptr.to_int()) continue;
                auto get_enpoint_tag = [&show_mutations] (Contig_Entry_CBPtr ce_cbptr, bool c_right, bool edge_end)
                {
                    ostringstream tmp;
                    tmp << ce_cbptr.to_int();
                    if (show_mutations)
                    {
                        tmp << ":B:" << (c_right == edge_end? 0 : ce_cbptr->mut_cont().nonconst_size());
                    }
                    tmp << "\t" << (c_right? "+" : "-");
                    return tmp.str();
                };
                os << "L\t" << get_enpoint_tag(ce_cbptr, c_right, false) << "\t"
                   << get_enpoint_tag(p.first, p.second == c_right, true) << "\t0M" << endl;
            }
        }
    }
} // Graph::export_gfa

void Graph::unmap_short_contigs(unsigned min_len, unsigned max_deg)
{
    LOG("Graph", info) << ptree("begin");
    for (bool done = false; not done; )
    {
        done = true;
        for (auto ce_bptr : ce_cont() | referenced)
        {
            if (not ce_bptr->is_normal()) continue;
            if (ce_bptr->len() >= min_len) continue;
            auto oc_left = ce_bptr->out_chunks_dir(false, 3, 1);
            auto oc_right = ce_bptr->out_chunks_dir(true, 3, 1);
            if (oc_left.size() <= max_deg and oc_right.size() <= max_deg) continue;
            // unmap contig entry
            done = false;
            /*
            vector< pair< Read_Entry_BPtr, Range_Type > > unmap_v;
            for (auto rc_bptr : ce_bptr->chunk_cont() | referenced)
            {
                unmap_v.emplace_back(rc_bptr->re_bptr(), Range_Type(rc_bptr->get_r_start(), rc_bptr->get_r_end()));
            }
            for (const auto& p : unmap_v)
            {
                unmap_re_region(p.first, p.second);
            }
            */
            Unmapper::re_set_type unmap_re_set;
            for (auto rc_bptr : ce_bptr->chunk_cont() | referenced)
            {
                unmap_re_set[rc_bptr->re_bptr()].insert(Range_Type(rc_bptr->get_r_start(), rc_bptr->get_r_end()));
            }
            unmap_re_regions(move(unmap_re_set));
            // restart loop
            break;
        }
    }
    LOG("Graph", info) << ptree("end");
} // Graph::unmap_short_contigs

void Graph::load_bwt(const string& bwt_prefix)
{
    // load BWT and SSA of the reads
    LOG("graph", info) << ptree("bwt").put("file", bwt_prefix + ".bwt");
    _index_set.pBWT = new BWT(bwt_prefix + ".bwt");
    LOG("graph", info) << ptree("ssa").put("file", bwt_prefix + ".ssa");
    _index_set.pSSA = new SampledSuffixArray(bwt_prefix + ".ssa");
    LOG("graph", info) << ptree("read_list").put("file", bwt_prefix + ".id.gz");
    {
        zstr::ifstream ifs(bwt_prefix + ".id.gz");
        string sid;
        while (getline(ifs, sid))
        {
            unsigned iid = _iid_to_sid_m.size();
            _sid_to_iid_m[sid] = iid;
            _iid_to_sid_m[iid] = sid;
        }
    }
    LOG("graph", info) << ptree("end");
} // Graph::load_bwt

void Graph::load_aux_bwt(const string& aux_bwt_file)
{
    // load BWT of Illumina reads
    LOG("graph", info) << ptree("bwt").put("file", aux_bwt_file);
    _aux_index_set.pBWT = new BWT(aux_bwt_file);
    LOG("graph", info) << ptree("end");
} // Graph::load_aux_bwt

Graph::find_reads_with_seq_type Graph::find_reads_with_seq(const Seq_Proxy_Type& seq, unsigned max_count) const
{
    find_reads_with_seq_type res;
    auto find_seq_strand = [&] (const Seq_Proxy_Type& s, bool strand) {
        auto interval = BWTAlgorithms::findInterval(_index_set, s);
        if (not interval.isValid()
            or (max_count > 0 and interval.size() > max_count))
        {
            return;
        }
        for (auto i = interval.lower; i <= interval.upper; ++i)
        {
            auto elem = _index_set.pSSA->calcSA(i, _index_set.pBWT);
            auto iid = elem.getID();
            res[strand].emplace_back(_iid_to_sid_m.at(iid), elem.getPos());
        }
    };
    find_seq_strand(seq, false);
    find_seq_strand(seq.revcomp(), true);
    return res;
} // Graph::find_reads_with_seq

Graph::find_read_entries_with_seq_type
Graph::find_read_entries_with_seq(const Seq_Proxy_Type& seq, unsigned max_count) const
{
    find_read_entries_with_seq_type res;
    res.second = false;
    auto r = find_reads_with_seq(seq, max_count);
    for (int st = 0; st < 2; ++st)
    {
        for (const auto& p : r[st])
        {
            auto re_cbptr = re_cont().find(p.first);
            if (not re_cbptr)
            {
                res.second = true;
                continue;
            }
            res.first[st].emplace_back(re_cbptr, p.second);
        }
    }
    return res;
} // Graph::find_read_entries_with_seq

void Graph::compute_mutation_uniqueness(Size_Type flank_len)
{
    LOG("graph", info) << ptree("begin");
    if (not (_index_set.pBWT and _index_set.pSSA))
    {
        LOG("graph", error) << "compute_mutation_uniqueness: BWT index is required" << endl;
        exit(EXIT_FAILURE);
    }
    for (auto ce_bptr : ce_cont() | referenced)
    {
        if (not ce_bptr->is_normal()) continue;
        if (not ce_bptr->separated_mutations(flank_len, true))
        {
            static bool warned = false;
            if (not warned)
            {
                warned = true;
                LOG("graph", warning) << ptree("non_separated_mutations")
                    .put("flank_len", flank_len)
                    .put("ce_ptr", ce_bptr.to_int());
            }
            continue;
        }
        for (auto mut_bptr : ce_bptr->mut_cont() | referenced)
        {
            Seq_Type flank_left = ce_bptr->substr(mut_bptr->rf_start() - flank_len, flank_len);
            Seq_Type flank_right = ce_bptr->substr(mut_bptr->rf_end(), flank_len);
            for (int al = 0; al < 2; ++al)
            {
                Seq_Type s = flank_left;
                s += (al == 0
                      ? Seq_Type(ce_bptr->substr(mut_bptr->rf_start(), mut_bptr->rf_len()))
                      : mut_bptr->seq());
                s += flank_right;
                auto r = find_read_entries_with_seq(s, 1000);
                if (r.second)
                {
                    // hits not in graph
                    mut_bptr->uniq(al) = 3;
                    continue;
                }
                if (r.first.empty())
                {
                    // too many hits
                    mut_bptr->uniq(al) = 4;
                    continue;
                }
                bool unique = true;
                for (int st = 0; st < 2; ++st)
                {
                    for (const auto& p : r.first[st])
                    {
                        auto rc_cbptr = p.first->chunk_cont().get_chunk_with_pos(p.second + flank_len);
                        if (not rc_cbptr
                            or rc_cbptr->ce_bptr() != ce_bptr
                            or rc_cbptr->get_c_start() >= mut_bptr->rf_end() + flank_len
                            or rc_cbptr->get_c_end() + flank_len <= mut_bptr->rf_start())
                        {
                            unique = false;
                            break;
                        }
                    }
                    if (not unique) break;
                }
                mut_bptr->uniq(al) = (unique? 1 : 2);
            } // for al
        } // for mut_bptr
    } // for ce_bptr
    LOG("graph", info) << ptree("begin");
} // Graph::compute_mutation_uniqueness

void Graph::compute_mutation_copy_num(Size_Type flank_len)
{
    LOG("graph", info) << ptree("begin");
    if (not aux_index_set().pBWT)
    {
        LOG("graph", error) << "compute_mutation_copy_num: auxiliary BWT index is required" << endl;
        exit(EXIT_FAILURE);
    }
    if (get_aux_coverage() < 0)
    {
        compute_aux_coverage(flank_len);
    }
    for (auto ce_bptr : ce_cont() | referenced)
    {
        if (not ce_bptr->is_normal()) continue;
        if (not ce_bptr->separated_mutations(flank_len, true))
        {
            static bool warned = false;
            if (not warned)
            {
                warned = true;
                LOG("graph", warning) << ptree("non_separated_mutations")
                    .put("flank_len", flank_len)
                    .put("ce_ptr", ce_bptr.to_int());
            }
            continue;
        }
        for (auto mut_bptr : ce_bptr->mut_cont() | referenced)
        {
            if (mut_bptr->uniq(0) != 1 or mut_bptr->uniq(1) != 1)
            {
                LOG("graph", debug) << ptree("non_unique_mutation")
                    .put("mut_bptr", mut_bptr.to_int())
                    .put("mut_ptr", mut_bptr.raw());
                continue;
            }
            Seq_Type flank_left = ce_bptr->substr(mut_bptr->rf_start() - flank_len, flank_len);
            Seq_Type flank_right = ce_bptr->substr(mut_bptr->rf_end(), flank_len);
            for (int al = 0; al < 2; ++al)
            {
                Seq_Type s = flank_left;
                s += (al == 0
                      ? Seq_Type(ce_bptr->substr(mut_bptr->rf_start(), mut_bptr->rf_len()))
                      : mut_bptr->seq());
                s += flank_right;

                auto cnt = BWTAlgorithms::countSequenceOccurrences(s, aux_index_set());
                int cn = round((double)cnt / get_aux_coverage());
                if (0 <= cn and cn <= 2)
                {
                    mut_bptr->copy_num(al) = cn;
                }
                else if (cn > 2)
                {
                    mut_bptr->copy_num(al) = -2;
                }
            } // for al
        } // for mut_bptr
    } // for ce_bptr
    LOG("graph", info) << ptree("end");
} // Graph::compute_mutation_copy_num

void Graph::compute_aux_coverage(Size_Type flank_len)
{
    LOG("graph", info) << ptree("begin");
    if (not aux_index_set().pBWT)
    {
        LOG("graph", error) << "compute_aux_coverage: auxiliary BWT index is required" << endl;
        exit(EXIT_FAILURE);
    }
    multiset< size_t > cov_s;
    for (auto ce_bptr : ce_cont() | referenced)
    {
        if (not ce_bptr->is_normal()) continue;
        if (not ce_bptr->separated_mutations(flank_len, true))
        {
            static bool warned = false;
            if (not warned)
            {
                warned = true;
                LOG("graph", warning) << ptree("non_separated_mutations")
                    .put("flank_len", flank_len)
                    .put("ce_ptr", ce_bptr.to_int());
            }
            continue;
        }
        for (auto mut_bptr : ce_bptr->mut_cont() | referenced)
        {
            if (mut_bptr->uniq(0) != 1 or mut_bptr->uniq(1) != 1)
            {
                LOG("graph", debug) << ptree("non_unique_mutation")
                    .put("mut_bptr", mut_bptr.to_int())
                    .put("mut_ptr", mut_bptr.raw());
                continue;
            }
            Seq_Type flank_left = ce_bptr->substr(mut_bptr->rf_start() - flank_len, flank_len);
            Seq_Type flank_right = ce_bptr->substr(mut_bptr->rf_end(), flank_len);
            for (int al = 0; al < 2; ++al)
            {
                Seq_Type s = flank_left;
                s += (al == 0
                      ? Seq_Type(ce_bptr->substr(mut_bptr->rf_start(), mut_bptr->rf_len()))
                      : mut_bptr->seq());
                s += flank_right;
                auto cnt = BWTAlgorithms::countSequenceOccurrences(s, aux_index_set());
                LOG("graph", debug) << ptree("allele_coverage")
                    .put("mut_bptr", mut_bptr.to_int())
                    .put("mut_ptr", mut_bptr.raw())
                    .put("al", al)
                    .put("cnt", cnt);
                cov_s.insert(cnt);
            } // for al
        } // for mut_bptr
    } // for ce_bptr
    // ignore 10% lowest and highest values
    size_t cov_cnt = 0;
    size_t cov_sum = 0;
    size_t i = 0;
    for (auto v : cov_s)
    {
        if (i >= cov_s.size() / 10 and i < cov_s.size() * 9 / 10)
        {
            cov_sum += v;
            ++cov_cnt;
        }
        ++i;
    }
    _aux_coverage = cov_cnt > 0? cov_sum / cov_cnt : 0;
    LOG("graph", info) << ptree("end").put("res", _aux_coverage);
} // Graph::compute_aux_coverage

} // namespace MAC
