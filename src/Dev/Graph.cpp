#include "Graph.hpp"

#include "globals.hpp"
#include "Cigar.hpp"
#include "indent.hpp"
#include "print_seq.hpp"
#include "../Util/Util.h"

using namespace std;


namespace MAC
{

void Graph::add_read(string&& name, Seq_Type&& seq)
{
    // create read entry and place it in container
    Read_Entry_BPtr re_bptr = Read_Entry_Fact::new_elem(std::move(name), seq.size());
    //cerr << indent::tab << "re:\n" << indent::inc << re << indent::dec;
    re_cont().insert(re_bptr);

    // create contig entry and place it in container
    Contig_Entry_BPtr ce_bptr = Contig_Entry_Fact::new_elem(std::move(seq));
    //cerr << indent::tab << "ce:\n" << indent::inc << ce << indent::dec;
    ce_cont().insert(ce_bptr);

    // create initial read chunk
    Read_Chunk_BPtr rc_bptr = Read_Chunk_Fact::new_elem(re_bptr, ce_bptr);
    // add it to RE and CE containers
    re_bptr->chunk_cont().insert(rc_bptr);
    ce_bptr->chunk_cont().insert(rc_bptr);

    ASSERT(check(set< Read_Entry_CBPtr >( { re_bptr })));
}

bool Graph::cut_contig_entry(Contig_Entry_BPtr ce_bptr, Size_Type c_brk, Mutation_CBPtr mut_left_cbptr)
{
    // there is nothing to cut if:
    if (// cut is at the start, and no insertion goes to the left
        (c_brk == 0 and mut_left_cbptr == nullptr)
        // cut is at the end, and:
        or (c_brk == ce_bptr->get_len()
            and (// there are no mutations
                ce_bptr->mut_cont().size() == 0
                // or the last one is not an insertion
                or not ce_bptr->mut_cont().rbegin()->is_ins()
                // or the last insertion is not at the end
                or ce_bptr->mut_cont().rbegin()->get_start() < ce_bptr->get_len())))
    {
        return false;
    }

    //cerr << "before cutting contig entry: " << (void*)ce_bptr << " at " << c_brk << " mut_left_cbptr=" << (void*)mut_left_cbptr << '\n' << *this;

    // create new contig entry object; set base sequences; add new one to container
    Contig_Entry_BPtr ce_new_bptr = Contig_Entry_Fact::new_elem(string(ce_bptr->seq().substr(c_brk)));
    ce_bptr->seq().resize(c_brk);
    ce_cont().insert(ce_new_bptr);

    // split any mutations that span c_pos
    {
        Mutation_BPtr mut_bptr;
        while ((mut_bptr = ce_bptr->mut_cont().find_span_pos(c_brk).unconst()) != nullptr)
        {
            ASSERT(mut_bptr->get_start() < c_brk and c_brk < mut_bptr->get_end());
            ASSERT(mut_bptr != mut_left_cbptr);
            ce_bptr->cut_mutation(mut_bptr, c_brk - mut_bptr->get_start(), 0);
        }
    }

    // split Mutation_Cont, save rhs in new Contig_Entry
    ce_new_bptr->mut_cont() = ce_bptr->mut_cont().split(c_brk, mut_left_cbptr);

    // unlink Read_Chunk objects from their RE containers
    ce_bptr->chunk_cont().erase_from_re_cont();

    // split Read_Chunk_Cont, save rhs in new Contig_Entry
    ce_new_bptr->chunk_cont() = ce_bptr->chunk_cont().split(c_brk, mut_left_cbptr);
    ASSERT(ce_bptr->chunk_cont().max_end() <= c_brk);
    ASSERT(c_brk <= ce_new_bptr->chunk_cont().begin()->get_c_start());

    // rebase all mutations and read chunks from the rhs to the breakpoint
    ce_new_bptr->mut_cont().shift(-int(c_brk));
    ce_new_bptr->chunk_cont().shift(-int(c_brk));

    // link back the chunks into their RE containers
    ce_bptr->chunk_cont().insert_into_re_cont();
    ce_new_bptr->chunk_cont().insert_into_re_cont();

    // remove unused Mutation objects
    ce_bptr->mut_cont().drop_unused();
    ce_new_bptr->mut_cont().drop_unused();

    // if either contig has no read chunks mapped to it, remove it
    if (ce_bptr->chunk_cont().size() == 0)
    {
        ce_cont().erase(ce_bptr);
        Contig_Entry_Fact::del_elem(ce_bptr);
    }
    if (ce_new_bptr->chunk_cont().size() == 0)
    {
        ce_cont().erase(ce_new_bptr);
        Contig_Entry_Fact::del_elem(ce_new_bptr);
    }

    ASSERT(check(set< Contig_Entry_CBPtr >( { ce_bptr, ce_new_bptr })));
    return true;
} // cut_contig_entry

bool Graph::cut_read_chunk(Read_Chunk_BPtr rc_bptr, Size_Type r_brk, bool force)
{
    ASSERT(rc_bptr->get_r_len() > 0);
    ASSERT(rc_bptr->get_r_start() <= r_brk and r_brk <= rc_bptr->get_r_end());
    if (rc_bptr->is_unmappable())
    {
        // never cut unmappable chunks
        return false;
    }
    if (r_brk == rc_bptr->get_r_start() or r_brk == rc_bptr->get_r_end())
    {
        // cut is at the edge; it must be forced
        ASSERT(force);
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
            if (rc_bptr->get_c_end() != rc_bptr->ce_bptr()->get_len())
            {
                // if chunk ends in insertion, keep it on lhs, along with the rest of the chunk
                Mutation_CBPtr mut_left_cbptr = nullptr;
                if (rc_bptr->mut_ptr_cont().size() > 0
                    and rc_bptr->mut_ptr_cont().rbegin()->mut_cbptr()->get_start() == rc_bptr->get_c_end())
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
        ASSERT(not force);
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
                                  min(pos.mut_offset, mut_bptr->get_len()),
                                  min(pos.mut_offset, mut_bptr->get_seq_len()));
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
            and pos.prev_mut().get_start() == pos.c_pos)
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

bool Graph::cut_read_entry(Read_Entry_BPtr re_bptr, Size_Type r_brk, bool force)
{
    if (r_brk == 0 or r_brk == re_bptr->get_len())
    {
        // cut is on the edge of the read; it must be forced
        ASSERT(force);
        Read_Chunk_BPtr rc_bptr = (r_brk == 0? &*re_bptr->chunk_cont().begin() : &*re_bptr->chunk_cont().rbegin());
        return cut_read_chunk(rc_bptr, r_brk, true);
    }
    else
    {
        // cut is inside the read, it must not be forced
        ASSERT(not force);
        Read_Chunk_BPtr rc_bptr = re_bptr->chunk_cont().get_chunk_with_pos(r_brk).unconst();
        // the chunk to be cut must exist
        ASSERT(rc_bptr->get_r_start() <= r_brk and r_brk < rc_bptr->get_r_end());
        return cut_read_chunk(rc_bptr, r_brk, false);
    }
}

/*
void Graph::remap_chunks(map< Read_Chunk_CPtr, shared_ptr< Read_Chunk > >& rc_map,
                         Mutation_Cont& extra_mut_cont)
{
    if (rc_map.size() == 0)
        return;
    const Contig_Entry* ce1_cptr = rc_map.begin()->second->get_ce_ptr();
    const Contig_Entry* ce2_cptr = rc_map.begin()->first->get_ce_ptr();

    // first construct map with extra Mutation objects as keys
    map< Mutation_CPtr, Mutation_CPtr > extra_mut_map;
    for (auto it = extra_mut_cont.begin(); it != extra_mut_cont.end(); ++it)
    {
        extra_mut_map[&*it] = NULL;
    }

    // next, go through new read chunks
    // each time an extra mutation is first used, add it to c1 and save translation in map
    for (auto it2 = rc_map.begin(); it2 != rc_map.end(); ++it2)
    {
        Read_Chunk_CPtr c2rc_cptr = it2->first;
        Read_Chunk* c1rc_ptr = it2->second.get();
        ASSERT(c2rc_cptr->get_ce_ptr() == ce2_cptr);
        ASSERT(c1rc_ptr->get_ce_ptr() == ce1_cptr);
        ASSERT(c2rc_cptr->get_re_ptr() == c1rc_ptr->get_re_ptr());
        for (size_t i = 0; i < c1rc_ptr->get_mut_ptr_cont().size(); ++i)
        {
            Mutation_CPtr extra_mut_cptr = c1rc_ptr->get_mut_ptr_cont()[i];
            if (extra_mut_map.count(extra_mut_cptr) == 0)
            {
                // this mutation is not from the extra Mutations container, so it already exists in c1
                continue;
            }
            if (extra_mut_map[extra_mut_cptr] == NULL)
            {
                // first time we see this mutation
                // we add it to ce1, and save the pointer
                modify_contig_entry(ce1_cptr, [&] (Contig_Entry& ce) {
                    extra_mut_map[extra_mut_cptr] = add_mut_to_cont(ce.mut_cont(), *extra_mut_cptr);
                });
            }
            ASSERT(extra_mut_map[extra_mut_cptr] != NULL);
            c1rc_ptr->mut_ptr_cont()[i] = extra_mut_map[extra_mut_cptr];
        }
        // done moving mutations
        // now replace the read chunk in its read entry, and update pointers
        modify_read_chunk(c2rc_cptr, [&] (Read_Chunk& rc) {
            rc = *c1rc_ptr;
        });
        // now c2rc_cptr has the c1rc object
        // add chunk to ce1
        modify_contig_entry(ce1_cptr, [&] (Contig_Entry& ce) {
            ce.add_chunk(c2rc_cptr);
        });
    }
}
*/

void Graph::merge_chunk_contigs(Read_Chunk_BPtr c1rc1_chunk_bptr, Read_Chunk_BPtr c2rc2_chunk_bptr, Cigar& rc1rc2_cigar)
{
    ASSERT(rc1rc2_cigar.check(c1rc1_chunk_bptr->get_seq(), c2rc2_chunk_bptr->get_seq()));
    ASSERT(c1rc1_chunk_bptr->get_c_start() == 0 and c1rc1_chunk_bptr->get_c_end() == c1rc1_chunk_bptr->ce_bptr()->get_len());
    ASSERT(c2rc2_chunk_bptr->get_c_start() == 0 and c2rc2_chunk_bptr->get_c_end() == c2rc2_chunk_bptr->ce_bptr()->get_len());
    // do not do anything if the chunks are already mapped to the same contig
    // NOTE: with this, we are ignoring alternate mappings
    Contig_Entry_BPtr c1_ce_bptr = c1rc1_chunk_bptr->ce_bptr();
    Contig_Entry_BPtr c2_ce_bptr = c2rc2_chunk_bptr->ce_bptr();
    if (c1_ce_bptr == c2_ce_bptr)
    {
        return;
    }

    //cerr << "before merging read chunks:\n" << *c1rc1_chunk_cptr << *c2rc2_chunk_cptr << rc1rc2_cigar << *this;

    // contruct read chunk for the rc1->rc2 mapping
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
    for (auto rc_bref : c2_ce_bptr->chunk_cont())
    {
        Read_Chunk_BPtr c2rc_bptr = &rc_bref;
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
    ASSERT(rc1_ce_bptr->chunk_cont().size() == 1);
    rc1_ce_bptr->chunk_cont().clear_and_dispose();
    rc1_ce_bptr->mut_cont().clear_and_dispose();
    Contig_Entry_Fact::del_elem(rc1_ce_bptr);

    // rc2c2_chunk_bptr
    Contig_Entry_BPtr rc2_ce_bptr = rc2c2_chunk_bptr->ce_bptr();
    ASSERT(rc2_ce_bptr->chunk_cont().size() == 1);
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
    ASSERT(check(set< Contig_Entry_CBPtr >( { c1_ce_bptr })));
}

vector< std::tuple< Read_Chunk_BPtr, Read_Chunk_BPtr, Cigar > >
Graph::chunker(Read_Entry_BPtr re1_bptr, Read_Entry_BPtr re2_bptr, Cigar& cigar)
{
    vector< std::tuple< Read_Chunk_BPtr, Read_Chunk_BPtr, Cigar > > rc_mapping;

    /// match start in re1
    Size_Type r1_start = cigar.get_rf_start();
    /// match length in re1
    Size_Type r1_len = cigar.get_rf_len();
    /// match start in re2
    Size_Type r2_start = cigar.get_qr_start();
    /// match length in re2
    Size_Type r2_len = cigar.get_qr_len();
    /// is the match between re1 and re2 reverse-complemented
    bool r2_rc = cigar.is_reversed();
    // we repeatedly cut the read entries of either read
    // until their chunks match in the way described by the cigar string
    // keep track of read chunk mapping, and cigar transformation between them
    bool done;
    while (true)
    {
        // after every graph modification, restart at the beginning
        done = true;
        rc_mapping.clear();

        /// re1 part matched so far: [r1_start, r1_pos)
        Size_Type r1_pos = r1_start;
        /// re2 part matched so far: [r2_start, r2_pos) or [r2_pos, r2_end) if rc
        Size_Type r2_pos = (not r2_rc? r2_start : r2_start + r2_len);
        /// cigar ops matched so far: [0, op_start)
        size_t op_start = 0;
        while (op_start < cigar.get_n_ops())
        {
            /// next chunk in re1
            Read_Chunk_BPtr rc1_bptr = re1_bptr->chunk_cont().get_chunk_with_pos(r1_pos).unconst();
            /// re1 position relative to rc1_r_start
            Size_Type rc1_offset = r1_pos - rc1_bptr->get_r_start();
            /// re1 length left in rc1 after partial match
            Size_Type rc1_remaining_len = rc1_bptr->get_r_len() - rc1_offset;
            /// next chunk in re2
            Read_Chunk_BPtr rc2_bptr = (not r2_rc?
                                        re2_bptr->chunk_cont().get_chunk_with_pos(r2_pos).unconst()
                                        : re2_bptr->chunk_cont().get_chunk_with_pos(r2_pos - 1).unconst());
            /// re2 position relative to rc2_r_start or rc2_r_end if rc
            Size_Type rc2_offset = (not r2_rc?
                                    r2_pos - rc2_bptr->get_r_start()
                                    : rc2_bptr->get_r_end() - r2_pos);
            /// re2 length left in rc2 after partial match
            Size_Type rc2_remaining_len = rc2_bptr->get_r_len() - rc2_offset;

            // invariant: we matched read 1 chunks before rc1
            // to read 2 chunks before/after rc2
            // using cigar ops before op_start
            ASSERT(r1_pos < r1_start + r1_len);
            ASSERT(r2_rc or r2_pos < r2_start + r2_len);
            ASSERT(not r2_rc or r2_pos > r2_start);
            ASSERT(rc1_offset == 0 or rc1_bptr->is_unmappable());
            ASSERT(rc2_offset == 0 or rc2_bptr->is_unmappable());
            ASSERT(rc1_bptr);
            ASSERT(rc1_bptr->get_r_start() + rc1_offset == r1_pos);
            ASSERT(rc1_remaining_len > 0);
            ASSERT(rc2_bptr);
            ASSERT(r2_rc or rc2_bptr->get_r_start() + rc2_offset == r2_pos);
            ASSERT(not r2_rc or rc2_bptr->get_r_end() - rc2_offset == r2_pos);
            ASSERT(rc2_remaining_len > 0);
            ASSERT(r1_pos == cigar.get_rf_offset(op_start));
            ASSERT(r2_pos == cigar.get_qr_offset(op_start));

            // advance past cigar ops until either chunk ends
            size_t op_end = op_start + 1;
            while (op_end < cigar.get_n_ops()
                    and (cigar.get_rf_sub_len(op_start, op_end) < rc1_remaining_len
                         or (cigar.get_rf_sub_len(op_start, op_end) == rc1_remaining_len and cigar.is_insertion(op_end)))
                    and (cigar.get_qr_sub_len(op_start, op_end) < rc2_remaining_len
                         or (cigar.get_qr_sub_len(op_start, op_end) == rc2_remaining_len and cigar.is_deletion(op_end))))
            {
                ++op_end;
            }
            // stop conditions: can be derived by inspecting loop
            ASSERT(// past all ops
                   op_end == cigar.get_n_ops()
                   // or current chunks perfectly matched by cigar
                   or (cigar.get_rf_sub_len(op_start, op_end) == rc1_remaining_len
                       and cigar.get_qr_sub_len(op_start, op_end) == rc2_remaining_len)
                   // or cigar op extends past rc1
                   or cigar.get_rf_sub_len(op_start, op_end) > rc1_remaining_len
                   // or cigar op extends past rc2
                   or cigar.get_qr_sub_len(op_start, op_end) > rc2_remaining_len
                   // or rc1 mapped fully, but rc2 has non-insertion op at the end
                   or (cigar.get_rf_sub_len(op_start, op_end) == rc1_remaining_len
                       and cigar.get_qr_sub_len(op_start, op_end) < rc2_remaining_len
                       and not cigar.is_insertion(op_end))
                   // or rc2 mapped fully, but rc1 has non-deletion op at the end
                   or (cigar.get_rf_sub_len(op_start, op_end) < rc1_remaining_len
                       and cigar.get_qr_sub_len(op_start, op_end) == rc2_remaining_len
                       and not cigar.is_deletion(op_end)));
            // we previously cut the chunks at the ends of the mapping
            // cuts don't succeed only if those chunks are unmappable
            ASSERT(op_end < cigar.get_n_ops()
                   or cigar.get_rf_sub_len(op_start, op_end) >= rc1_remaining_len
                   or rc1_bptr->is_unmappable());
            ASSERT(op_end < cigar.get_n_ops()
                   or cigar.get_qr_sub_len(op_start, op_end) >= rc2_remaining_len
                   or rc2_bptr->is_unmappable());
            if (cigar.get_rf_sub_len(op_start, op_end) > rc1_remaining_len)
            {
                // the first inequality trivially holds with <=
                // but it can be shown it holds in fact with <
                ASSERT(cigar.get_rf_offset(op_end - 1) < rc1_bptr->get_r_end()
                       and rc1_bptr->get_r_end() < cigar.get_rf_offset(op_end));
            }
            if (cigar.get_qr_sub_len(op_start, op_end) > rc2_remaining_len)
            {
                // as with rf, the inequalities involving op_end-1 tivially hold with <=
                // but it can be shown they hold with <
                ASSERT(r2_rc or (cigar.get_qr_offset(op_end - 1) < rc2_bptr->get_r_end()
                       and rc2_bptr->get_r_end() < cigar.get_qr_offset(op_end)));
                ASSERT(not r2_rc or (cigar.get_qr_offset(op_end) < rc2_bptr->get_r_start()
                       and rc2_bptr->get_r_start() < cigar.get_qr_offset(op_end - 1)));
            }

            // check if either chunk ended during the last cigar op
            if (cigar.get_rf_sub_len(op_start, op_end) > rc1_remaining_len
                or cigar.get_qr_sub_len(op_start, op_end) > rc2_remaining_len)
            {
                // find out which of the 2 ended earlier, cut the cigar op at that position
                Size_Type r1_break_len = 0;
                if (cigar.get_rf_sub_len(op_start, op_end) > rc1_remaining_len and not cigar.is_insertion(op_end - 1))
                {
                    r1_break_len = cigar.get_rf_op_prefix_len(op_end - 1, rc1_bptr->get_r_end());
                    ASSERT(0 < r1_break_len and r1_break_len < cigar.get_rf_op_len(op_end - 1));
                }

                Size_Type r2_break_len = 0;
                if (cigar.get_qr_sub_len(op_start, op_end) > rc2_remaining_len and not cigar.is_deletion(op_end -1))
                {
                    r2_break_len = cigar.get_qr_op_prefix_len(op_end - 1, (not r2_rc?
                                                                           rc2_bptr->get_r_end()
                                                                           : rc2_bptr->get_r_start()));
                    ASSERT(0 < r2_break_len and r2_break_len < cigar.get_qr_op_len(op_end - 1));
                }

                ASSERT(r1_break_len > 0 or r2_break_len > 0);
                cigar.cut_op(op_end - 1, (r1_break_len == 0?
                                          r2_break_len
                                          : (r2_break_len == 0?
                                             r1_break_len
                                             : min(r1_break_len, r2_break_len))));
            }

            ASSERT(cigar.get_rf_sub_len(op_start, op_end) <= rc1_remaining_len
                   and cigar.get_qr_sub_len(op_start, op_end) <= rc2_remaining_len);
            // now we are sure the op ends on at least one of the read chunk boundaries
            // unless these are the last chunks and they are both unmappable
            if (op_end == cigar.get_n_ops() and rc1_bptr->is_unmappable() and rc2_bptr->is_unmappable())
            {
                ASSERT(done);
                break;
            }
            ASSERT(cigar.get_rf_sub_len(op_start, op_end) == rc1_remaining_len
                   or cigar.get_qr_sub_len(op_start, op_end) == rc2_remaining_len);
            // if it doesn't end on both, we might need to cut the other
            if (cigar.get_qr_sub_len(op_start, op_end) < rc2_remaining_len)
            {
                // it follows from main stop condition that the next op is not an insertion
                ASSERT(op_end == cigar.get_n_ops() or not cigar.is_insertion(op_end));
                if (cigar.get_qr_sub_len(op_start, op_end) == 0)
                {
                    // no progress on rc2: rc1 is mapped entirely to a deletion
                    // move on to next chunk on read 1
                    ASSERT(cigar.get_rf_sub_len(op_start, op_end) == rc1_remaining_len
                           and cigar.get_rf_sub_len(op_start, op_end) > 0);
                    r1_pos = rc1_bptr->get_r_end();
                    op_start = op_end;
                    continue;
                }
                else if (rc2_bptr->is_unmappable())
                {
                    // rc2 cannot be cut
                    if (not rc1_bptr->is_unmappable())
                    {
                        // rc1 is mapped to an unmappable chunk; we unmap it
                        unmap_chunk(rc1_bptr);
                        done = false;
                        break;
                    }
                    else // both rc1 and rc2 are unmappable
                    {
                        // we advance both
                        r1_pos = rc1_bptr->get_r_end();
                        r2_pos = (not r2_rc?
                                  r2_pos + cigar.get_qr_sub_len(op_start, op_end)
                                  : r2_pos - cigar.get_qr_sub_len(op_start, op_end));
                        op_start = op_end;
                        continue;
                    }
                }
                else
                {
                    cut_read_entry(re2_bptr, cigar.get_qr_offset(op_end), false);
                    done = false;
                    break;
                }
            }
            if (cigar.get_rf_sub_len(op_start, op_end) < rc1_remaining_len)
            {
                // it follows from main stop condition that the next op is not a deletion
                ASSERT(op_end == cigar.get_n_ops() or not cigar.is_deletion(op_end));
                if (cigar.get_rf_sub_len(op_start, op_end) == 0)
                {
                    // no progress on rc1: rc2 mapped entirely to an insertion
                    // move on to next chunk on read 2
                    ASSERT(cigar.get_qr_sub_len(op_start, op_end) == rc2_remaining_len
                           and cigar.get_qr_sub_len(op_start, op_end) > 0);
                    r2_pos = (not r2_rc? rc2_bptr->get_r_end() : rc2_bptr->get_r_start());
                    op_start = op_end;
                    continue;
                }
                else if (rc1_bptr->is_unmappable())
                {
                    // rc1 cannot be cut
                    if (not rc2_bptr->is_unmappable())
                    {
                        // rc2 is mapped to an unmappable chunk; we unmap it
                        unmap_chunk(rc2_bptr);
                        done = false;
                        break;
                    }
                    else // both rc1 and rc2 are unmappable
                    {
                        // we advance both
                        r1_pos += cigar.get_rf_sub_len(op_start, op_end);
                        r2_pos = (not r2_rc? rc2_bptr->get_r_end() : rc2_bptr->get_r_start());
                        op_start = op_end;
                        continue;
                    }
                }
                else
                {
                    cut_read_entry(re1_bptr, cigar.get_rf_offset(op_end), false);
                    done = false;
                    break;
                }
            }
            // reached when both rc1 and rc2 end at the current cigar op
            ASSERT(cigar.get_rf_sub_len(op_start, op_end) == rc1_remaining_len
                   and cigar.get_qr_sub_len(op_start, op_end) == rc2_remaining_len);
            if (rc1_bptr->is_unmappable() and not rc2_bptr->is_unmappable())
            {
                unmap_chunk(rc2_bptr);
                done = false;
                break;
            }
            if (rc2_bptr->is_unmappable() and not rc1_bptr->is_unmappable())
            {
                unmap_chunk(rc1_bptr);
                done = false;
                break;
            }
            ASSERT(rc1_bptr->is_unmappable() == rc2_bptr->is_unmappable());
            if (not rc1_bptr->is_unmappable())
            {
                // add read chunk mapping
                rc_mapping.push_back(std::make_tuple(rc1_bptr, rc2_bptr, cigar.substring(op_start, op_end)));
            }
            // advance both chunks and cigar
            r1_pos = rc1_bptr->get_r_end();
            r2_pos = (not r2_rc? rc2_bptr->get_r_end() : rc2_bptr->get_r_start());
            op_start = op_end;
        } // while (op_start < cigar.get_n_ops())
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
    // construct cigar object
    Cigar cigar(cigar_string, r2_rc, r1_start, r2_start);
    ASSERT(r1_len == cigar.get_rf_len());
    ASSERT(r2_len == cigar.get_qr_len());

    // discard indels at either end of the cigar string
    while (cigar.get_n_ops() > 0 and not cigar.is_match(0))
    {
        cigar = cigar.substring(1, cigar.get_n_ops() - 1);
    }
    while (cigar.get_n_ops() > 0 and not cigar.is_match(cigar.get_n_ops() - 1))
    {
        cigar = cigar.substring(0, cigar.get_n_ops() - 1);
    }
    r1_start = cigar.get_rf_start();
    r1_len = cigar.get_rf_len();
    r2_start = cigar.get_qr_start();
    r2_len = cigar.get_qr_len();

    // if no match is left, nothing to do
    if (r1_len == 0)
    {
        return;
    }

    Read_Entry_BPtr re1_bptr = re_cont().find(r1_name).unconst();
    ASSERT(re1_bptr != nullptr);
    Read_Entry_BPtr re2_bptr = re_cont().find(r2_name).unconst();
    ASSERT(re2_bptr != nullptr);

    string r1_seq = re1_bptr->get_seq();
    string r2_seq = re2_bptr->get_seq();
    /*
    cerr << indent::tab << "adding overlap:" << indent::inc
         << indent::nl << "re1: " << r1_seq.substr(r1_start, r1_len)
         << indent::nl << "re2: " << (not cigar.is_reversed()? r2_seq.substr(r2_start, r2_len) : reverseComplement(r2_seq.substr(r2_start, r2_len)))
         << indent::nl << "initial cigar:\n" << indent::inc << cigar << indent::dec;
    */
    cigar.disambiguate(r1_seq.substr(r1_start, r1_len), r2_seq.substr(r2_start, r2_len));
    //cerr << indent::tab << "disambiguated cigar:\n" << indent::inc << cigar << indent::dec << indent::dec;

    // cut r1 & r2 at the ends of the match region
    // NOTE: unmappable chunks are never cut
    cut_read_entry(re1_bptr, r1_start, true);
    cut_read_entry(re1_bptr, r1_start + r1_len, true);
    cut_read_entry(re2_bptr, r2_start, true);
    cut_read_entry(re2_bptr, r2_start + r2_len, true);

    auto rc_mapping = chunker(re1_bptr, re2_bptr, cigar);

    // reached when we have a complete rc map
    for (auto& tmp : rc_mapping)
    {
        Read_Chunk_BPtr rc1_bptr;
        Read_Chunk_BPtr rc2_bptr;
        Cigar rc1rc2_cigar;
        std::tie(rc1_bptr, rc2_bptr, rc1rc2_cigar) = std::move(tmp);
        merge_chunk_contigs(rc1_bptr, rc2_bptr, rc1rc2_cigar);
    }

    //TODO
/*
    // check for unmappable chunks in the contigs recently merged
    scan_read_for_unmappable_chunks(re1_bptr, r1_start, r1_start + r1_len);

    ASSERT(check(set< Read_Entry_CBPtr >( { re1_bptr, re2_bptr })));

    if (global::merge_contigs_at_each_step)
    {
        //cerr << "before merging:\n" << *this;
        merge_read_contigs(re1_bptr);
        //cerr << "after merging:\n" << *this;
        ASSERT(check(set< Read_Entry_CBPtr >( { re1_bptr, re2_bptr })));
    }
*/
}

bool Graph::cat_contigs(Contig_Entry_BPtr ce_bptr, bool c_right)
{
    auto tmp = ce_bptr->can_cat_dir(c_right);
    Contig_Entry_BPtr ce_next_bptr = std::get<0>(tmp).unconst();
    bool same_orientation = std::get<1>(tmp);
    vector< Read_Chunk_CBPtr > rc_cbptr_cont = std::move(std::get<2>(tmp));
    if (not ce_next_bptr)
    {
        return false;
    }
    if (ce_bptr->is_unmappable() != ce_next_bptr->is_unmappable())
    {
        // don't merge mappable with unmappable
        return false;
    }

    if (not same_orientation)
    {
        // we reverse one of the contigs
        // reversal modifications are done in-place, so chunk vectors are still valid
        ce_bptr->reverse();
        c_right = not c_right;
        same_orientation = true;
    }
    if (not c_right)
    {
        swap(ce_bptr, ce_next_bptr);
        c_right = true;
        auto tmp2 = ce_bptr->can_cat_dir(true);
        rc_cbptr_cont = std::move(std::get<2>(tmp2));
        ASSERT(not rc_cbptr_cont.empty());
    }
    // at this point:
    // - contigs are in the same orientation
    // - ce_bptr is the start of the merge and ce_next_bptr is the tail
    // - rc_cbptr_cont contains chunks from ce_bptr that go into ce_next_bptr

    //cerr << "merging contigs\n" << *ce_cptr << "and\n" << *ce_next_cptr;

    ce_cont().erase(ce_next_bptr);
    Contig_Entry::cat_c_right(ce_bptr, ce_next_bptr, rc_cbptr_cont);

    //cerr << "merging contigs result\n" << *ce_cptr;

    return true;
} // cat_contigs

/*
void Graph::merge_read_contigs(const Read_Entry* re_cptr)
{
    bool done = false;
    while (not done)
    {
        done = true;
        for (auto rc_it = re_cptr->get_chunk_cont().begin(); rc_it != re_cptr->get_chunk_cont().end(); ++rc_it)
        {
            Read_Chunk_CPtr rc_cptr = &*rc_it;
            if (rc_cptr->is_unmappable())
                continue;
            const Contig_Entry* ce_cptr = rc_cptr->get_ce_ptr();
            // first try to merge first contig to the left or read start
            if (rc_it == re_cptr->get_chunk_cont().begin())
            {
                if (try_merge_contig(ce_cptr, rc_cptr->get_rc()))
                {
                    done = false;
                    break;
                }
            }
            // then every contig to the right of its chunk
            if (try_merge_contig(ce_cptr, not rc_cptr->get_rc()))
            {
                done = false;
                break;
            }
        }
    }
}

void Graph::merge_all_read_contigs()
{
    for (auto re_it = _re_cont.begin(); re_it != _re_cont.end(); ++re_it)
    {
        const Read_Entry* re_cptr = &*re_it;
        merge_read_contigs(re_cptr);
    }
}
*/

void Graph::unmap_chunk(Read_Chunk_BPtr rc_bptr)
{
    // we save read entry and offset, and recompute chunk after graph modifications
    Read_Entry_BPtr re_bptr = rc_bptr->re_bptr();
    Size_Type rc_start = rc_bptr->get_r_start();
    Size_Type rc_end = rc_bptr->get_r_end();
    Contig_Entry_BPtr ce_bptr = rc_bptr->ce_bptr();

    // trim contig entry to the extent of the read chunk
    cut_contig_entry(ce_bptr, rc_bptr->get_c_start(), NULL);
    //rc_bptr = re_bptr->chunk_cont().get_chunk_with_pos(rc_start);
    ASSERT(rc_bptr->get_r_start() == rc_start);
    ASSERT(rc_bptr->get_r_end() == rc_end);
    ce_bptr = rc_bptr->ce_bptr();
    Mutation_CBPtr last_ins_cbptr = NULL;
    if (rc_bptr->mut_ptr_cont().size() > 0
        and rc_bptr->mut_ptr_cont().rbegin()->mut_cbptr()->is_ins()
        and rc_bptr->mut_ptr_cont().rbegin()->mut_cbptr()->get_start() == rc_bptr->get_c_end())
    {
        last_ins_cbptr = rc_bptr->mut_ptr_cont().rbegin()->mut_cbptr();
    }
    cut_contig_entry(ce_bptr, rc_bptr->get_c_end(), last_ins_cbptr);
    //rc_cptr = re_cptr->get_chunk_with_pos(rc_start);
    ASSERT(rc_bptr->get_r_start() == rc_start);
    ASSERT(rc_bptr->get_r_end() == rc_end);
    ce_bptr = rc_bptr->ce_bptr();

    // at this point, the chunk should span the entire contig
    ASSERT(rc_bptr->get_c_start() == 0 and rc_bptr->get_c_end() == ce_bptr->get_len());

    // chunks mapping to this contig are remapped to individual contigs, with unmappable flag set
    // we also save them in a list using Read_Entry and offsets
    vector< std::tuple< Read_Entry_BPtr, Size_Type, Size_Type > > l;
    ce_bptr->chunk_cont().clear_and_dispose([&] (Read_Chunk_BPtr rc_bptr) {
        l.push_back(std::make_tuple(rc_bptr->re_bptr(), rc_bptr->get_r_start(), rc_bptr->get_r_end()));
        Contig_Entry_BPtr ce_new_bptr = Contig_Entry_Fact::new_elem(rc_bptr->get_seq());
        rc_bptr->mut_ptr_cont().clear_and_dispose();
        rc_bptr->ce_bptr() = ce_new_bptr;
        rc_bptr->set_unmappable();
        ce_new_bptr->chunk_cont().insert(rc_bptr);
        ce_cont().insert(ce_new_bptr);
    });
    // done with old Contig_Entry
    ce_bptr->mut_cont().clear_and_dispose();
    ce_cont().erase(ce_bptr);

    // for all read chunks made unmappable, try to merge them with nearby unmappable chunks
    set< Read_Entry_CBPtr > re_set;
    for (auto& e : l) {
        std::tie(re_bptr, rc_start, rc_end) = e;
        re_set.insert(re_bptr);
        rc_bptr = re_bptr->chunk_cont().get_chunk_with_pos(rc_start).unconst();
        ASSERT(rc_bptr->is_unmappable());
        ASSERT(rc_bptr->get_r_start() <= rc_start and rc_end <= rc_bptr->get_r_end());
        if (rc_bptr->get_r_start() != rc_start or rc_bptr->get_r_end() != rc_end)
        {
            // already been merged with other unmappable chunks
            continue;
        }
        extend_unmapped_chunk(re_bptr, rc_start, rc_end);
    }
    ASSERT(check(re_set));
} // unmap_chunk

void Graph::extend_unmapped_chunk_dir(Read_Entry_BPtr re_bptr, Size_Type pos, bool r_right)
{
    bool r_left = not r_right;
    Size_Type leftover_bp = (r_left? pos : re_bptr->get_len() - pos);
    while (leftover_bp > 0)
    {
        // calls in previous iterations might extend the unmapped region
        // we first recompute the chunks based solely on: re_cptr & unmap_start
        Read_Chunk_BPtr rc_bptr = re_bptr->chunk_cont().get_chunk_with_pos(r_left? pos : pos - 1).unconst();
        ASSERT(rc_bptr->is_unmappable());
        ASSERT(not r_left or rc_bptr->get_r_start() <= pos);
        ASSERT(not r_right or pos <= rc_bptr->get_r_end());
        pos = (r_left? rc_bptr->get_r_start() : rc_bptr->get_r_end());
        leftover_bp = (r_left? pos : re_bptr->get_len() - pos);
        if (leftover_bp == 0)
        {
            break;
        }
        Read_Chunk_BPtr next_rc_bptr = re_bptr->chunk_cont().get_sibling(rc_bptr, true, r_right).unconst();

        if (next_rc_bptr->is_unmappable())
        {
            // consecutive unmappable chunks; merge contigs
            ASSERT(not rc_bptr->get_rc());
            ASSERT(not next_rc_bptr->get_rc());
            pos = (r_left? next_rc_bptr->get_r_start() : next_rc_bptr->get_r_end());
            bool success = cat_contigs(r_left? next_rc_bptr->ce_bptr() : rc_bptr->ce_bptr(), true);
            ASSERT(success);
            continue;
        }
        else // next_rc_bptr not unmappable
        {
            // next chunk is mappable
            if (leftover_bp <= global::unmap_trigger_len)
            {
                unmap_chunk(next_rc_bptr);
                continue;
            }
            else // global::unmap_trigger_len < leftover_bp
            {
                // check there is a minimum of mappable bp
                // if not, unmap small intermediate chunks and merge them
                // skip small mappable chunks
                list< std::tuple< Size_Type, Size_Type > > skipped_chunks;
                Size_Type skipped_len = 0;
                while (not next_rc_bptr->is_unmappable()
                       and skipped_len + next_rc_bptr->get_r_len() <= global::unmap_trigger_len)
                {
                    skipped_chunks.push_back(std::make_tuple(next_rc_bptr->get_r_start(), next_rc_bptr->get_r_end()));
                    skipped_len += next_rc_bptr->get_r_len();
                    next_rc_bptr = re_bptr->chunk_cont().get_sibling(next_rc_bptr, true, r_right).unconst();
                    // since skipped_len < unmap_trigger_len < leftover_bp
                    ASSERT(next_rc_bptr);
                }
                ASSERT(skipped_len <= global::unmap_trigger_len);
                if (next_rc_bptr->is_unmappable())
                {
                    // unmap all skipped chunks (if they are not modified in previous iterations)
                    for (auto it = skipped_chunks.begin(); it != skipped_chunks.end(); ++it)
                    {
                        Size_Type rc_start;
                        Size_Type rc_end;
                        std::tie(rc_start, rc_end) = *it;
                        rc_bptr = re_bptr->chunk_cont().get_chunk_with_pos(rc_start).unconst();
                        if (rc_bptr->get_r_start() != rc_start or rc_bptr->get_r_end() != rc_end or rc_bptr->is_unmappable())
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
    ASSERT(check(set< Read_Entry_CBPtr >( { re_bptr })));
} // extend_unmapped_chunk_dir

void Graph::extend_unmapped_chunk(Read_Entry_BPtr re_bptr, Size_Type rc_start, Size_Type rc_end)
{
    extend_unmapped_chunk_dir(re_bptr, rc_start, false);
    extend_unmapped_chunk_dir(re_bptr, rc_end, true);
}

/*
void Graph::scan_read_for_unmappable_chunks(const Read_Entry* re_cptr, Size_Type rc_start, Size_Type rc_end)
{
    Size_Type pos = rc_start;
    while (pos < rc_end)
    {
        Read_Chunk_CPtr rc_cptr = re_cptr->get_chunk_with_pos(pos);
        ASSERT(rc_cptr != NULL);
        ASSERT(rc_cptr->get_r_start() <= pos);
        ASSERT(pos < rc_cptr->get_r_end());
        pos = rc_cptr->get_r_end();
        scan_contig_for_unmappable_chunks(rc_cptr->get_ce_ptr());
    }
}

void Graph::scan_contig_for_unmappable_chunks(const Contig_Entry* ce_cptr)
{
    // if 2 chunks of the same read are mapped to overlapping regions of this contig, unmap both
    list< const Read_Entry* > re_list;
    list< std::tuple< Size_Type, Size_Type > > pos_list;
    for (size_t i = 0; i < ce_cptr->get_chunk_cptr_cont().size(); ++i)
    {
        Read_Chunk_CPtr rc1_cptr = ce_cptr->get_chunk_cptr_cont()[i];
        for (size_t j = 0; j < i; ++j)
        {
            Read_Chunk_CPtr rc2_cptr = ce_cptr->get_chunk_cptr_cont()[j];
            if (rc1_cptr->get_re_ptr() == rc2_cptr->get_re_ptr())
                / *
                and ((rc1_cptr->get_c_start() <= rc2_cptr->get_c_start() and rc2_cptr->get_c_start() < rc1_cptr->get_c_end())
                     or (rc2_cptr->get_c_start() <= rc1_cptr->get_c_start() and rc1_cptr->get_c_start() < rc2_cptr->get_c_end())
                     or (rc1_cptr->get_c_end() == rc2_cptr->get_c_start()
                         and rc1_cptr->get_mut_ptr_cont().size() > 0
                         and rc2_cptr->get_mut_ptr_cont().size() > 0
                         and rc1_cptr->get_mut_ptr_cont().back()->is_ins()
                         and rc1_cptr->get_mut_ptr_cont().back() == rc2_cptr->get_mut_ptr_cont().front())
                     or (rc2_cptr->get_c_end() == rc1_cptr->get_c_start()
                         and rc2_cptr->get_mut_ptr_cont().size() > 0
                         and rc1_cptr->get_mut_ptr_cont().size() > 0
                         and rc2_cptr->get_mut_ptr_cont().back()->is_ins()
                         and rc2_cptr->get_mut_ptr_cont().back() == rc1_cptr->get_mut_ptr_cont().front())))
                * /
            {
                // overlapping chunks from the same read
                re_list.push_back(rc1_cptr->get_re_ptr());
                pos_list.push_back(std::make_tuple(rc1_cptr->get_r_start(), rc1_cptr->get_r_end()));
                //re_list.push_back(rc2_cptr->get_re_ptr());
                //pos_list.push_back(std::make_tuple(rc2_cptr->get_r_start(), rc2_cptr->get_r_end()));
                / *
                cerr << "overlapping chunks from same read: " << rc1_cptr->get_re_ptr()->get_name() << ' '
                << rc1_cptr->get_r_start() << ' ' << rc1_cptr->get_r_end() << ' '
                << rc2_cptr->get_r_start() << ' ' << rc2_cptr->get_r_end() << '\n'
                << *rc1_cptr << *rc2_cptr << *ce_cptr;
                * /
            }
        }
    }
    auto re_it = re_list.begin();
    auto pos_it = pos_list.begin();
    for (; re_it != re_list.end() && pos_it != pos_list.end(); ++re_it, ++pos_it)
    {
        const Read_Entry* re_cptr;
        Size_Type rc_start;
        Size_Type rc_end;
        re_cptr = *re_it;
        std::tie(rc_start, rc_end) = *pos_it;
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
}

void Graph::clear_contig_colours()
{
    for (auto ce_it = _ce_cont.begin(); ce_it != _ce_cont.end(); ++ce_it)
    {
        const Contig_Entry* ce_cptr = &*ce_it;
        modify_contig_entry(ce_cptr, [] (Contig_Entry& ce) {
            ce.colour() = 0;
        });
    }
}

void Graph::clear_contig_visit()
{
    for (auto ce_it = _ce_cont.begin(); ce_it != _ce_cont.end(); ++ce_it)
    {
        const Contig_Entry* ce_cptr = &*ce_it;
        modify_contig_entry(ce_cptr, [] (Contig_Entry& ce) {
            ce.colour() &= ~0x3;
        });
    }
}

std::tuple< Size_Type, bool> Graph::visit_contig(const Contig_Entry* ce_cptr, bool dir)
{
    if ((ce_cptr->get_colour() & 0x1) != 0)
    {
        return std::make_tuple(0, (ce_cptr->get_colour() & (dir == false? 0x4 : 0x8)) != 0);
    }
    auto neighbours_sptr = ce_cptr->get_neighbours(dir);
    ASSERT(neighbours_sptr->size() > 0);
    if (neighbours_sptr->size() > 1)
    {
        modify_contig_entry(ce_cptr, [&] (Contig_Entry& ce) {
            ce.colour() |= (dir == false? 0x4 : 0x8);
        });
        return std::make_tuple(0, true);
    }
    ASSERT(neighbours_sptr->size() == 1);
    modify_contig_entry(ce_cptr, [] (Contig_Entry& ce) {
        ce.colour() |= 0x1;
    });
    neighbours_sptr = ce_cptr->get_neighbours(not dir);
    Size_Type tmp = 0;
    bool dead_end = false;
    if (neighbours_sptr->size() == 1)
    {
        const Contig_Entry* ce_next_cptr;
        bool flip;
        std::tie(ce_next_cptr, flip) = neighbours_sptr->begin()->first;
        std::tie(tmp, dead_end) = visit_contig(ce_next_cptr, (not dir) == flip);
    }
    if (neighbours_sptr->size() != 1 or dead_end)
    {
        modify_contig_entry(ce_cptr, [&] (Contig_Entry& ce) {
            ce.colour() |= (not dir == false? 0x4 : 0x8);
        });
    }
    return std::make_tuple(ce_cptr->get_len() + tmp, false);
}

void Graph::dfs_scontig(const Contig_Entry* ce_cptr, bool ce_endpoint, bool dir,
                        list< std::tuple< const Contig_Entry*, size_t, size_t, bool > >& l, bool& cycle)
{
    while (true)
    {
        if ((ce_cptr->get_colour() & 0x1) != 0)
        {
            cycle = (ce_cptr->get_colour() & (ce_endpoint == false? 0x4 : 0x8)) == 0;
            break;
        }
        auto neighbours_sptr = ce_cptr->get_neighbours(ce_endpoint);
        ASSERT(neighbours_sptr->size() > 0);
        if (neighbours_sptr->size() > 1)
        {
            modify_contig_entry(ce_cptr, [&] (Contig_Entry& ce) {
                ce.colour() |= (ce_endpoint == false? 0x4 : 0x8);
            });
            cycle = false;
            break;
        }
        ASSERT(neighbours_sptr->size() == 1);
        modify_contig_entry(ce_cptr, [] (Contig_Entry& ce) {
            ce.colour() |= 0x1;
        });
        {
            auto t = std::make_tuple(ce_cptr, ce_cptr->get_len(), ce_cptr->get_len(), dir == ce_endpoint);
            not dir? l.push_front(t) : l.push_back(t);
        }
        neighbours_sptr = ce_cptr->get_neighbours(not ce_endpoint);
        if (neighbours_sptr->size() == 1)
        {
            const Contig_Entry* ce_next_cptr;
            bool flip;
            unsigned int tmp_support;
            Size_Type min_skipped_len;
            Size_Type max_skipped_len;
            std::tie(ce_next_cptr, flip) = neighbours_sptr->begin()->first;
            std::tie(tmp_support, min_skipped_len, max_skipped_len) = neighbours_sptr->begin()->second;
            if (max_skipped_len > 0)
            {
                auto t = std::make_tuple((const Contig_Entry*)NULL, min_skipped_len, max_skipped_len, false);
                not dir? l.push_front(t) : l.push_back(t);
            }
            ce_cptr = ce_next_cptr;
            ce_endpoint = not ce_endpoint == flip;
        }
        else
        {
            modify_contig_entry(ce_cptr, [&] (Contig_Entry& ce) {
                ce.colour() |= (not ce_endpoint == false? 0x4 : 0x8);
            });
            cycle = false;
            break;
        }
    }
}

void Graph::print_supercontig_lengths(ostream& os)
{
    for (auto ce_it = _ce_cont.begin(); ce_it != _ce_cont.end(); ++ce_it)
    {
        const Contig_Entry* ce_cptr = &*ce_it;
        if (ce_cptr->is_unmappable())
        {
            continue;
        }
        ASSERT((ce_cptr->get_colour() & 0x2) == 0);
        if ((ce_cptr->get_colour() & 0x1) != 0)
        {
            continue;
        }
        // mark contig as visited
        modify_contig_entry(ce_cptr, [] (Contig_Entry& ce) {
            ce.colour() |= 0x1;
        });
        Size_Type supercontig_left_len = 0;
        Size_Type supercontig_right_len = 0;
        const Contig_Entry* ce_next_cptr;
        bool flip;
        bool dead_end = false;
        auto neighbours_sptr = ce_cptr->get_neighbours(false);
        if (neighbours_sptr->size() == 1)
        {
            std::tie(ce_next_cptr, flip) = neighbours_sptr->begin()->first;
            std::tie(supercontig_left_len, dead_end) = visit_contig(ce_next_cptr, false == flip);
        }
        if (neighbours_sptr->size() != 1 or dead_end)
        {
            modify_contig_entry(ce_cptr, [] (Contig_Entry& ce) {
                ce.colour() |= 0x4;
            });
        }
        neighbours_sptr = ce_cptr->get_neighbours(true);
        dead_end = false;
        if (neighbours_sptr->size() == 1)
        {
            std::tie(ce_next_cptr, flip) = neighbours_sptr->begin()->first;
            std::tie(supercontig_right_len, dead_end) = visit_contig(ce_next_cptr, true == flip);
        }
        if (neighbours_sptr->size() != 1 or dead_end)
        {
            modify_contig_entry(ce_cptr, [] (Contig_Entry& ce) {
                ce.colour() |= 0x8;
            });
        }
        os << supercontig_left_len + supercontig_right_len + ce_cptr->get_len() << '\n';
    }
}

void Graph::print_supercontig_lengths_2(ostream& os)
{
    for (auto ce_it = _ce_cont.begin(); ce_it != _ce_cont.end(); ++ce_it)
    {
        const Contig_Entry* ce_cptr = &*ce_it;
        if (ce_cptr->is_unmappable())
        {
            continue;
        }
        ASSERT((ce_cptr->get_colour() & 0x2) == 0);
        if ((ce_cptr->get_colour() & 0x1) != 0)
        {
            continue;
        }
        // mark contig as visited
        modify_contig_entry(ce_cptr, [] (Contig_Entry& ce) {
            ce.colour() |= 0x1;
        });
        list< std::tuple< const Contig_Entry*, Size_Type, Size_Type, bool > > l[2];
        bool cycle[2];
        shared_ptr< map< std::tuple< const Contig_Entry*, bool >, std::tuple< unsigned int, Size_Type, Size_Type > > >
        neighbours_sptr[2];
        for (int side = 0; side < 2; ++side)
        {
            neighbours_sptr[side] = ce_cptr->get_neighbours(side == 1);
            if (neighbours_sptr[side]->size() != 1)
            {
                modify_contig_entry(ce_cptr, [&] (Contig_Entry& ce) {
                    ce.colour() |= (side == 0? 0x4 : 0x8);
                });
            }
        }
        for (int side = 0; side < 2; ++side)
        {
            cycle[side] = false;
            if (neighbours_sptr[side]->size() == 1)
            {
                const Contig_Entry* ce_next_cptr;
                bool flip;
                std::tie(ce_next_cptr, flip) = neighbours_sptr[side]->begin()->first;
                unsigned int tmp_support;
                Size_Type min_skipped_len;
                Size_Type max_skipped_len;
                std::tie(tmp_support, min_skipped_len, max_skipped_len) = neighbours_sptr[side]->begin()->second;
                if (max_skipped_len > 0)
                {
                    auto t = std::make_tuple((const Contig_Entry*)NULL, min_skipped_len, max_skipped_len, false);
                    side == 0? l[side].push_front(t) : l[side].push_back(t);
                }
                dfs_scontig(ce_next_cptr, (side == 1) == flip, side == 1, l[side], cycle[side]);
            }
        }
        ASSERT(cycle[0] == cycle[1]);
        Size_Type scon_len = ce_cptr->get_len();
        Size_Type scon_min_skip_len = 0;
        Size_Type scon_max_skip_len = 0;
        for (int side = 0; side < 2; ++side)
        {
            for (auto it = l[side].begin(); it != l[side].end(); ++it)
            {
                const Contig_Entry* tmp_ce_cptr;
                Size_Type min_skipped_len;
                Size_Type max_skipped_len;
                bool flip;
                std::tie(tmp_ce_cptr, min_skipped_len, max_skipped_len, flip) = *it;
                if (tmp_ce_cptr != NULL)
                {
                    scon_len += tmp_ce_cptr->get_len();
                }
                else
                {
                    scon_min_skip_len += min_skipped_len;
                    scon_max_skip_len += max_skipped_len;
                }
            }
        }
        os << scon_len << '\t' << scon_min_skip_len << '-' << scon_max_skip_len << '\t';
        for (int side = 0; side < 2; ++side)
        {
            if (side == 1)
            {
                if (l[0].size() > 0)
                {
                    os << ',';
                }
                os << '(' << ce_cptr->get_contig_id() << ',' << ce_cptr->get_len() << ",0)";
            }
            for (auto it = l[side].begin(); it != l[side].end(); ++it)
            {
                if (side == 1 or (side == 0 and it != l[side].begin()))
                {
                    os << ',';
                }
                const Contig_Entry* tmp_ce_cptr;
                Size_Type min_skipped_len;
                Size_Type max_skipped_len;
                bool flip;
                std::tie(tmp_ce_cptr, min_skipped_len, max_skipped_len, flip) = *it;
                if (tmp_ce_cptr != NULL)
                {
                    os << '(' << tmp_ce_cptr->get_contig_id() << ',' << tmp_ce_cptr->get_len() << ',' << int(flip) << ')';
                }
                else
                {
                    os << "(.," << min_skipped_len << ',' << max_skipped_len << ")";
                }
            }
        }
        os << '\n';
    }
}

void Graph::set_contig_ids()
{
    size_t contig_id = 0;
    for (auto ce_it = _ce_cont.begin(); ce_it != _ce_cont.end(); ++ce_it)
    {
        const Contig_Entry* ce_cptr = &*ce_it;
        modify_contig_entry(ce_cptr, [&] (Contig_Entry& ce) {
            ce.contig_id() = contig_id;
        });
        ++contig_id;
    }
}

void Graph::unmap_single_chunks()
{
    for (auto re_it = _re_cont.begin(); re_it != _re_cont.end(); ++re_it)
    {
        const Read_Entry* re_cptr = &*re_it;
        bool done = false;
        while (not done)
        {
            done = true;
            for (auto rc_it = re_cptr->get_chunk_cont().begin(); rc_it != re_cptr->get_chunk_cont().end(); ++rc_it)
            {
                Read_Chunk_CPtr rc_cptr = &*rc_it;
                const Contig_Entry* ce_cptr = rc_cptr->get_ce_ptr();
                if (not ce_cptr->is_unmappable() and ce_cptr->get_chunk_cptr_cont().size() == 1)
                {
                    unmap_chunk(rc_cptr);
                    done = false;
                    break;
                }
            }
        }
    }
}

void Graph::print_separated_het_mutations(ostream& os, size_t min_support_report, Size_Type min_separation) const
{
    for (const auto& ce : _ce_cont)
    {
        ce.print_separated_het_mutations(os, min_support_report, min_separation);
    }
}

void Graph::print_unmappable_contigs(ostream& os) const
{
    for (const auto& ce : _ce_cont)
    {
        if (ce.is_unmappable())
        {
            os << '>' << ce.get_contig_id() << '\n' << ce.get_seq() << '\n';
        }
    }
}

bool Graph::check_colours() const
{
    for (auto ce_it = _ce_cont.begin(); ce_it != _ce_cont.end(); ++ce_it)
    {
        const Contig_Entry* ce_cptr = &*ce_it;
        ASSERT(ce_cptr->check_colour(false));
        ASSERT(ce_cptr->check_colour(true));
    }
    return true;
}
*/

bool Graph::check_all() const
{
    size_t chunks_count_1 = 0;
    size_t chunks_count_2 = 0;
    // check read entry objects
    for (const auto& re_cbref : _re_cont)
    {
        Read_Entry_CBPtr re_cbptr = &re_cbref;
        ASSERT(re_cbptr->check());
        chunks_count_1 += re_cbptr->chunk_cont().size();
    }
    // check contig entry objects
    for (const auto& ce_cbref : _ce_cont)
    {
        Contig_Entry_CBPtr ce_cbptr = &ce_cbref;
        ASSERT(ce_cbptr->check());
        chunks_count_2 += ce_cbptr->chunk_cont().size();
    }
    ASSERT(chunks_count_1 == chunks_count_2);
    return true;
}

bool Graph::check(const set< Read_Entry_CBPtr >& re_set, const set< Contig_Entry_CBPtr >& ce_set) const
{
#ifdef NO_GRAPH_CHECKS
    return true;
#endif
    // compute contig entries referenced by the selected read entries
    set< Contig_Entry_CBPtr > ce_extra_set;
    for (const auto& re_cbptr : re_set)
    {
        for (const auto& rc_cbref : re_cbptr->chunk_cont())
        {
            ce_extra_set.insert(rc_cbref.raw().ce_bptr());
        }
    }

    // compute read entries referenced by selected contig entries
    set< Read_Entry_CBPtr > re_extra_set;
    for (const auto& ce_cbptr : ce_set)
    {
        for (const auto& rc_cbref : ce_cbptr->chunk_cont())
        {
            re_extra_set.insert(rc_cbref.raw().re_bptr());
        }
    }

    // check read entry objects
    for (const auto& re_cbptr : re_set)
    {
        ASSERT(re_cbptr->check());
    }
    for (const auto& re_cbptr : re_extra_set)
    {
        if (re_set.count(re_cbptr) > 0)
        {
            continue;
        }
        ASSERT(re_cbptr->check());
    }

    // check contig entry objects
    for (const auto& ce_cbptr : ce_set)
    {
        ASSERT(ce_cbptr->check());
    }
    for (const auto& ce_cbptr : ce_extra_set)
    {
        if (ce_set.count(ce_cbptr) > 0)
        {
            continue;
        }
        ASSERT(ce_cbptr->check());
    }

    return true;
}

/*
void Graph::dump_detailed_counts(ostream& os) const
{
    // First read stats
    os << "RE\tname\tlen\tnum.chunks\tchunk.lens\tcontigs\n";
    for (auto re_it = _re_cont.begin(); re_it != _re_cont.end(); ++re_it)
    {
        const Read_Entry* re_cptr = &*re_it;
        os << "RE\t" << re_cptr->get_name() << '\t' << re_cptr->get_len() << '\t' << re_cptr->get_chunk_cont().size() << '\t';
        for (auto rc_it = re_cptr->get_chunk_cont().begin(); rc_it != re_cptr->get_chunk_cont().end(); ++rc_it)
        {
            Read_Chunk_CPtr rc_cptr = &*rc_it;
            if (rc_it != re_cptr->get_chunk_cont().begin())
                os << ',';
            os << (rc_cptr->is_unmappable()? "*" : "") << rc_cptr->get_r_len();
        }
        os << '\t';
        for (auto rc_it = re_cptr->get_chunk_cont().begin(); rc_it != re_cptr->get_chunk_cont().end(); ++rc_it)
        {
            Read_Chunk_CPtr rc_cptr = &*rc_it;
            if (rc_it != re_cptr->get_chunk_cont().begin())
                os << ',';
            os << rc_cptr->get_ce_ptr()->get_contig_id();
        }
        os << '\n';
    }
    // next, contig stats
    os << "CE\tid\tlen\tnum.chunks\tbp.chunks\tnum.mut\tnum.mut.chunks"
       << "\tnum.snp\tnum.ins\tnum.del\tnum.mnp\tbp.mut"
       << "\tcovg.left\tdeg.left\tcovg.right\tdeg.right\tunmappable\tdeg.left.skip\tdeg.right.skip\tcontigs.left.skip\tcontigs.right.skip\n";
    for (auto ce_it = _ce_cont.begin(); ce_it != _ce_cont.end(); ++ce_it)
    {
        const Contig_Entry* ce_cptr = &*ce_it;
        os << "CE\t" << ce_cptr->get_contig_id() << '\t' << ce_cptr->get_len() << '\t' << ce_cptr->get_chunk_cptr_cont().size() << '\t';
        size_t reads_bp = 0;
        size_t total_muts_reads = 0;
        for (auto rc_cptr_it = ce_cptr->get_chunk_cptr_cont().begin(); rc_cptr_it != ce_cptr->get_chunk_cptr_cont().end(); ++rc_cptr_it)
        {
            Read_Chunk_CPtr rc_cptr = *rc_cptr_it;
            reads_bp += rc_cptr->get_r_len();
            total_muts_reads += rc_cptr->get_mut_ptr_cont().size();
        }
        os << reads_bp << '\t' << ce_cptr->get_mut_cont().size() << '\t' << total_muts_reads << '\t';
        size_t n_snp = 0;
        size_t n_ins = 0;
        size_t n_del = 0;
        size_t n_mnp = 0;
        size_t total_mut_bp = 0;
        for (auto mut_it = ce_cptr->get_mut_cont().begin(); mut_it != ce_cptr->get_mut_cont().end(); ++ mut_it)
        {
            Mutation_CPtr mut_cptr = &*mut_it;
            if (mut_cptr->is_snp())
                ++n_snp;
            else if (mut_cptr->is_ins())
                ++n_ins;
            else if (mut_cptr->is_del())
                ++n_del;
            else
                ++n_mnp;
            total_mut_bp += mut_cptr->get_len() + mut_cptr->get_seq_len();
        }
        os << n_snp << '\t' << n_ins << '\t' << n_del << '\t' << n_mnp << '\t' << total_mut_bp << '\t';
        size_t cnt_0;
        size_t uniq_0;
        size_t cnt_1;
        size_t uniq_1;
        std::tie(cnt_0, uniq_0, cnt_1, uniq_1) = ce_cptr->get_out_degrees();
        os << cnt_0 << '\t' << uniq_0 << '\t' << cnt_1 << '\t' << uniq_1 << '\t';
        os << (int)ce_cptr->is_unmappable() << '\t';
        if (ce_cptr->is_unmappable())
        {
            os << ".\t.\t.\t.";
        }
        else
        {
            auto neighbours_left_sptr = ce_cptr->get_neighbours(false);
            auto neighbours_right_sptr = ce_cptr->get_neighbours(true);
            os << neighbours_left_sptr->size() << '\t' << neighbours_right_sptr->size() << '\t';
            for (auto it = neighbours_left_sptr->begin(); it != neighbours_left_sptr->end(); ++it)
            {
                const Contig_Entry* tmp_ce_cptr;
                bool tmp_flip;
                if (it != neighbours_left_sptr->begin())
                {
                    os << ',';
                }
                std::tie(tmp_ce_cptr, tmp_flip) = it->first;
                unsigned int tmp_support;
                Size_Type min_skipped_len;
                Size_Type max_skipped_len;
                std::tie(tmp_support, min_skipped_len, max_skipped_len) = it->second;
                os << '(' << tmp_ce_cptr->get_contig_id() << ',' << int(tmp_flip) << ',' << tmp_support
                   << ',' << min_skipped_len << ',' << max_skipped_len << ')';
            }
            if (neighbours_left_sptr->size() == 0)
            {
                os << '.';
            }
            os << '\t';
            for (auto it = neighbours_right_sptr->begin(); it != neighbours_right_sptr->end(); ++it)
            {
                const Contig_Entry* tmp_ce_cptr;
                bool tmp_flip;
                if (it != neighbours_right_sptr->begin())
                {
                    os << ',';
                }
                std::tie(tmp_ce_cptr, tmp_flip) = it->first;
                unsigned int tmp_support;
                Size_Type min_skipped_len;
                Size_Type max_skipped_len;
                std::tie(tmp_support, min_skipped_len, max_skipped_len) = it->second;
                os << '(' << tmp_ce_cptr->get_contig_id() << ',' << int(tmp_flip) << ',' << tmp_support
                   << ',' << min_skipped_len << ',' << max_skipped_len << ')';
            }
            if (neighbours_right_sptr->size() == 0)
            {
                os << '.';
            }
        }
        os << '\n';
    }
}
*/

ostream& operator << (ostream& os, const Graph& g)
{
    os << indent::tab << "(Graph\n" << indent::inc
       << indent::tab << "Read_Entry_Cont:\n" << indent::inc;
    print_seq(os, g._re_cont, "");
    os << indent::dec << indent::tab << "Contig_Entry_Cont:\n" << indent::inc;
    print_seq(os, g._ce_cont, "");
    os << indent::dec << indent::dec << indent::tab << ")\n";
    return os;
}

} // namespace MAC
