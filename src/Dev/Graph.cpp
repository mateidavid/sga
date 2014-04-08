#include "Graph.hpp"

#include "globals.hpp"
#include "Cigar.hpp"
#include "indent.hpp"
#include "print_seq.hpp"
#include "../Util/Util.h"

using namespace std;


namespace MAC
{

/*
void Graph::erase_contig_entry(const Contig_Entry* ce_cptr)
{
    for (auto it = ce_cptr->get_chunk_cptr_cont().begin(); it != ce_cptr->get_chunk_cptr_cont().end(); ++it)
    {
        // the read chunks no longer point here
        ASSERT((*it)->get_ce_ptr() != ce_cptr);
    }
    _ce_cont.erase(_ce_cont.iterator_to(*ce_cptr));
}
*/

void Graph::add_read(const string* name_ptr, Seq_Type* seq_ptr)
{
    // first, create read entry and place it in container
    Read_Entry_BPtr re_bptr = Read_Entry_Fact::new_elem(name_ptr, seq_ptr->size());
    //cerr << indent::tab << "re:\n" << indent::inc << re << indent::dec;
    _re_cont.insert(re_bptr);

    // create contig entry and place it in container
    Contig_Entry_BPtr ce_bptr = Contig_Entry_Fact::new_elem(seq_ptr);
    //cerr << indent::tab << "ce:\n" << indent::inc << ce << indent::dec;
    ce_bptr->add_chunk(&*re_bptr->chunk_cont().begin());
    //cerr << indent::tab << "ce with chunk:\n" << indent::inc << ce << indent::dec;
    _ce_cont.insert(ce_bptr);

    // fix initial rc: assing it to contig entry
    re_bptr->chunk_cont().begin()->assign_to_contig(ce_bptr, 0, seq_ptr->size(), false, Mutation_Ptr_Cont());

    ASSERT(check(set< Read_Entry_CBPtr >( { re_bptr })));
}

bool Graph::cut_contig_entry(Contig_Entry_BPtr ce_bptr, Size_Type c_brk, Mutation_CBPtr mut_left_cbptr)
{
    (void)ce_bptr;
    (void)c_brk;
    (void)mut_left_cbptr;
/*
    if ((c_brk == 0 and mut_left_cptr == NULL)
            or (c_brk == ce_cptr->get_len() and (ce_cptr->get_mut_cont().size() == 0
                    or not ce_cptr->get_mut_cont().rbegin()->is_ins()
                    or ce_cptr->get_mut_cont().rbegin()->get_start() < ce_cptr->get_len())))
    {
        return false;
    }

    //cerr << "before cutting contig entry: " << (void*)ce_cptr << " at " << c_brk << " mut_left_cptr=" << (void*)mut_left_cptr << '\n' << *this;

    // first, split any mutations that span c_pos
    vector< const Mutation* > v = ce_cptr->get_mutations_spanning_pos(c_brk);
    for (size_t i = 0; i < v.size(); ++i)
    {
        ASSERT(v[i]->get_start() < c_brk and c_brk < v[i]->get_end());
        ASSERT(v[i] != mut_left_cptr);
        cut_mutation(ce_cptr, v[i], c_brk - v[i]->get_start(), 0);
    }

    // create new contig entry object; set string; add it to container
    const Contig_Entry* ce_new_cptr;
    ce_new_cptr = insert_contig_entry(Contig_Entry(new string(ce_cptr->get_seq().substr(c_brk))));

    // aquire second-half mutations from original contig entry object into the new one; compute mutation pointer map
    map< const Mutation*, const Mutation* > mut_cptr_map;
    auto ce_modifier = [&] (Contig_Entry& ce) {
        mut_cptr_map = ce.acquire_second_half_mutations(ce_cptr, c_brk, mut_left_cptr);
    };
    modify_contig_entry(ce_new_cptr, ce_modifier);

    // save the list of read chunk objects that must be modified
    Read_Chunk_CPtr_Cont chunk_cptr_cont = ce_cptr->get_chunk_cptr_cont();

    for (auto rc_cptr_it = chunk_cptr_cont.begin(); rc_cptr_it != chunk_cptr_cont.end(); ++rc_cptr_it)
    {
        Read_Chunk_CPtr rc_cptr = (*rc_cptr_it);
        bool move_to_rhs;
        shared_ptr< Read_Chunk > rc_new_sptr;

        // construct&apply read chunk modifier to implement split
        auto rc_modifier = [&] (Read_Chunk& rc) {
            tie(move_to_rhs, rc_new_sptr) = rc.split(c_brk, mut_cptr_map, ce_new_cptr);
        };
        modify_read_chunk(rc_cptr, rc_modifier);

        Read_Chunk_CPtr rc_rhs_cptr = NULL;
        // if a new read chunk is created, insert it in read entry object
        if (rc_new_sptr)
        {
            modify_read_entry(rc_cptr->get_re_ptr(), [&] (Read_Entry& re) {
                rc_rhs_cptr = re.add_read_chunk(rc_new_sptr.get());
            });
        }
        // if the rc must move to the rhs ce, remove it from lhs ce
        else if (move_to_rhs)
        {
            modify_contig_entry(ce_cptr, [&] (Contig_Entry& ce) {
                ce.remove_chunk(rc_cptr);
            });
            rc_rhs_cptr = rc_cptr;
        }

        // if a rc must be added to rhs ce, do it now
        if (rc_rhs_cptr != NULL)
        {
            modify_contig_entry(ce_new_cptr, [&] (Contig_Entry& ce) {
                ce.add_chunk(rc_rhs_cptr);
            });
        }
    }

    // drop mapped mutations from first part
    modify_contig_entry(ce_cptr, [&] (Contig_Entry& ce) {
        ce.drop_mutations(mut_cptr_map);
    });
    // drop the second part of the base sequence from the first part
    modify_contig_entry(ce_cptr, [&] (Contig_Entry& ce) {
        ce.drop_base_seq(c_brk);
    });
    // fix contig length for chunks left on LHS
    / * NOOO
    for (auto rc_cptr_it = ce_cptr->get_chunk_cptr_cont().begin(); rc_cptr_it != ce_cptr->get_chunk_cptr_cont().end(); ++rc_cptr_it)
    {
        modify_read_chunk(*rc_cptr_it, [&] (Read_Chunk& rc) { rc.c_len() = ce_cptr->get_len(); });
    }
    * /

    // drop unused mutations from either contig
    modify_contig_entry(ce_cptr, [] (Contig_Entry& ce) {
        ce.drop_unused_mutations();
    });
    modify_contig_entry(ce_new_cptr, [] (Contig_Entry& ce) {
        ce.drop_unused_mutations();
    });

    if (ce_new_cptr->get_chunk_cptr_cont().size() == 0)
    {
        erase_contig_entry(ce_new_cptr);
        ASSERT(check(set< const Contig_Entry* >( { ce_cptr })));
    }
    else if (ce_cptr->get_chunk_cptr_cont().size() == 0)
    {
        erase_contig_entry(ce_cptr);
        ASSERT(check(set< const Contig_Entry* >( { ce_new_cptr })));
    }
    else
    {
        //cerr << "after cutting contig entry: " << (void*)ce_cptr << '\n' << *this;
        ASSERT(check(set< const Contig_Entry* >( { ce_cptr, ce_new_cptr })));
    }
*/
    return true;
}

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
}

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

void Graph::merge_read_chunks(Read_Chunk_CPtr c1rc1_chunk_cptr, Read_Chunk_CPtr c2rc2_chunk_cptr, Cigar& rc1rc2_cigar)
{
    ASSERT(rc1rc2_cigar.check(c1rc1_chunk_cptr->get_seq(), c2rc2_chunk_cptr->get_seq()));
    ASSERT(c1rc1_chunk_cptr->get_c_start() == 0 and c1rc1_chunk_cptr->get_c_end() == c1rc1_chunk_cptr->get_ce_ptr()->get_len());
    ASSERT(c2rc2_chunk_cptr->get_c_start() == 0 and c2rc2_chunk_cptr->get_c_end() == c2rc2_chunk_cptr->get_ce_ptr()->get_len());
    // do not do anything if the chunks are already mapped to the same contig
    // NOTE: with this, we are ignoring alternate mappings
    const Contig_Entry* c1_ce_cptr = c1rc1_chunk_cptr->get_ce_ptr();
    const Contig_Entry* c2_ce_cptr = c2rc2_chunk_cptr->get_ce_ptr();
    if (c1_ce_cptr == c2_ce_cptr)
        return;

    //cerr << "before merging read chunks:\n" << *c1rc1_chunk_cptr << *c2rc2_chunk_cptr << rc1rc2_cigar << *this;

    // step 1: contruct read chunk for the rc1->rc2 mapping
    shared_ptr< Read_Chunk > rc1rc2_chunk_sptr;
    shared_ptr< Contig_Entry > rc1_ce_sptr;
    std::tie(rc1rc2_chunk_sptr, rc1_ce_sptr) = Read_Chunk::make_chunk_from_cigar_and_chunks(
                rc1rc2_cigar, *c1rc1_chunk_cptr, *c2rc2_chunk_cptr);

    // initialize Read_Chunk translation map, and container for new mutations
    map< Read_Chunk_CPtr, shared_ptr< Read_Chunk > > rc_map;
    Mutation_Cont extra_c1_mut_cont;

    // construct c1->rc2 mapping directly from c1->rc1 & rc1->rc2
    rc_map[c2rc2_chunk_cptr] = c1rc1_chunk_cptr->collapse_mapping(*rc1rc2_chunk_sptr, extra_c1_mut_cont);
    Read_Chunk_CPtr c1rc2_chunk_cptr = rc_map[c2rc2_chunk_cptr].get();

    // next, construct r2->c2 mapping by reversing c2->rc2
    shared_ptr< Read_Chunk > rc2c2_chunk_sptr;
    shared_ptr< Contig_Entry > rc2_ce_sptr;
    shared_ptr< Mutation_Trans_Cont > c2r2_to_r2c2_mut_trans_cont_sptr;
    std::tie(rc2c2_chunk_sptr, rc2_ce_sptr, c2r2_to_r2c2_mut_trans_cont_sptr) = c2rc2_chunk_cptr->reverse();

    // construct c1->c2 mapping
    shared_ptr< Read_Chunk > c1c2_chunk_sptr = c1rc2_chunk_cptr->collapse_mapping(*rc2c2_chunk_sptr, extra_c1_mut_cont);

    // for the remaining chunks mapped to c2, remap them to c1 through c1->c2 mapping
    for (auto it = c2_ce_cptr->get_chunk_cptr_cont().begin(); it != c2_ce_cptr->get_chunk_cptr_cont().end(); ++it)
    {
        Read_Chunk_CPtr rc_cptr = *it;
        if (rc_cptr == c2rc2_chunk_cptr)
            continue;
        rc_map[rc_cptr] = c1c2_chunk_sptr->collapse_mapping(*rc_cptr, extra_c1_mut_cont);
        rc_map[rc_cptr]->check();
    }
    // at this point, all read chunks mapped to c2 are translated
    // with pointers in rc_map and new mutations in extra_mut_cont
    ASSERT(rc_map.size() == c2_ce_cptr->get_chunk_cptr_cont().size());

    remap_chunks(rc_map, extra_c1_mut_cont);
    modify_contig_entry(c1_ce_cptr, [] (Contig_Entry& ce) {
        ce.drop_unused_mutations();
    });
    erase_contig_entry(c2_ce_cptr);

    //cerr << "after merging read chunks:\n" << *this;
    ASSERT(check(set< const Contig_Entry* >( { c1_ce_cptr })));
}

shared_ptr< vector< std::tuple< Read_Chunk_CPtr, Read_Chunk_CPtr, Cigar > > > Graph::chunker(
    const Read_Entry* re1_cptr, const Read_Entry* re2_cptr, Cigar& cigar)
{
    shared_ptr< vector< std::tuple< Read_Chunk_CPtr, Read_Chunk_CPtr, Cigar > > >
    rc_mapping_sptr(new vector< std::tuple< Read_Chunk_CPtr, Read_Chunk_CPtr, Cigar > >());

    Size_Type r1_start = cigar.get_rf_start();
    Size_Type r1_len = cigar.get_rf_len();
    Size_Type r2_start = cigar.get_qr_start();
    Size_Type r2_len = cigar.get_qr_len();
    bool r2_rc = cigar.is_reversed();
    // we repeatedly cut the read entries of either read
    // until their chunks match in the way described by the cigar string
    // keep track of read chunk mapping, and cigar transformation between them
    bool done;
    while (true)
    {
        // after every graph modification, restart at the beginning
        done = true;
        rc_mapping_sptr->clear();

        // current positions
        Size_Type r1_pos = r1_start;
        Size_Type r2_pos = (not r2_rc? r2_start : r2_start + r2_len);
        size_t op_start = 0;
        while (op_start < cigar.get_n_ops())
        {
            ASSERT(r1_start < re1_cptr->get_len());
            Read_Chunk_CPtr rc1_cptr = re1_cptr->get_chunk_with_pos(r1_pos);
            Size_Type rc1_offset = r1_pos - rc1_cptr->get_r_start();
            Size_Type rc1_remaining_len = rc1_cptr->get_r_len() - rc1_offset;
            ASSERT(rc1_offset == 0 or rc1_cptr->is_unmappable());
            ASSERT(not r2_rc or r2_pos > 0);
            Read_Chunk_CPtr rc2_cptr = (not r2_rc? re2_cptr->get_chunk_with_pos(r2_pos) : re2_cptr->get_chunk_with_pos(r2_pos - 1));
            Size_Type rc2_offset = (not r2_rc? r2_pos - rc2_cptr->get_r_start() : rc2_cptr->get_r_end() - r2_pos);
            ASSERT(rc2_offset == 0 or rc2_cptr->is_unmappable());
            Size_Type rc2_remaining_len = rc2_cptr->get_r_len() - rc2_offset;

            // invariant: we matched read 1 chunks before rc1
            // to read 2 chunks before/after rc2
            // using cigar ops before op_start
            ASSERT(r1_pos < r1_start + r1_len);
            ASSERT(r2_rc or r2_pos < r2_start + r2_len);
            ASSERT(not r2_rc or r2_pos > r2_start);

            ASSERT(rc1_cptr != NULL);
            ASSERT(rc1_cptr->get_r_start() + rc1_offset == r1_pos);
            ASSERT(rc1_remaining_len > 0);
            ASSERT(rc2_cptr != NULL);
            ASSERT(r2_rc or rc2_cptr->get_r_start() + rc2_offset == r2_pos);
            ASSERT(not r2_rc or rc2_cptr->get_r_end() - rc2_offset == r2_pos);
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
            ASSERT(op_end == cigar.get_n_ops()
                   or (cigar.get_rf_sub_len(op_start, op_end) == rc1_remaining_len
                       and cigar.get_qr_sub_len(op_start, op_end) == rc2_remaining_len)
                   or cigar.get_rf_sub_len(op_start, op_end) > rc1_remaining_len
                   or cigar.get_qr_sub_len(op_start, op_end) > rc2_remaining_len
                   or (cigar.get_rf_sub_len(op_start, op_end) == rc1_remaining_len
                       and cigar.get_qr_sub_len(op_start, op_end) < rc2_remaining_len
                       and not cigar.is_insertion(op_end))
                   or (cigar.get_rf_sub_len(op_start, op_end) < rc1_remaining_len
                       and cigar.get_qr_sub_len(op_start, op_end) == rc2_remaining_len
                       and not cigar.is_deletion(op_end))
                   or (cigar.get_rf_sub_len(op_start, op_end) == rc1_remaining_len
                       and cigar.get_qr_sub_len(op_start, op_end) == rc2_remaining_len));
            // we previously cut the chunks at the ends of the mapping
            // cuts don't succeed only if those chunks are unmappable
            ASSERT(op_end < cigar.get_n_ops()
                   or cigar.get_rf_sub_len(op_start, op_end) >= rc1_remaining_len
                   or rc1_cptr->is_unmappable());
            ASSERT(op_end < cigar.get_n_ops()
                   or cigar.get_qr_sub_len(op_start, op_end) >= rc2_remaining_len
                   or rc2_cptr->is_unmappable());
            if (cigar.get_rf_sub_len(op_start, op_end) > rc1_remaining_len)
            {
                // the first inequality trivially holds with <=
                // but it can be shown it holds in fact with <
                ASSERT(cigar.get_rf_offset(op_end - 1) < rc1_cptr->get_r_end() and rc1_cptr->get_r_end() < cigar.get_rf_offset(op_end));
            }
            if (cigar.get_qr_sub_len(op_start, op_end) > rc2_remaining_len)
            {
                // as with rf, the inequalities involving op_end-1 tivially hold with <=
                // but it can be shown they hold with <
                ASSERT(r2_rc or (cigar.get_qr_offset(op_end - 1) < rc2_cptr->get_r_end() and rc2_cptr->get_r_end() < cigar.get_qr_offset(op_end)));
                ASSERT(not r2_rc or (cigar.get_qr_offset(op_end) < rc2_cptr->get_r_start() and rc2_cptr->get_r_start() < cigar.get_qr_offset(op_end - 1)));
            }

            // check if either chunk ended during the last cigar op
            if (cigar.get_rf_sub_len(op_start, op_end) > rc1_remaining_len
                    or cigar.get_qr_sub_len(op_start, op_end) > rc2_remaining_len)
            {
                // find out which of the 2 ended earlier, cut the cigar op at that position
                Size_Type r1_break_len = 0;
                if (cigar.get_rf_sub_len(op_start, op_end) > rc1_remaining_len and not cigar.is_insertion(op_end - 1))
                {
                    r1_break_len = cigar.get_rf_op_prefix_len(op_end - 1, rc1_cptr->get_r_end());
                    ASSERT(0 < r1_break_len and r1_break_len < cigar.get_rf_op_len(op_end - 1));
                }

                Size_Type r2_break_len = 0;
                if (cigar.get_qr_sub_len(op_start, op_end) > rc2_remaining_len and not cigar.is_deletion(op_end -1))
                {
                    r2_break_len = cigar.get_qr_op_prefix_len(op_end - 1, (not r2_rc? rc2_cptr->get_r_end() : rc2_cptr->get_r_start()));
                    ASSERT(0 < r2_break_len and r2_break_len < cigar.get_qr_op_len(op_end - 1));
                }

                ASSERT(r1_break_len > 0 or r2_break_len > 0);
                cigar.cut_op(op_end - 1, (r1_break_len == 0? r2_break_len : (r2_break_len == 0? r1_break_len : min(r1_break_len, r2_break_len))));
            }

            ASSERT(cigar.get_rf_sub_len(op_start, op_end) <= rc1_remaining_len
                   and cigar.get_qr_sub_len(op_start, op_end) <= rc2_remaining_len);
            // now we are sure the op ends on at least one of the read chunk boundaries
            // unless these are the last chunks and they are both unmappable
            if (op_end == cigar.get_n_ops() and rc1_cptr->is_unmappable() and rc2_cptr->is_unmappable())
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
                    r1_pos = rc1_cptr->get_r_end();
                    op_start = op_end;
                    continue;
                }
                else if (rc2_cptr->is_unmappable())
                {
                    // rc2 cannot be cut
                    if (not rc1_cptr->is_unmappable())
                    {
                        // rc1 is mapped to an unmappable chunk; we unmap it
                        unmap_chunk(rc1_cptr);
                        done = false;
                        break;
                    }
                    else // both rc1 and rc2 are unmappable
                    {
                        // we advance both
                        r1_pos = rc1_cptr->get_r_end();
                        r2_pos = (not r2_rc? r2_pos + cigar.get_qr_sub_len(op_start, op_end) : r2_pos - cigar.get_qr_sub_len(op_start, op_end));
                        op_start = op_end;
                        continue;
                    }
                }
                else
                {
                    cut_read_entry(re2_cptr, cigar.get_qr_offset(op_end), false);
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
                    r2_pos = (not r2_rc? rc2_cptr->get_r_end() : rc2_cptr->get_r_start());
                    op_start = op_end;
                    continue;
                }
                else if (rc1_cptr->is_unmappable())
                {
                    // rc1 cannot be cut
                    if (not rc2_cptr->is_unmappable())
                    {
                        // rc2 is mapped to an unmappable chunk; we unmap it
                        unmap_chunk(rc2_cptr);
                        done = false;
                        break;
                    }
                    else // both rc1 and rc2 are unmappable
                    {
                        // we advance both
                        r1_pos += cigar.get_rf_sub_len(op_start, op_end);
                        r2_pos = (not r2_rc? rc2_cptr->get_r_end() : rc2_cptr->get_r_start());
                        op_start = op_end;
                        continue;
                    }
                }
                else
                {
                    cut_read_entry(re1_cptr, cigar.get_rf_offset(op_end), false);
                    done = false;
                    break;
                }
            }
            // reached when both rc1 and rc2 end at the current cigar op
            ASSERT(cigar.get_rf_sub_len(op_start, op_end) == rc1_remaining_len
                   and cigar.get_qr_sub_len(op_start, op_end) == rc2_remaining_len);
            if (rc1_cptr->is_unmappable() and not rc2_cptr->is_unmappable())
            {
                unmap_chunk(rc2_cptr);
                done = false;
                break;
            }
            if (rc2_cptr->is_unmappable() and not rc1_cptr->is_unmappable())
            {
                unmap_chunk(rc1_cptr);
                done = false;
                break;
            }
            ASSERT(rc1_cptr->is_unmappable() == rc2_cptr->is_unmappable());
            if (not rc1_cptr->is_unmappable())
            {
                // add read chunk mapping
                rc_mapping_sptr->push_back(std::make_tuple(rc1_cptr, rc2_cptr, cigar.substring(op_start, op_end)));
            }
            // advance both chunks and cigar
            r1_pos = rc1_cptr->get_r_end();
            r2_pos = (not r2_rc? rc2_cptr->get_r_end() : rc2_cptr->get_r_start());
            op_start = op_end;
        }
        if (done) break;
    }
    return rc_mapping_sptr;
}
*/

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

    /*
    shared_ptr< vector< std::tuple< Read_Chunk_CPtr, Read_Chunk_CPtr, Cigar > > >
    rc_mapping_sptr = chunker(re1_cptr, re2_cptr, cigar);

    // reached when we have a complete rc map
    for (size_t i = 0; i < rc_mapping_sptr->size(); ++i)
    {
        Read_Chunk_CPtr rc1_cptr;
        Read_Chunk_CPtr rc2_cptr;
        Cigar rc1rc2_cigar;
        std::tie(rc1_cptr, rc2_cptr, rc1rc2_cigar) = (*rc_mapping_sptr)[i];
        merge_read_chunks(rc1_cptr, rc2_cptr, rc1rc2_cigar);
    }

    // check for unmappable chunks in the contigs recently merged
    scan_read_for_unmappable_chunks(re1_cptr, r1_start, r1_start + r1_len);

    ASSERT(check(set< const Read_Entry* >( { re1_cptr, re2_cptr })));

    if (global::merge_contigs_at_each_step)
    {
        //cerr << "before merging:\n" << *this;
        merge_read_contigs(re1_cptr);
        //cerr << "after merging:\n" << *this;
        ASSERT(check(set< const Read_Entry* >( { re1_cptr, re2_cptr })));
    }
    */
}

/*
void Graph::reverse_contig(const Contig_Entry* ce_cptr)
{
    // construct read chunk modifier that applies reverse()
    auto rc_reverse_mod = [&] (Read_Chunk_CPtr rc_cptr) {
        modify_read_chunk(rc_cptr, [] (Read_Chunk& rc) {
            rc.reverse();
        });
    };

    // reverse the contig entry sequence and mutations (in place)
    modify_contig_entry(ce_cptr, [&rc_reverse_mod] (Contig_Entry& ce) {
        ce.reverse(rc_reverse_mod);
    });

    // reverse the mapping chunks
    / *
    auto rc_modifier = [] (Read_Chunk& rc) { rc.reverse(); };
    for (auto rc_cptr_it = ce_cptr->get_chunk_cptr_cont().begin(); rc_cptr_it != ce_cptr->get_chunk_cptr_cont().end(); ++rc_cptr_it)
    {
        modify_read_chunk(*rc_cptr_it, rc_modifier);
    }
    * /
}

bool Graph::try_merge_contig(const Contig_Entry* ce_cptr, bool dir)
{
    auto chunks_out_cont_sptr = ce_cptr->is_mergeable(dir);
    if (not chunks_out_cont_sptr)
    {
        return false;
    }

    Read_Chunk_CPtr rc_cptr = *chunks_out_cont_sptr->begin();
    Read_Chunk_CPtr rc_next_cptr = ce_cptr->get_next_chunk(dir, rc_cptr);
    const Contig_Entry* ce_next_cptr = rc_next_cptr->get_ce_ptr();
    if (rc_cptr->is_unmappable() != rc_next_cptr->is_unmappable())
    {
        // don't merge mappable with unmappable
        return false;
    }
    bool same_orientation = (rc_cptr->get_rc() == rc_next_cptr->get_rc());
    if (not same_orientation)
    {
        // we first reverse this contig
        // reversal modifications are done in-place, so chunk vectors are still valid
        reverse_contig(ce_cptr);
        dir = not dir;
        same_orientation = true;
    }
    if (not dir)
    {
        swap(ce_cptr, ce_next_cptr);
        dir = true;
        chunks_out_cont_sptr = ce_cptr->is_mergeable_one_way(true);
        ASSERT(chunks_out_cont_sptr);
    }
    // at this point:
    // - contigs are in the same orientation
    // - ce_cptr is the start of the merge and ce_next_cptr is the tail
    // - chunks_out_cont_sptr contains chunks from ce_cptr that go into ce_next_cptr

    //cerr << "merging contigs\n" << *ce_cptr << "and\n" << *ce_next_cptr;

    Size_Type prefix_len = ce_cptr->get_len();
    // construct modifier that allows merge_forward() to rebase chunks into ce_cptr
    auto rc_rebase_mod = [&] (Read_Chunk_CPtr rc_cptr, const Mutation_Trans_Cont& mut_map) {
        modify_read_chunk(rc_cptr, [&] (Read_Chunk& rc) {
            rc.rebase(ce_cptr, mut_map, prefix_len);
        });
    };
    // merge ce_next_cptr into ce_cptr using merge_forward()
    modify_contig_entry(ce_cptr, [&] (Contig_Entry& ce) {
        ce.merge_forward(ce_next_cptr, rc_rebase_mod);
    });
    // we are done with the right contig
    erase_contig_entry(ce_next_cptr);
    // merge the read chunks that were previously crossing the contig break
    set< Read_Chunk_CPtr > erased_chunks_set;
    Mutation::add_mut_mod_type add_mut_mod = [&] (const Mutation& m) {
        Mutation_CPtr res;
        modify_contig_entry(ce_cptr, [&] (Contig_Entry& ce) {
            res = ce.add_mutation(m);
        });
        return res;
    };
    for (auto rc_cptr_it = chunks_out_cont_sptr->begin(); rc_cptr_it != chunks_out_cont_sptr->end(); ++rc_cptr_it)
    {
        Read_Chunk_CPtr rc_cptr = *rc_cptr_it;
        Read_Chunk_CPtr rc_next_cptr = ce_cptr->get_next_chunk(true, rc_cptr);
        if (not rc_cptr->get_rc())
        {
            erased_chunks_set.insert(rc_next_cptr);
            modify_read_entry(rc_cptr->get_re_ptr(), [&] (Read_Entry& re) {
                re.merge_next_chunk(rc_cptr, add_mut_mod);
            });
        }
        else
        {
            erased_chunks_set.insert(rc_cptr);
            modify_read_entry(rc_next_cptr->get_re_ptr(), [&] (Read_Entry& re) {
                re.merge_next_chunk(rc_next_cptr, add_mut_mod);
            });
        }
    }
    // remove erased chunks from contig entry
    modify_contig_entry(ce_cptr, [&] (Contig_Entry& ce) {
        ce.remove_chunks(erased_chunks_set);
    });
    // drop any unused mutations
    modify_contig_entry(ce_cptr, [&] (Contig_Entry& ce) {
        ce.drop_unused_mutations();
    });

    //cerr << "merging contigs result\n" << *ce_cptr;

    return true;
}

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

void Graph::unmap_chunk(Read_Chunk_CPtr rc_cptr)
{
    // we save read entry and offset, and recompute chunk after graph modifications
    const Read_Entry* re_cptr = rc_cptr->get_re_ptr();
    Size_Type rc_start = rc_cptr->get_r_start();
    Size_Type rc_end = rc_cptr->get_r_end();
    const Contig_Entry* ce_cptr = rc_cptr->get_ce_ptr();

    // trim contig entry to the extent of the read chunk
    cut_contig_entry(ce_cptr, rc_cptr->get_c_start(), NULL);
    rc_cptr = re_cptr->get_chunk_with_pos(rc_start);
    ASSERT(rc_cptr->get_r_end() == rc_end);
    ce_cptr = rc_cptr->get_ce_ptr();
    Mutation_CPtr last_ins_cptr = NULL;
    if (rc_cptr->get_mut_ptr_cont().size() > 0
            and rc_cptr->get_mut_ptr_cont().back()->is_ins()
            and rc_cptr->get_mut_ptr_cont().back()->get_start() == rc_cptr->get_c_end())
    {
        last_ins_cptr = rc_cptr->get_mut_ptr_cont().back();
    }
    cut_contig_entry(ce_cptr, rc_cptr->get_c_end(), last_ins_cptr);
    rc_cptr = re_cptr->get_chunk_with_pos(rc_start);
    ASSERT(rc_cptr->get_r_end() == rc_end);
    ce_cptr = rc_cptr->get_ce_ptr();

    // at this point, the chunk should span the entire contig
    ASSERT(rc_cptr->get_c_start() == 0 and rc_cptr->get_c_end() == ce_cptr->get_len());

    // chunks mapping to this contig are remapped to individual contigs, with unmappable flag set
    // we also save them in a list using Read_Entry and offsets
    list< const Read_Entry* > re_list;
    list< std::tuple< Size_Type, Size_Type > > pos_list;
    for (auto rc_cptr_it = ce_cptr->get_chunk_cptr_cont().begin(); rc_cptr_it != ce_cptr->get_chunk_cptr_cont().end(); ++rc_cptr_it)
    {
        rc_cptr = *rc_cptr_it;
        re_list.push_back(rc_cptr->get_re_ptr());
        pos_list.push_back(std::make_tuple(rc_cptr->get_r_start(), rc_cptr->get_r_end()));
        const Contig_Entry* new_ce_cptr = insert_contig_entry(Contig_Entry(rc_cptr));
        modify_read_chunk(rc_cptr, [&] (Read_Chunk& rc) {
            rc.assign_to_contig(new_ce_cptr, 0, new_ce_cptr->get_len(), false, vector<Mutation_CPtr>());
            rc.set_unmappable();
        });
    }

    // done with old Contig_Entry
    erase_contig_entry(ce_cptr);

    // for all read chunks made unmappable, try to merge them with nearby unmappable chunks
    auto re_it = re_list.begin();
    auto pos_it = pos_list.begin();
    for (; re_it != re_list.end(); ++re_it, ++pos_it)
    {
        re_cptr = *re_it;
        std::tie(rc_start, rc_end) = *pos_it;
        rc_cptr = re_cptr->get_chunk_with_pos(rc_start);
        ASSERT(rc_cptr->is_unmappable());
        ASSERT(rc_cptr->get_r_start() <= rc_start and rc_end <= rc_cptr->get_r_end());
        if (rc_cptr->get_r_start() != rc_start or rc_cptr->get_r_end() != rc_end)
        {
            // already been merged with other unmappable chunks
            continue;
        }
        extend_unmapped_chunk(re_cptr, rc_start, rc_end);
    }
    ASSERT(check(set< const Read_Entry* >(re_list.begin(), re_list.end())));
}

void Graph::extend_unmapped_chunk_dir(const Read_Entry* re_cptr, Size_Type pos, bool dir)
{
    Size_Type leftover_bp = (not dir? pos : re_cptr->get_len() - pos);
    while (leftover_bp > 0)
    {
        // calls in previous iterations might extend the unmapped region
        // we first recompute the chunks based solely on: re_cptr & unmap_start
        Read_Chunk_CPtr rc_cptr = re_cptr->get_chunk_with_pos(not dir? pos : pos - 1);
        ASSERT(rc_cptr->is_unmappable());
        ASSERT(dir or rc_cptr->get_r_start() <= pos);
        ASSERT(not dir or rc_cptr->get_r_end() >= pos);
        pos = (not dir? min(pos, rc_cptr->get_r_start()) : max(pos, rc_cptr->get_r_end()));
        leftover_bp = (not dir? pos : re_cptr->get_len() - pos);
        if (leftover_bp == 0)
        {
            break;
        }
        Read_Chunk_CPtr next_rc_cptr = re_cptr->get_chunk_with_pos(not dir? pos - 1 : pos);

        if (not next_rc_cptr->is_unmappable() and leftover_bp <= global::unmap_trigger_len)
        {
            unmap_chunk(next_rc_cptr);
            continue;
        }
        if (next_rc_cptr->is_unmappable())
        {
            pos = (not dir? next_rc_cptr->get_r_start() : next_rc_cptr->get_r_end());
            merge_unmappable_chunks((not dir? next_rc_cptr : rc_cptr),
                                    (not dir? rc_cptr : next_rc_cptr));
            continue;
        }

        // next chunk is mappable
        // check there is a minimum of mappable bp
        // if not, unmap small intermediate chunks and merge them
        ASSERT(not next_rc_cptr->is_unmappable());
        ASSERT(leftover_bp > global::unmap_trigger_len);
        // skip small mappable chunks
        list< std::tuple< Size_Type, Size_Type > > skipped_chunks;
        Size_Type skipped_len = 0;
        while (not next_rc_cptr->is_unmappable()
                and skipped_len + next_rc_cptr->get_r_len() <= global::unmap_trigger_len)
        {
            skipped_chunks.push_back(std::make_tuple(next_rc_cptr->get_r_start(), next_rc_cptr->get_r_end()));
            skipped_len += next_rc_cptr->get_r_len();
            next_rc_cptr = re_cptr->get_sibling(next_rc_cptr, dir);
            // since skipped_len < unmap_trigger_len <= leftover_bp
            ASSERT(next_rc_cptr != NULL);
        }
        ASSERT(skipped_len <= global::unmap_trigger_len);
        if (next_rc_cptr->is_unmappable())
        {
            // unmap all skipped chunks (if they are not modified in previous iterations)
            for (auto it = skipped_chunks.begin(); it != skipped_chunks.end(); ++it)
            {
                Size_Type rc_start;
                Size_Type rc_end;
                std::tie(rc_start, rc_end) = *it;
                rc_cptr = re_cptr->get_chunk_with_pos(rc_start);
                if (rc_cptr->get_r_start() != rc_start or rc_cptr->get_r_end() != rc_end or rc_cptr->is_unmappable())
                {
                    // already been merged with other unmappable chunks
                    continue;
                }
                unmap_chunk(rc_cptr);
            }
        }
        else
        {
            break;
        }
    }
    ASSERT(check(set< const Read_Entry* >( { re_cptr })));
}

void Graph::extend_unmapped_chunk(const Read_Entry* re_cptr, Size_Type rc_start, Size_Type rc_end)
{
    extend_unmapped_chunk_dir(re_cptr, rc_start, false);
    extend_unmapped_chunk_dir(re_cptr, rc_end, true);
}

void Graph::merge_unmappable_chunks(Read_Chunk_CPtr rc1_cptr, Read_Chunk_CPtr rc2_cptr)
{
    ASSERT(rc1_cptr->is_unmappable() and rc2_cptr->is_unmappable());
    ASSERT(rc1_cptr->get_re_ptr() == rc2_cptr->get_re_ptr());
    ASSERT(rc1_cptr->get_r_end() == rc2_cptr->get_r_start());
    bool success;
    success = try_merge_contig(rc1_cptr->get_ce_ptr(), true);
    ASSERT(success);
}

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
