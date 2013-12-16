#include "Graph.hpp"
#include "Cigar.hpp"
#include "indent.hpp"
#include "print_seq.hpp"
#include "../Util/Util.h"

using namespace std;


namespace MAC
{
    const Read_Entry* Graph::insert_read_entry(const Read_Entry& re)
    {
        Read_Entry_Cont::iterator re_it;
        bool success;
        tie(re_it, success) = _re_cont.insert(re);
        assert(success);
        return &(*re_it);
    }

    const Contig_Entry* Graph::insert_contig_entry(const Contig_Entry& ce)
    {
        Contig_Entry_Cont::iterator ce_it;
        bool success;
        tie(ce_it, success) = _ce_cont.insert(_ce_cont.end(), ce);
        assert(success);
        return &(*ce_it);
    }

    void Graph::erase_contig_entry(const Contig_Entry* ce_cptr)
    {
        for (auto it = ce_cptr->get_chunk_cptr_cont().begin(); it != ce_cptr->get_chunk_cptr_cont().end(); ++it)
        {
            // the read chunks no longer point here
            assert((*it)->get_ce_ptr() != ce_cptr);
        }
        _ce_cont.erase(_ce_cont.iterator_to(*ce_cptr));
    }

    void Graph::add_read(const string* name_ptr, Seq_Type* seq_ptr)
    {
        // first, create read entry and place it in container
        Read_Entry re(name_ptr, seq_ptr->size());
        //cerr << indent::tab << "re:\n" << indent::inc << re << indent::dec;
        const Read_Entry* re_cptr = insert_read_entry(re);

        // create contig entry and place it in container
        Contig_Entry ce(seq_ptr);
        //cerr << indent::tab << "ce:\n" << indent::inc << ce << indent::dec;
        ce.add_chunk(&(*re_cptr->get_chunk_cont().begin()));
        //cerr << indent::tab << "ce with chunk:\n" << indent::inc << ce << indent::dec;
        const Contig_Entry* ce_cptr = insert_contig_entry(ce);

        // fix initial rc: assing it to contig entry
        modify_read_chunk(&(*re_cptr->get_chunk_cont().begin()),
                          [&] (Read_Chunk& rc) { rc.assign_to_contig(ce_cptr, 0, seq_ptr->size(), false, vector< const Mutation* >()); });

        assert(check(set< const Read_Entry* >({ re_cptr })));
    }

    void Graph::cut_mutation(const Contig_Entry* ce_cptr, const Mutation* mut_cptr, Size_Type c_offset, Size_Type r_offset)
    {
        // construct & apply contig entry modifier which:
        // - cuts given mutation
        // - inserts the remaing part in its container
        // - returns pointer to new mutation
        const Mutation* m_new_cptr;
        modify_contig_entry(ce_cptr, [&] (Contig_Entry& ce) { m_new_cptr = ce.cut_mutation(mut_cptr, c_offset, r_offset); });

        // construct read chunk modifier which adds new mutation if original one was present
        auto rc_modifier = [&] (Read_Chunk& rc) { rc.cond_add_mutation(mut_cptr, m_new_cptr); };

        // fetch read chunks which need to change, and apply modifier to each of them
        vector< Read_Chunk_CPtr > v = ce_cptr->get_chunks_with_mutation(mut_cptr);
        for (size_t i = 0; i < v.size(); ++i)
        {
            modify_read_chunk(v[i], rc_modifier);
        }
    }

    void Graph::cut_contig_entry(const Contig_Entry* ce_cptr, Size_Type c_brk, const Mutation* mut_left_cptr)
    {
        if ((c_brk == 0 and mut_left_cptr == NULL)
            or (c_brk == ce_cptr->get_len() and (ce_cptr->get_mut_cont().size() == 0 or not ce_cptr->get_mut_cont().rbegin()->is_ins())))
        {
            return;
        }

        //cerr << "before cutting contig entry: " << (void*)ce_cptr << " at " << c_brk << " mut_left_cptr=" << (void*)mut_left_cptr << '\n' << *this;

        // first, split any mutations that span c_pos
        vector< const Mutation* > v = ce_cptr->get_mutations_spanning_pos(c_brk);
        for (size_t i = 0; i < v.size(); ++i)
        {
            assert(v[i]->get_start() < c_brk and c_brk < v[i]->get_end());
            assert(v[i] != mut_left_cptr);
            cut_mutation(ce_cptr, v[i], c_brk - v[i]->get_start(), 0);
        }

        // create new contig entry object; set string; add it to container
        const Contig_Entry* ce_new_cptr;
        ce_new_cptr = insert_contig_entry(Contig_Entry(new string(ce_cptr->get_seq().substr(c_brk))));

        // aquire second-half mutations from original contig entry object into the new one; compute mutation pointer map
        map< const Mutation*, const Mutation* > mut_cptr_map;
        auto ce_modifier = [&] (Contig_Entry& ce) { mut_cptr_map = ce.acquire_second_half_mutations(ce_cptr, c_brk, mut_left_cptr); };
        modify_contig_entry(ce_new_cptr, ce_modifier);

        // save the list of read chunk objects that must be modified
        Read_Chunk_CPtr_Cont chunk_cptr_cont = ce_cptr->get_chunk_cptr_cont();

        for (auto rc_cptr_it = chunk_cptr_cont.begin(); rc_cptr_it != chunk_cptr_cont.end(); ++rc_cptr_it)
        {
            Read_Chunk_CPtr rc_cptr = (*rc_cptr_it);
            bool move_to_rhs;
            shared_ptr< Read_Chunk > rc_new_sptr;

            // construct&apply read chunk modifier to implement split
            auto rc_modifier = [&] (Read_Chunk& rc) { tie(move_to_rhs, rc_new_sptr) = rc.apply_contig_split(c_brk, mut_cptr_map, ce_new_cptr); };
            modify_read_chunk(rc_cptr, rc_modifier);

            Read_Chunk_CPtr rc_rhs_cptr = NULL;
            // if a new read chunk is created, insert it in read entry object
            if (rc_new_sptr)
            {
                modify_read_entry(rc_cptr->get_re_ptr(), [&] (Read_Entry& re) { rc_rhs_cptr = re.add_read_chunk(rc_new_sptr.get()); });
            }
            // if the rc must move to the rhs ce, remove it from lhs ce
            else if (move_to_rhs)
            {
                modify_contig_entry(ce_cptr, [&] (Contig_Entry& ce) { ce.remove_chunk(rc_cptr); });
                rc_rhs_cptr = rc_cptr;
            }

            // if a rc must be added to rhs ce, do it now
            if (rc_rhs_cptr != NULL)
            {
                modify_contig_entry(ce_new_cptr, [&] (Contig_Entry& ce) { ce.add_chunk(rc_rhs_cptr); });
            }
        }

        // drop mapped mutations from first part
        modify_contig_entry(ce_cptr, [&] (Contig_Entry& ce) { ce.drop_mutations(mut_cptr_map); });
        // drop the second part of the base sequence from the first part
        modify_contig_entry(ce_cptr, [&] (Contig_Entry& ce) { ce.drop_base_seq(c_brk); });
        // fix contig length for chunks left on LHS
        for (auto rc_cptr_it = ce_cptr->get_chunk_cptr_cont().begin(); rc_cptr_it != ce_cptr->get_chunk_cptr_cont().end(); ++rc_cptr_it)
        {
            modify_read_chunk(*rc_cptr_it, [&] (Read_Chunk& rc) { rc.c_len() = ce_cptr->get_len(); });
        }

        // drop unused mutations from either contig
        modify_contig_entry(ce_cptr, [] (Contig_Entry& ce) { ce.drop_unused_mutations(); });
        modify_contig_entry(ce_new_cptr, [] (Contig_Entry& ce) { ce.drop_unused_mutations(); });

        //cerr << "after cutting contig entry: " << (void*)ce_cptr << '\n' << *this;
        assert(check(set< const Contig_Entry* >({ ce_cptr, ce_new_cptr })));
    }

    void Graph::cut_read_chunk(Read_Chunk_CPtr rc_cptr, Size_Type r_brk)
    {
        const Mutation* mut_cptr;
        Size_Type c_offset;
        Size_Type r_offset;
        tie(mut_cptr, c_offset, r_offset) = rc_cptr->get_mutation_to_cut(r_brk, false);
        if (mut_cptr != NULL)
        {
            cut_mutation(rc_cptr->get_ce_ptr(), mut_cptr, c_offset, r_offset);
        }
        // now we are certain the breakpoint no longer falls inside a mutation (insertion/mnp)

        Size_Type c_pos;
        Size_Type r_pos;
        size_t i;
        tie(c_pos, r_pos, i) = rc_cptr->get_cut_data(r_brk, false);
        assert(r_pos == r_brk);

        // check if an insertion at c_pos has to remain on the left of the cut
        const Mutation* mut_left_cptr = NULL;
        if (i > 0
            and rc_cptr->get_mut_ptr_cont()[i-1]->get_start() == c_pos
            and rc_cptr->get_mut_ptr_cont()[i-1]->is_ins())
        {
            mut_left_cptr = rc_cptr->get_mut_ptr_cont()[i-1];
        }

        // cut contig at given c_offset
        // all insertions at c_pos go to the right except for at most one, passed as argument
        cut_contig_entry(rc_cptr->get_ce_ptr(), c_pos, mut_left_cptr);
    }

    void Graph::cut_read_entry(const Read_Entry* re_cptr, Size_Type r_brk, bool force)
    {
        // if cut is at the very start or end, and the force flag is on, cut the underlying contig directly
        if (r_brk == 0 or r_brk == re_cptr->get_len())
        {
             if (not force)
             {
                 return;
             }

             // compute contig cut coordinates: ce_cptr & c_brk
             Read_Chunk_CPtr rc_cptr = (r_brk == 0?
                                        &*re_cptr->get_chunk_cont().begin()
                                        : &*re_cptr->get_chunk_cont().rbegin());
             Size_Type c_brk = ((r_brk == 0) == (not rc_cptr->get_rc())?
                                rc_cptr->get_c_start()
                                : rc_cptr->get_c_end());
             // if read chunk ends in insertion, that must remain on the left of the cut
             const Mutation* mut_left_cptr = NULL;
             if (c_brk == rc_cptr->get_c_end()
                 and rc_cptr->get_mut_ptr_cont().size() > 0
                 and (*rc_cptr->get_mut_ptr_cont().rbegin())->is_ins()
                 and (*rc_cptr->get_mut_ptr_cont().rbegin())->get_start() == c_brk)
             {
                 mut_left_cptr = *(rc_cptr->get_mut_ptr_cont().rbegin());
             }
             cut_contig_entry(rc_cptr->get_ce_ptr(), c_brk, mut_left_cptr);
        }
        else
        {
            // find read chunk that must be cut
            Read_Chunk_CPtr rc_cptr = re_cptr->get_chunk_with_pos(r_brk);

            // if cut is at the start of the chunk, there is nothing to do
            if (rc_cptr == NULL or rc_cptr->get_r_start() == r_brk)
                return;
            assert(rc_cptr->get_r_start() < r_brk and r_brk < rc_cptr->get_r_end());

            cut_read_chunk(rc_cptr, r_brk);
        }
    }

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
            assert(c2rc_cptr->get_ce_ptr() == ce2_cptr);
            assert(c1rc_ptr->get_ce_ptr() == ce1_cptr);
            assert(c2rc_cptr->get_re_ptr() == c1rc_ptr->get_re_ptr());
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
                assert(extra_mut_map[extra_mut_cptr] != NULL);
                c1rc_ptr->mut_ptr_cont()[i] = extra_mut_map[extra_mut_cptr];
            }
            // done moving mutations
            // now replace the read chunk in its read entry, and update pointers
            modify_read_chunk(c2rc_cptr, [&] (Read_Chunk& rc) { rc = *c1rc_ptr; });
            // now c2rc_cptr has the c1rc object
            // add chunk to ce1
            modify_contig_entry(ce1_cptr, [&] (Contig_Entry& ce) { ce.add_chunk(c2rc_cptr); });
        }
    }

    void Graph::merge_read_chunks(Read_Chunk_CPtr c1rc1_chunk_cptr, Read_Chunk_CPtr c2rc2_chunk_cptr, Cigar& rc1rc2_cigar)
    {
        assert(rc1rc2_cigar.check(c1rc1_chunk_cptr->get_seq(), c2rc2_chunk_cptr->get_seq()));
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
            if (*it == c2rc2_chunk_cptr)
            {
                continue;
            }
            rc_map[*it] = c1c2_chunk_sptr->collapse_mapping(**it, extra_c1_mut_cont);
        }
        // at this point, all read chunks mapped to c2 are translated
        // with pointers in rc_map and new mutations in extra_mut_cont
        assert(rc_map.size() == c2_ce_cptr->get_chunk_cptr_cont().size());

        remap_chunks(rc_map, extra_c1_mut_cont);
        modify_contig_entry(c1_ce_cptr, [] (Contig_Entry& ce) { ce.drop_unused_mutations(); });
        erase_contig_entry(c2_ce_cptr);

        //cerr << "after merging read chunks:\n" << *this;
        assert(check(set< const Contig_Entry* >({ c1_ce_cptr })));
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
            done = true;
            rc_mapping_sptr->clear();

            // current positions
            Size_Type r1_pos = r1_start;
            Size_Type r2_pos = (not r2_rc? r2_start : r2_start + r2_len);
            size_t op_start = 0;
            while (op_start < cigar.get_n_ops())
            {
                assert(r1_start < re1_cptr->get_len());
                Read_Chunk_CPtr rc1_cptr = re1_cptr->get_chunk_with_pos(r1_pos);
                assert(not r2_rc or r2_pos > 0);
                Read_Chunk_CPtr rc2_cptr = (not r2_rc? re2_cptr->get_chunk_with_pos(r2_pos) : re2_cptr->get_chunk_with_pos(r2_pos - 1));

                // invariant: we matched read 1 chunks before rc1
                // to read 2 chunks before/after rc2
                // using cigar ops before op_start
                assert(r1_pos < r1_start + r1_len);
                assert(r2_rc or r2_pos < r2_start + r2_len);
                assert(not r2_rc or r2_pos > r2_start);

                assert(rc1_cptr != NULL);
                assert(rc1_cptr->get_r_start() == r1_pos);
                assert(rc1_cptr->get_r_len() > 0);
                assert(rc2_cptr != NULL);
                assert(r2_rc or rc2_cptr->get_r_start() == r2_pos);
                assert(not r2_rc or rc2_cptr->get_r_end() == r2_pos);
                assert(rc2_cptr->get_r_len() > 0);

                assert(r1_pos == cigar.get_rf_offset(op_start));
                assert(r2_pos == cigar.get_qr_offset(op_start));

                // advance past cigar ops until either chunk ends
                size_t op_end = op_start + 1;
                while (op_end < cigar.get_n_ops()
                       and (cigar.get_rf_sub_len(op_start, op_end) < rc1_cptr->get_r_len()
                            or (cigar.get_rf_sub_len(op_start, op_end) == rc1_cptr->get_r_len() and cigar.is_insertion(op_end)))
                       and (cigar.get_qr_sub_len(op_start, op_end) < rc2_cptr->get_r_len()
                            or (cigar.get_qr_sub_len(op_start, op_end) == rc2_cptr->get_r_len() and cigar.is_deletion(op_end))))
                {
                    ++op_end;
                }
                // stop conditions: can be derived by inspecting loop
                assert(op_end == cigar.get_n_ops()
                       or (cigar.get_rf_sub_len(op_start, op_end) == rc1_cptr->get_r_len()
                           and cigar.get_qr_sub_len(op_start, op_end) == rc2_cptr->get_r_len())
                       or cigar.get_rf_sub_len(op_start, op_end) > rc1_cptr->get_r_len()
                       or cigar.get_qr_sub_len(op_start, op_end) > rc2_cptr->get_r_len()
                       or (cigar.get_rf_sub_len(op_start, op_end) == rc1_cptr->get_r_len()
                           and cigar.get_qr_sub_len(op_start, op_end) < rc2_cptr->get_r_len()
                           and not cigar.is_insertion(op_end))
                       or (cigar.get_rf_sub_len(op_start, op_end) < rc1_cptr->get_r_len()
                           and cigar.get_qr_sub_len(op_start, op_end) == rc2_cptr->get_r_len()
                           and not cigar.is_deletion(op_end)));
                if (cigar.get_rf_sub_len(op_start, op_end) > rc1_cptr->get_r_len())
                {
                    // the first inequality trivially holds with <=
                    // but it can be shown it holds in fact with <
                    assert(cigar.get_rf_offset(op_end - 1) < rc1_cptr->get_r_end() and rc1_cptr->get_r_end() < cigar.get_rf_offset(op_end));
                }
                if (cigar.get_qr_sub_len(op_start, op_end) > rc2_cptr->get_r_len())
                {
                    // as with rf, the inequalities involving op_end-1 tivially hold with <=
                    // but it can be shown they hold with <
                    assert(r2_rc or (cigar.get_qr_offset(op_end - 1) < rc2_cptr->get_r_end() and rc2_cptr->get_r_end() < cigar.get_qr_offset(op_end)));
                    assert(not r2_rc or (cigar.get_qr_offset(op_end) < rc2_cptr->get_r_start() and rc2_cptr->get_r_start() < cigar.get_qr_offset(op_end - 1)));
                }

                // check if either chunk ended during the last cigar op
                if (cigar.get_rf_sub_len(op_start, op_end) > rc1_cptr->get_r_len()
                    or cigar.get_qr_sub_len(op_start, op_end) > rc2_cptr->get_r_len())
                {
                    // find out which of the 2 ended earlier, cut the cigar op at that position
                    Size_Type r1_break_len = 0;
                    if (cigar.get_rf_sub_len(op_start, op_end) > rc1_cptr->get_r_len() and not cigar.is_insertion(op_end - 1))
                    {
                        r1_break_len = cigar.get_rf_op_prefix_len(op_end - 1, rc1_cptr->get_r_end());
                        assert(0 < r1_break_len and r1_break_len < cigar.get_rf_op_len(op_end - 1));
                    }

                    Size_Type r2_break_len = 0;
                    if (cigar.get_qr_sub_len(op_start, op_end) > rc2_cptr->get_r_len() and not cigar.is_deletion(op_end -1))
                    {
                        r2_break_len = cigar.get_qr_op_prefix_len(op_end - 1, (not r2_rc? rc2_cptr->get_r_end() : rc2_cptr->get_r_start()));
                        assert(0 < r2_break_len and r2_break_len < cigar.get_qr_op_len(op_end - 1));
                    }

                    assert(r1_break_len > 0 or r2_break_len > 0);
                    cigar.cut_op(op_end - 1, (r1_break_len == 0? r2_break_len : (r2_break_len == 0? r1_break_len : min(r1_break_len, r2_break_len))));
                }
                // now we are sure the op ends on at least one of the read chunk boundaries
                assert(cigar.get_rf_sub_len(op_start, op_end) <= rc1_cptr->get_r_len()
                       and cigar.get_qr_sub_len(op_start, op_end) <= rc2_cptr->get_r_len());
                assert(cigar.get_rf_sub_len(op_start, op_end) == rc1_cptr->get_r_len()
                       or cigar.get_qr_sub_len(op_start, op_end) == rc2_cptr->get_r_len());
                // if it doesn't end on both, we might need to cut the other
                if (cigar.get_qr_sub_len(op_start, op_end) < rc2_cptr->get_r_len())
                {
                    // it follows from main stop condition that the next op is not an insertion
                    assert(op_end < cigar.get_n_ops() and not cigar.is_insertion(op_end));
                    if (cigar.get_qr_sub_len(op_start, op_end) == 0)
                    {
                        // no progress on rc2: rc1 is mapped entirely to a deletion
                        // move on to next chunk on read 1
                        assert(cigar.get_rf_sub_len(op_start, op_end) == rc1_cptr->get_r_len()
                               and cigar.get_rf_sub_len(op_start, op_end) > 0);
                        r1_pos = rc1_cptr->get_r_end();
                        op_start = op_end;
                        continue;
                    }
                    else
                    {
                        cut_read_entry(re2_cptr, cigar.get_qr_offset(op_end), false);
                        done = false;
                        break;
                    }
                }
                if (cigar.get_rf_sub_len(op_start, op_end) < rc1_cptr->get_r_len())
                {
                    // it follows from main stop condition that the next op is not a deletion
                    assert(op_end < cigar.get_n_ops() and not cigar.is_deletion(op_end));
                    if (cigar.get_rf_sub_len(op_start, op_end) == 0)
                    {
                        // no progress on rc1: rc2 mapped entirely to an insertion
                        // move on to next chunk on read 2
                        assert(cigar.get_qr_sub_len(op_start, op_end) == rc2_cptr->get_r_len()
                               and cigar.get_qr_sub_len(op_start, op_end) > 0);
                        r2_pos = (not r2_rc? rc2_cptr->get_r_end() : rc2_cptr->get_r_start());
                        op_start = op_end;
                        continue;
                    }
                    else
                    {
                        cut_read_entry(re1_cptr, cigar.get_rf_offset(op_end), false);
                        done = false;
                        break;
                    }
                }
                // reached when both rc1 and rc2 end at the current cigar op
                assert(cigar.get_rf_sub_len(op_start, op_end) == rc1_cptr->get_r_len()
                       and cigar.get_qr_sub_len(op_start, op_end) == rc2_cptr->get_r_len());
                // add read chunk mapping
                rc_mapping_sptr->push_back(std::make_tuple(rc1_cptr, rc2_cptr, cigar.substring(op_start, op_end)));
                // advance both chunks and cigar
                r1_pos = rc1_cptr->get_r_end();
                r2_pos = (not r2_rc? rc2_cptr->get_r_end() : rc2_cptr->get_r_start());
                op_start = op_end;
            }
            if (done) break;
        }
        return rc_mapping_sptr;
    }

    void Graph::add_overlap(const string& r1_name, const string& r2_name,
                            Size_Type r1_start, Size_Type r1_len,
                            Size_Type r2_start, Size_Type r2_len, bool r2_rc,
                            const string& cigar_string)
    {
        // construct cigar object
        Cigar cigar(cigar_string, r2_rc, r1_start, r2_start);
        assert(r1_len == cigar.get_rf_len());
        assert(r2_len == cigar.get_qr_len());

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
            return;

        const Read_Entry* re1_cptr = get_read_entry(r1_name);
        assert(re1_cptr != NULL);
        const Read_Entry* re2_cptr = get_read_entry(r2_name);
        assert(re2_cptr != NULL);

        string r1_seq = re1_cptr->get_seq();
        string r2_seq = re2_cptr->get_seq();
        /*
        cerr << indent::tab << "adding overlap:" << indent::inc
             << indent::nl << "re1: " << r1_seq.substr(r1_start, r1_len)
             << indent::nl << "re2: " << (not cigar.is_reversed()? r2_seq.substr(r2_start, r2_len) : reverseComplement(r2_seq.substr(r2_start, r2_len)))
             << indent::nl << "initial cigar:\n" << indent::inc << cigar << indent::dec;
        */
        cigar.disambiguate(r1_seq.substr(r1_start, r1_len), r2_seq.substr(r2_start, r2_len));
        /*
        cerr << indent::tab << "disambiguated cigar:\n" << indent::inc << cigar << indent::dec << indent::dec;
        */

        // cut r1 & r2 at the ends of the match region
        cut_read_entry(re1_cptr, r1_start, true);
        cut_read_entry(re1_cptr, r1_start + r1_len, true);
        cut_read_entry(re2_cptr, r2_start, true);
        cut_read_entry(re2_cptr, r2_start + r2_len, true);

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

        // check if we can merge any adjacent contigs
        //merge_contigs(re1_cptr);

        assert(check(set< const Read_Entry* >({ re1_cptr, re2_cptr })));
    }

    void Graph::reverse_contig(const Contig_Entry* ce_cptr)
    {
        // construct read chunk modifier that applies reverse()
        auto rc_reverse_mod = [&] (Read_Chunk_CPtr rc_cptr) {
            modify_read_chunk(rc_cptr, [] (Read_Chunk& rc) { rc.reverse(); });
        };

        // reverse the contig entry sequence and mutations (in place)
        modify_contig_entry(ce_cptr, [&rc_reverse_mod] (Contig_Entry& ce) { ce.reverse(rc_reverse_mod); });

        // reverse the mapping chunks
        /*
        auto rc_modifier = [] (Read_Chunk& rc) { rc.reverse(); };
        for (auto rc_cptr_it = ce_cptr->get_chunk_cptr_cont().begin(); rc_cptr_it != ce_cptr->get_chunk_cptr_cont().end(); ++rc_cptr_it)
        {
            modify_read_chunk(*rc_cptr_it, rc_modifier);
        }
        */
    }

    bool Graph::try_merge_contig(const Contig_Entry* ce_cptr, bool dir)
    {
        auto chunks_out_cont_sptr = ce_cptr->is_mergeable(dir);
        if (not chunks_out_cont_sptr)
            return false;
        Read_Chunk_CPtr rc_cptr = *chunks_out_cont_sptr->begin();
        Read_Chunk_CPtr rc_next_cptr = ce_cptr->get_next_chunk(dir, rc_cptr);
        const Contig_Entry* ce_next_cptr = rc_next_cptr->get_ce_ptr();
        bool same_orientation = (rc_cptr->get_rc() == rc_next_cptr->get_rc());
        if (not same_orientation)
        {
            // we first reverse this contig
            reverse_contig(ce_cptr);
            // reversal modifications are done in-place, so chunk vectors are still valid
            same_orientation = true;
        }
        if (not dir)
        {
            swap(ce_cptr, ce_next_cptr);
            dir = true;
            chunks_out_cont_sptr = ce_next_cptr->is_mergeable_one_way(true);
            assert(chunks_out_cont_sptr);
        }
        // at this point:
        // - contigs are in the same orientation
        // - ce_cptr is the start of the merge and ce_next_cptr is the tail
        // - chunks_out_cont_sptr contains chunks from ce_cptr that go into ce_next_cptr

        Size_Type prefix_len = ce_cptr->get_len();
        // merge ce_next_cptr into ce_cptr using merge_forward()
        shared_ptr< Mutation_Trans_Cont > mut_map_sptr;
        modify_contig_entry(ce_cptr, [&ce_next_cptr,&mut_map_sptr] (Contig_Entry& ce) { mut_map_sptr = ce.merge_forward(ce_next_cptr); });
        // apply rebase modifier to all chunks from the right contig
        for (auto rc_cptr_it = ce_next_cptr->get_chunk_cptr_cont().begin(); rc_cptr_it != ce_next_cptr->get_chunk_cptr_cont().end(); ++rc_cptr_it)
        {
            rc_cptr = *rc_cptr_it;
            modify_read_chunk(rc_cptr, [&] (Read_Chunk& rc) { rc.rebase(ce_cptr, prefix_len, *mut_map_sptr); });
        }
        // we are done with the right contig
        erase_contig_entry(ce_next_cptr);
        // finally, merge the read chunks that were previously crossing the contig break
        for (auto rc_cptr_it = chunks_out_cont_sptr->begin(); rc_cptr_it != chunks_out_cont_sptr->end(); ++rc_cptr_it)
        {
            Read_Chunk_CPtr rc_cptr = *rc_cptr_it;
            if (not rc_cptr->get_rc())
            {
                modify_read_entry(rc_cptr->get_re_ptr(), [&rc_cptr] (Read_Entry& re) { re.merge_next_chunk(rc_cptr); });
            }
            else
            {
                Read_Chunk_CPtr rc_next_cptr = ce_cptr->get_next_chunk(true, rc_cptr);
                modify_read_entry(rc_next_cptr->get_re_ptr(), [&rc_next_cptr] (Read_Entry& re) { re.merge_next_chunk(rc_next_cptr); });
            }
        }
        return true;
    }

    void Graph::merge_read_contigs(const Read_Entry* re_cptr)
    {
        bool done;
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
        while (not done);
    }

    bool Graph::check_all() const
    {
        size_t chunks_count_1 = 0;
        size_t chunks_count_2 = 0;
        // check read entry objects
        for (auto re_it = _re_cont.begin(); re_it != _re_cont.end(); ++re_it)
        {
            assert(re_it->check());
            chunks_count_1 += re_it->get_chunk_cont().size();
        }
        // check contig entry objects
        for (auto ce_it = _ce_cont.begin(); ce_it != _ce_cont.end(); ++ce_it)
        {
            assert(ce_it->check());
            chunks_count_2 += ce_it->get_chunk_cptr_cont().size();
        }
        assert(chunks_count_1 == chunks_count_2);
        return true;
    }

    bool Graph::check(const set< const Read_Entry* >& re_set, const set< const Contig_Entry* >& ce_set) const
    {
#ifdef NO_GRAPH_CHECKS
        return true;
#endif
        // compute contig entries referenced by read entries
        set< const Contig_Entry* > ce_extra_set;
        for (auto re_cptr_it = re_set.begin(); re_cptr_it != re_set.end(); ++re_cptr_it)
        {
            const Read_Entry* re_cptr = *re_cptr_it;
            for (auto rc_it = re_cptr->get_chunk_cont().begin(); rc_it != re_cptr->get_chunk_cont().end(); ++rc_it)
            {
                ce_extra_set.insert(rc_it->get_ce_ptr());
            }
        }

        // compute read entries referenced by contig entries
        set< const Read_Entry* > re_extra_set;
        for (auto ce_cptr_it = ce_set.begin(); ce_cptr_it != ce_set.end(); ++ce_cptr_it)
        {
            const Contig_Entry* ce_cptr = *ce_cptr_it;
            for (auto rc_cptr_it = ce_cptr->get_chunk_cptr_cont().begin(); rc_cptr_it != ce_cptr->get_chunk_cptr_cont().end(); ++rc_cptr_it)
            {
                Read_Chunk_CPtr rc_cptr = *rc_cptr_it;
                re_extra_set.insert(rc_cptr->get_re_ptr());
            }
        }

        // check read entry objects
        for (auto re_cptr_it = re_set.begin(); re_cptr_it != re_set.end(); ++re_cptr_it)
        {
            const Read_Entry* re_cptr = *re_cptr_it;
            assert(re_cptr->check());
        }
        for (auto re_cptr_it = re_extra_set.begin(); re_cptr_it != re_extra_set.end(); ++re_cptr_it)
        {
            const Read_Entry* re_cptr = *re_cptr_it;
            if (re_set.count(re_cptr) > 0)
                continue;
            assert(re_cptr->check());
        }

        // check contig entry objects
        for (auto ce_cptr_it = ce_set.begin(); ce_cptr_it != ce_set.end(); ++ce_cptr_it)
        {
            const Contig_Entry* ce_cptr = *ce_cptr_it;
            assert(ce_cptr->check());
        }
        for (auto ce_cptr_it = ce_extra_set.begin(); ce_cptr_it != ce_extra_set.end(); ++ce_cptr_it)
        {
            const Contig_Entry* ce_cptr = *ce_cptr_it;
            if (ce_set.count(ce_cptr) > 0)
                continue;
            assert(ce_cptr->check());
        }

        return true;
    }

    void Graph::dump_detailed_counts(ostream& os) const
    {
        // First read stats
        os << "RE\tname\tlen\tnum.chunks\tchunk.lens\n";
        for (auto re_it = _re_cont.begin(); re_it != _re_cont.end(); ++re_it)
        {
            const Read_Entry* re_cptr = &*re_it;
            os << "RE\t" << re_cptr->get_name() << '\t' << re_cptr->get_len() << '\t' << re_cptr->get_chunk_cont().size() << '\t';
            for (auto rc_it = re_cptr->get_chunk_cont().begin(); rc_it != re_cptr->get_chunk_cont().end(); ++rc_it)
            {
                Read_Chunk_CPtr rc_cptr = &*rc_it;
                if (rc_it != re_cptr->get_chunk_cont().begin())
                    os << ',';
                os << rc_cptr->get_r_len();
            }
            os << '\n';
        }
        // next, contig stats
        os << "CE\tlen\tnum.chunks\tbp.chunks\tnum.mut\tnum.mut.chunks"
           << "\tnum.snp\tnum.ins\tnum.del\tnum.mnp\tbp.mut"
           << "\tcovg.left\tdeg.left\tcovg.right\tdeg.right\n";
        for (auto ce_it = _ce_cont.begin(); ce_it != _ce_cont.end(); ++ce_it)
        {
            const Contig_Entry* ce_cptr = &*ce_it;
            os << "CE\t" << ce_cptr->get_len() << '\t' << ce_cptr->get_chunk_cptr_cont().size() << '\t';
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
            os << cnt_0 << '\t' << uniq_0 << '\t' << cnt_1 << '\t' << uniq_1;
            os << '\n';
        }
    }

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
}
