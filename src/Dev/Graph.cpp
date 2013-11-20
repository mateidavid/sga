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

    void Graph::add_read(const string* name_ptr, Seq_Type* seq_ptr)
    {
        // first, create read entry and place it in container
        Read_Entry re(name_ptr, seq_ptr->size());
        cerr << "re:" << indent::inc << indent::nl << re << indent::dec << indent::nl;
        const Read_Entry* re_cptr = insert_read_entry(re);

        // create contig entry and place it in container
        Contig_Entry ce(seq_ptr);
        cerr << "ce:" << indent::inc << indent::nl << ce << indent::dec << indent::nl;
        ce.add_chunk(&(*re_cptr->get_chunk_cont().begin()));
        cerr << "ce with chunk:" << indent::inc << indent::nl << ce << indent::dec << indent::nl;
        const Contig_Entry* ce_cptr = insert_contig_entry(ce);

        // fix initial rc: assing it to contig entry
        modify_read_chunk(&(*re_cptr->get_chunk_cont().begin()),
                          [&] (Read_Chunk& rc) { rc.assign_to_contig(ce_cptr, 0, seq_ptr->size(), false, vector< const Mutation* >()); });
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
        if (c_brk == 0 or c_brk == ce_cptr->get_len())
        {
            return;
        }

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
        // drop those mutations from first part
        modify_contig_entry(ce_cptr, [&] (Contig_Entry& ce) { ce.drop_mutations(mut_cptr_map); });

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

        // finally, drop the second part of the base sequence from the first part
        modify_contig_entry(ce_cptr, [&] (Contig_Entry& ce) { ce.drop_base_seq(c_brk); });
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
                 and (*rc_cptr->get_mut_ptr_cont().rbegin())->is_ins())
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

    void merge_read_chunks(Read_Chunk_CPtr rc1_cptr, Read_Chunk_CPtr rc2_cptr, Cigar& cigar)
    {
        (void)rc1_cptr;
        (void)rc2_cptr;
        (void)cigar;

        // do not do anything if the chunks are already mapped to the same contig
        // NOTE: with this, we are ignoring alternate mappings
        if (rc1_cptr->get_ce_ptr() == rc2_cptr->get_ce_ptr())
            return;

        // step 1: contruct read chunk for the rc1->rc2 mapping
        shared_ptr< Read_Chunk > rc1rc2_chunk_sptr;
        shared_ptr< Mutation_Cont > rc1rc2_mut_cont_sptr;
        std::tie(rc1rc2_chunk_sptr, rc1rc2_mut_cont_sptr) = Read_Chunk::make_chunk_from_cigar(cigar);

        // translate rc1 mutations observed by rc2 into c1 mutations
        shared_ptr< Mutation_Trans_Cont > rc1rc2_to_c1rc2_mut_trans_cont_sptr;
        shared_ptr< Mutation_Cont > c1rc2_mut_cont_sptr;
        std::tie(rc1rc2_to_c1rc2_mut_trans_cont_sptr, c1rc2_mut_cont_sptr) = rc1_cptr->lift_read_mutations(*rc1rc2_mut_cont_sptr);

        shared_ptr< vector< std::tuple< Mutation_CPtr, Size_Type, Size_Type, bool > > > c1rc2_mut_list_sptr;
        c1rc2_mut_list_sptr = rc1_cptr->get_mutations_under_mapping(
            rc1rc2_chunk_sptr->get_mut_ptr_cont(), *rc1rc2_to_c1rc2_mut_trans_cont_sptr, c1rc2_mut_cont_sptr);

        //TODO
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
        clog << indent::nl << "adding overlap:" << indent::inc
             << indent::nl << "re1: " << r1_seq.substr(r1_start, r1_len)
             << indent::nl << "re2: " << (not cigar.is_reversed()? r2_seq.substr(r2_start, r2_len) : reverseComplement(r2_seq.substr(r2_start, r2_len)))
             << indent::nl << "initial cigar:" << indent::inc << cigar << indent::dec;
        cigar.disambiguate(r1_seq.substr(r1_start, r1_len), r2_seq.substr(r2_start, r2_len));
        clog << indent::nl << "disambiguated cigar:" << indent::inc << cigar << indent::dec << indent::dec;

        // cut r1 & r2 at the ends of the match region
        cut_read_entry(re1_cptr, r1_start, true);
        cut_read_entry(re1_cptr, r1_start + r1_len, true);
        cut_read_entry(re2_cptr, r2_start, true);
        cut_read_entry(re2_cptr, r2_start + r2_len, true);

        // initially, we repeatedly cut the read entries of either read
        // until their chunks match in the way described by the cigar string
        // keep track of read chunk mapping, and cigar transformation between them
        vector< tuple< Read_Chunk_CPtr, Read_Chunk_CPtr, Cigar > > rc_mapping;
        bool done;
        while (true)
        {
            done = true;
            rc_mapping.clear();

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
                while (cigar.get_rf_sub_len(op_start, op_end) < rc1_cptr->get_r_len()
                       and cigar.get_qr_sub_len(op_start, op_end) < rc2_cptr->get_r_len())
                {
                    ++op_end;
                }
                assert(cigar.get_rf_offset(op_end - 1) < rc1_cptr->get_r_end() and rc1_cptr->get_r_end() <= cigar.get_rf_offset(op_end));
                assert(r2_rc or (cigar.get_qr_offset(op_end - 1) < rc2_cptr->get_r_end() and rc2_cptr->get_r_end() <= cigar.get_qr_offset(op_end)));
                assert(not r2_rc or (cigar.get_qr_offset(op_end) <= rc2_cptr->get_r_start() and rc2_cptr->get_r_start() < cigar.get_qr_offset(op_end - 1)));

                // check if either chunk ended during the last cigar op
                if (cigar.get_rf_sub_len(op_start, op_end) > rc1_cptr->get_r_len()
                    or cigar.get_qr_sub_len(op_start, op_end) > rc2_cptr->get_r_len())
                {
                    // find out which of the 2 ended earlier, cut the cigar op at that position
                    Size_Type r1_break_len = 0;
                    if (cigar.get_rf_sub_len(op_start, op_end) > rc1_cptr->get_r_len())
                    {
                        r1_break_len = cigar.get_rf_op_prefix_len(op_end, rc1_cptr->get_r_end());
                    }

                    Size_Type r2_break_len = 0;
                    if (cigar.get_qr_sub_len(op_start, op_end) > rc2_cptr->get_r_len())
                    {
                        r2_break_len = cigar.get_qr_op_prefix_len(op_end, (not r2_rc? rc2_cptr->get_r_end() : rc2_cptr->get_r_start()));
                    }

                    assert(r1_break_len > 0 or r2_break_len > 0);
                    cigar.cut_op(op_end, (r1_break_len == 0? r2_break_len : (r2_break_len == 0? r1_break_len : min(r1_break_len, r2_break_len))));
                    ++op_end;
                }
                // now we are sure the op ends on at least one of the read chunk boundaries
                assert(cigar.get_rf_sub_len(op_start, op_end) <= rc1_cptr->get_r_len()
                       and cigar.get_qr_sub_len(op_start, op_end) <= rc2_cptr->get_r_len());
                assert(cigar.get_rf_sub_len(op_start, op_end) == rc1_cptr->get_r_len()
                       or cigar.get_qr_sub_len(op_start, op_end) == rc2_cptr->get_r_len());
                // if it doesn't end on both, we might need to cut the other
                if (cigar.get_qr_sub_len(op_start, op_end) < rc2_cptr->get_r_len())
                {
                    // rc1 ends; check if there are insertions we can append to the match
                    while (op_end < cigar.get_n_ops()
                           and cigar.is_insertion(op_end)
                           and cigar.get_qr_sub_len(op_start, op_end) < rc2_cptr->get_r_len())
                    {
                        ++op_end;
                    }
                    if (op_end < cigar.get_n_ops()
                        and cigar.is_insertion(op_end)
                        and cigar.get_qr_sub_len(op_start, op_end) > rc2_cptr->get_r_len())
                    {
                        // we went past rc2 end on insertions only
                        cigar.cut_op(op_end, cigar.get_qr_op_prefix_len(op_end, (not r2_rc? rc2_cptr->get_r_end() : rc2_cptr->get_r_start())));
                        ++op_end;
                    }
                    // if at this point we don't touch rc2 end, we have to cut rc2 or discard rc1
                    if (cigar.get_qr_sub_len(op_start, op_end) < rc2_cptr->get_r_len())
                    {
                        if (cigar.get_qr_sub_len(op_start, op_end) == 0)
                        {
                            // no progress on rc2: rc1 is mapped entirely to a deletion
                            // move on to next chunk on read 1
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
                }
                if (cigar.get_rf_sub_len(op_start, op_end) < rc1_cptr->get_r_len())
                {
                    // rc2 ends; check if there are deletions we can append to the match
                    while (op_end < cigar.get_n_ops()
                           and cigar.is_deletion(op_end)
                           and cigar.get_rf_sub_len(op_start, op_end) < rc1_cptr->get_r_len())
                    {
                        ++op_end;
                    }
                    if (op_end < cigar.get_n_ops()
                        and cigar.is_deletion(op_end)
                        and cigar.get_rf_sub_len(op_start, op_end) > rc1_cptr->get_r_len())
                    {
                        // we went past rc1 end on deletions only
                        cigar.cut_op(op_end, cigar.get_rf_op_prefix_len(op_end, rc1_cptr->get_r_end()));
                        ++op_end;
                    }
                    // if at this point we don't touch rc1 end, we have to cut rc1 or discard rc2
                    if (cigar.get_rf_sub_len(op_start, op_end) < rc1_cptr->get_r_len())
                    {
                        if (cigar.get_rf_sub_len(op_start, op_end) == 0)
                        {
                            // no progress on rc1: rc2 mapped entirely to an insertion
                            // move on to next chunk on read 2
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
                }
                // reached when both rc1 and rc2 end at the current cigar op
                assert(cigar.get_rf_sub_len(op_start, op_end) == rc1_cptr->get_r_len()
                       and cigar.get_qr_sub_len(op_start, op_end) == rc2_cptr->get_r_len());
                // add read chunk mapping
                rc_mapping.push_back(make_tuple(rc1_cptr, rc2_cptr, cigar.substring(op_start, op_end)));
                // advance both chunks and cigar
                r1_pos = rc1_cptr->get_r_end();
                r2_pos = (not r2_rc? rc2_cptr->get_r_end() : rc2_cptr->get_r_start());
                op_start = op_end;
            }
            if (done) break;
        }
        // reached when we have a complete rc map
        for (size_t i = 0; i < rc_mapping.size(); ++i)
        {
            merge_read_chunks(get<0>(rc_mapping[i]), get<1>(rc_mapping[i]), get<2>(rc_mapping[i]));
        }
    }

    void Graph::check() const
    {
        size_t chunks_count_1 = 0;
        size_t chunks_count_2 = 0;
        // check read entry objects
        for (auto re_it = _re_cont.begin(); re_it != _re_cont.end(); ++re_it)
        {
            re_it->check();
            chunks_count_1 += re_it->get_chunk_cont().size();
        }
        // check contig entry objects
        for (auto ce_it = _ce_cont.begin(); ce_it != _ce_cont.end(); ++ce_it)
        {
            ce_it->check();
            chunks_count_2 += ce_it->get_chunk_cptr_cont().size();
        }
        assert(chunks_count_1 == chunks_count_2);
    }

    ostream& operator << (ostream& os, const Graph& g)
    {
        os << "Read_Entry_Cont=" << indent::inc;
        print_seq(os, g._re_cont, indent::nl, indent::nl);
        os << indent::dec << indent::nl << "Contig_Entry_Cont=" << indent::inc;
        print_seq(os, g._ce_cont, indent::nl, indent::nl);
        os << indent::dec;
        return os;
    }
}
