#include "Graph.hpp"
#include "Cigar.hpp"
#include "indent.hpp"
#include "print_seq.hpp"

using namespace std;
using boost::tie;


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
        boost::tie(mut_cptr, c_offset, r_offset) = rc_cptr->get_mutation_to_cut(r_brk, false);
        if (mut_cptr != NULL)
        {
            cut_mutation(rc_cptr->get_ce_ptr(), mut_cptr, c_offset, r_offset);
        }
        // now we are certain the breakpoint no longer falls inside a mutation (insertion/mnp)

        Size_Type c_pos;
        Size_Type r_pos;
        size_t i;
        boost::tie(c_pos, r_pos, i) = rc_cptr->get_cut_data(r_brk, false);
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

    void Graph::add_overlap(const string& r1_name, const string& r2_name,
                            Size_Type r1_start, Size_Type r1_len,
                            Size_Type r2_start, Size_Type r2_len, bool r2_rc,
                            const string& cigar_string)
    {
        const Read_Entry* re1_cptr = get_read_entry(r1_name);
        assert(re1_cptr != NULL);
        const Read_Entry* re2_cptr = get_read_entry(r2_name);
        assert(re2_cptr != NULL);

        // construct 2 read chunk objects corresponding to the alignment
        Cigar cigar(cigar_string, r2_rc);
        assert(r1_len == cigar.get_rf_len());
        assert(r2_len == cigar.get_qr_len());
        vector< pair< Read_Chunk, Mutation_Cont > > v = Read_Chunk::make_chunks_from_cigar(r1_start, r2_start, cigar);

        // cut r1 & r2 at the ends of the match region
        cut_read_entry(re1_cptr, r1_start, true);
        cut_read_entry(re1_cptr, r1_start + r1_len, true);
        cut_read_entry(re2_cptr, r2_start, true);
        cut_read_entry(re2_cptr, r2_start + r2_len, true);
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
