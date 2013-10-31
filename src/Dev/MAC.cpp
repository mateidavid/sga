#include "MAC.hpp"
#include "Cigar.hpp"

using namespace std;
using boost::tie;


namespace MAC
{
    void add_read(Read_Entry_Cont& re_cont, Contig_Entry_Cont& ce_cont,
                  const string* name_ptr, const Seq_Type* seq_ptr)
    {
        // first, create read entry & contig entry
        Read_Entry re(name_ptr, seq_ptr->size());
        Contig_Entry ce(seq_ptr);
        ce.add_chunk(&(*re.get_chunk_cont().begin()));

        // insert them in their containers
        Read_Entry_Cont::iterator re_it;
        bool success;
        tie(re_it, success) = re_cont.insert(re);
        assert(success);
        Contig_Entry_Cont::iterator ce_it;
        tie(ce_it, success) = ce_cont.insert(ce_cont.end(), ce);
        assert(success);

        // set the contig entry pointer inside the single read chunk
        auto rc_modifier = [&] (Read_Chunk& rc) { rc.assign_to_contig(&(*ce_it), 0, seq_ptr->size(), false, vector<const Mutation*>()); };
        auto re_modifier = [&] (Read_Entry& re) { re.modify_read_chunk(&(*re_it->get_chunk_cont().begin()), rc_modifier); };
        re_cont.modify(re_it, re_modifier);
    }

    void cut_mutation(Read_Entry_Cont& re_cont, Contig_Entry_Cont& ce_cont,
                      const Contig_Entry* ce_cptr, const Mutation* mut_cptr, Size_Type c_offset, Size_Type r_offset)
    {
        // construct contig entry modifier which:
        // - cuts given mutation
        // - inserts the remaing part in its container
        // - returns pointer to new mutation
        const Mutation* m_new_cptr;
        auto ce_modifier = [&] (Contig_Entry& ce) { m_new_cptr = ce.cut_mutation(mut_cptr, c_offset, r_offset); };

        // apply modifier
        bool success;
        success = ce_cont.modify(ce_cont.iterator_to(*ce_cptr), ce_modifier);
        assert(success);

        // construct read chunk modifier which adds new mutation if original one was present
        auto rc_modifier = [&] (Read_Chunk& rc) { rc.cond_add_mutation(mut_cptr, m_new_cptr); };

        // fetch read chunks which need to change
        vector< Read_Chunk_CPtr > v = ce_cptr->get_chunks_with_mutation(mut_cptr);

        for (size_t i = 0; i < v.size(); ++i)
        {
            // construct read entry modifier which applies read chunk modifier
            auto re_modifier = [&] (Read_Entry& re) { re.modify_read_chunk(v[i], rc_modifier); };

            // apply re & rc modifier, thus adding mutation 2nd part
            success = re_cont.modify(re_cont.iterator_to(*(v[i]->get_re_ptr())), re_modifier);
            assert(success);
        }
    }

    void cut_contig_entry(Read_Entry_Cont& re_cont, Contig_Entry_Cont& ce_cont,
                          const Contig_Entry* ce_cptr, Size_Type c_pos,
                          Read_Chunk_CPtr rc_cptr, Size_Type r_pos)
    {
        (void)rc_cptr;
        (void)r_pos;

        // first, split any mutations that span c_pos
        vector< const Mutation* > v = ce_cptr->get_mutations_spanning_pos(c_pos);
        for (size_t i = 0; i < v.size(); ++i)
        {
            cut_mutation(re_cont, ce_cont, ce_cptr, v[i], c_pos - v[i]->get_start(), 0);
        }

        //
        
    }

    void cut_read_chunk(Read_Entry_Cont& re_cont, Contig_Entry_Cont& ce_cont, Read_Chunk_CPtr rc_cptr, Size_Type r_brk)
    {
        const Mutation* mut_cptr;
        Size_Type c_offset;
        Size_Type r_offset;
        boost::tie(mut_cptr, c_offset, r_offset) = rc_cptr->find_mutation_to_cut(r_brk);
        if (mut_cptr != NULL)
        {
            cut_mutation(re_cont, ce_cont, rc_cptr->get_ce_ptr(), mut_cptr, c_offset, r_offset);
        }
        // now we are certain the breakpoint no longer falls inside a mutation (insertion/mnp)

        Size_Type r_pos;
        Size_Type c_pos;
        size_t i;
        boost::tie(r_pos, c_pos, i) = rc_cptr->get_read_split_data(r_brk);
        assert(r_pos == r_brk);

        // cut contig at given c_offset & remeber to cut current chunk at r_offset
        cut_contig_entry(re_cont, ce_cont,
                         rc_cptr->get_ce_ptr(), c_pos,
                         rc_cptr, r_brk);
    }

    void cut_read_entry(Read_Entry_Cont& re_cont, Contig_Entry_Cont& ce_cont, const Read_Entry* re_cptr, Size_Type r_pos)
    {
        // find read chunk that must be cut
        Read_Chunk_CPtr rc_cptr = re_cptr->get_chunk_with_pos(r_pos);

        // if cut is at the start of the chunk, there is nothing to do
        if (rc_cptr == NULL or rc_cptr->get_r_start() == r_pos)
            return;
        assert(rc_cptr->get_r_start() < r_pos and r_pos < rc_cptr->get_r_end());

        cut_read_chunk(re_cont, ce_cont, rc_cptr, r_pos);
    }

    void add_overlap(Read_Entry_Cont& re_cont, Contig_Entry_Cont& ce_cont,
                     const string& r1_name, const string& r2_name,
                     Size_Type r1_start, Size_Type r1_len,
                     Size_Type r2_start, Size_Type r2_len, bool r2_rc,
                     const string& cigar_string)
    {
        (void)ce_cont;

        Read_Entry_Cont::iterator re1_it = re_cont.find(r1_name);
        assert(re1_it != re_cont.end());
        Read_Entry_Cont::iterator re2_it = re_cont.find(r2_name);
        assert(re2_it != re_cont.end());

        // construct 2 read chunk objects corresponding to the alignment
        Cigar cigar(cigar_string, r2_rc);
        assert(r1_len == cigar.get_rf_len());
        assert(r2_len == cigar.get_qr_len());
        vector< pair< Read_Chunk, Mutation_Cont > > v = Read_Chunk::make_chunks_from_cigar(r1_start, r2_start, cigar);

        // cut r1 & r2 at the ends of the match region

    }
}
