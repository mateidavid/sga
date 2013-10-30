//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Contig_Entry.hpp"
#include "Mutation.hpp"
#include "Read_Chunk.hpp"
#include "Read_Entry.hpp"
#include "print_seq.hpp"

using namespace std;


namespace MAC
{
    const Mutation* Contig_Entry::cut_mutation(const Mutation* mut_cptr, Size_Type c_offset, Size_Type r_offset)
    {
        // find mutation in container
        Mutation_Cont::iterator it;
        it = _mut_cont.iterator_to(*mut_cptr);
        assert(it != _mut_cont.end());

        // cut existing mutation, saving remaining part
        Mutation m_new;
        auto mut_modifier = [&] (Mutation& m) { m_new = m.cut(c_offset, r_offset); };
        bool success;
        success = _mut_cont.modify(it, mut_modifier);
        assert(success);

        // insert remaining part
        Mutation_Cont::iterator it_new;
        boost::tie(it_new, success) = _mut_cont.insert(m_new);

        return &(*it_new);
    }

    vector< Read_Chunk_CPtr > Contig_Entry::get_chunks_with_mutation(const Mutation* mut_cptr) const
    {
        vector< Read_Chunk_CPtr > res;
        for (auto it = _chunk_cptr_cont.begin(); it != _chunk_cptr_cont.end(); ++it)
        {
            if ((*it)->have_mutation(mut_cptr))
                res.push_back(*it);
        }
        return res;
    }

    vector< const Mutation* > Contig_Entry::get_mutations_spanning_pos(Size_Type c_pos) const
    {
        vector< const Mutation* > res;
        for (auto it = _mut_cont.begin(); it != _mut_cont.end() and it->get_start() < c_pos; ++it)
        {
            if (c_pos < it->get_end())
                res.push_back(&(*it));
        }
        return res;
    }

    std::ostream& operator << (std::ostream& os, const Contig_Entry& rhs)
    {
        os << "(seq=" << rhs.get_seq() << ",\nmut_list=\n  ";
        print_seq(os, rhs._mut_cont, "\n  ");
        os << "\nchunk_list=\n  ";
        print_ptr_seq(os, rhs._chunk_cptr_cont, "\n  ");
        os << "\n)\n";
        return os;
    }
}
