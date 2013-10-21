//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Contig_Entry.hpp"
#include "Mutation.hpp"
#include "Read_Chunk.hpp"
#include "Read_Entry.hpp"


namespace MAC
{
    std::ostream& operator << (std::ostream& os, const Contig_Entry& rhs)
    {
        os << "(seq=" << rhs.seq << ",\nmut_list=\n";
        std::for_each(rhs.mut_list.begin(), rhs.mut_list.end(), [&] (const Mutation& m) {
            os << "  " << m << "\n";
        });
        os << "chunk_list=\n";
        std::for_each(rhs.chunk_dlist.begin(), rhs.chunk_dlist.end(), [&] (const Read_Chunk_Desc& rc_d) {
            os << "  (r_name=" << rc_d->get_read_entry_desc()->get_name() << "," << *rc_d << ")\n";
        });
        os << ")\n";
        return os;
    }
}
