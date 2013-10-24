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


namespace MAC
{
    std::ostream& operator << (std::ostream& os, const Contig_Entry& rhs)
    {
        os << "(seq=" << rhs.get_seq() << ",\nmut_list=\n  ";
        print_seq(os, rhs.mut_cont, "\n  ");
        os << "\nchunk_list=\n  ";
        print_ptr_seq(os, rhs.chunk_cptr_cont, "\n  ");
        os << "\n)\n";
        return os;
    }
}
