//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Read_Chunk.hpp"
#include "Mutation.hpp"


namespace MAC
{
    std::ostream& operator << (std::ostream& os, const Read_Chunk& rhs)
    {
        os << "(r_start=" << rhs.r_start << ",r_len=" << rhs.r_len
        << ",c_start=" << rhs.c_start << ",c_len=" << rhs.c_len << ",rc=" << (int)rhs.rc << ",mut_list=(";
        std::for_each(rhs.mut_dlist.begin(), rhs.mut_dlist.end(), [&] (const Mutation_Desc& m) {
            if (m != *rhs.mut_dlist.begin()) os << ",";
                      os << *m;
        });
        os << "))";
        return os;
    }
}
