//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Read_Chunk.hpp"

#include "Mutation.hpp"
#include "print_seq.hpp"

using namespace std;


namespace MAC
{
    ostream& operator << (ostream& os, const Read_Chunk& rhs)
    {
        os << "(r_start=" << rhs.r_start << ",r_len=" << rhs.r_len
        << ",c_start=" << rhs.c_start << ",c_len=" << rhs.c_len << ",rc=" << (int)rhs.rc << ",mut_list=(";
        print_ptr_seq(os, rhs.mut_ptr_cont);
        os << "))";
        return os;
    }
}
