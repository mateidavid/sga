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
        os << "(r_start=" << rhs.get_r_start() << ",r_len=" << rhs.get_r_len()
        << ",c_start=" << rhs.get_c_start() << ",c_len=" << rhs.get_c_len() << ",rc=" << (int)rhs.get_rc() << ",mut_list=(";
        print_ptr_seq(os, rhs._mut_ptr_cont);
        os << "))";
        return os;
    }
}
