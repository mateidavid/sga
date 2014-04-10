//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Read_Chunk_Cont.hpp"

using namespace std;


namespace MAC
{

Read_Chunk_CE_Cont Read_Chunk_CE_Cont::split(Size_Type c_brk, Mutation_CBPtr mut_left_cbptr)
{
    ASSERT(mut_left_cbptr == nullptr or (mut_left_cbptr->get_start() == c_brk and mut_left_cbptr->is_ins()));
    Read_Chunk_CE_Cont lhs_cont;
    Read_Chunk_CE_Cont rhs_cont;
    
    while (size() > 0)
    {
        Read_Chunk_BPtr left_rc_bptr;
        Read_Chunk_BPtr right_rc_bptr;
        Read_Chunk_BPtr rc_bptr = &*begin();
        erase(rc_bptr);
        std::tie(left_rc_bptr, right_rc_bptr) = Read_Chunk::split(rc_bptr, c_brk, mut_left_cbptr);
        lhs_cont.insert(left_rc_bptr);
        rhs_cont.insert(right_rc_bptr);
    }

    *this = std::move(lhs_cont);
    return rhs_cont;
}

} // namespace MAC
