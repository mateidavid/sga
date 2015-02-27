//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Read_Chunk_Cont.hpp"
#include "Read_Entry.hpp"


namespace MAC
{

void Read_Chunk_CE_Cont::splice(Read_Chunk_CE_Cont& other_cont,
                                Size_Type c_brk, Mutation_CBPtr mut_left_cbptr, bool strict)
{
    ASSERT(not mut_left_cbptr
           or (mut_left_cbptr->rf_start() == c_brk and mut_left_cbptr->is_ins()));
    Read_Chunk_CE_Cont lhs_cont;
    Read_Chunk_CE_Cont rhs_cont;
    while (not other_cont.empty())
    {
        Read_Chunk_BPtr left_rc_bptr;
        Read_Chunk_BPtr right_rc_bptr;
        Read_Chunk_BPtr rc_bptr = &*other_cont.begin();
        other_cont.erase(rc_bptr);
        tie(left_rc_bptr, right_rc_bptr) = Read_Chunk::split(rc_bptr, c_brk, mut_left_cbptr, strict);
        if (left_rc_bptr)
        {
            lhs_cont.insert(left_rc_bptr);
        }
        if (right_rc_bptr)
        {
            rhs_cont.insert(right_rc_bptr);
        }
    }
    other_cont = move(lhs_cont);
    *this = move(rhs_cont);
}

void Read_Chunk_CE_Cont::erase_from_re_cont() const
{
    for (auto rc_cbptr : *this | referenced)
    {
        Read_Entry_BPtr re_bptr = rc_cbptr->re_bptr().unconst();
        if (re_bptr)
        {
            re_bptr->chunk_cont().erase(rc_cbptr);
        }
    }
}

void Read_Chunk_CE_Cont::insert_into_re_cont() const
{
    for (auto rc_cbptr : *this | referenced)
    {
        Read_Entry_BPtr re_bptr = rc_cbptr->re_bptr().unconst();
        re_bptr->chunk_cont().insert(rc_cbptr.unconst());
    }
}

} // namespace MAC
