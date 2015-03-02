//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Read_Chunk_Cont.hpp"
#include "Read_Entry.hpp"


namespace MAC
{

map< Read_Chunk_CBPtr, Read_Chunk_CBPtr >
Read_Chunk_CE_Cont::splice(Read_Chunk_CE_Cont& other_cont,
                           Size_Type c_brk, Mutation_CBPtr mut_left_cbptr, bool strict)
{
    ASSERT(not mut_left_cbptr
           or (mut_left_cbptr->rf_start() == c_brk and mut_left_cbptr->is_ins()));
    map< Read_Chunk_CBPtr, Read_Chunk_CBPtr > res;
    Read_Chunk_CE_Cont lhs_cont;
    Read_Chunk_CE_Cont rhs_cont;
    while (not other_cont.empty())
    {
        Read_Chunk_BPtr left_rc_bptr;
        Read_Chunk_BPtr right_rc_bptr;
        Read_Chunk_BPtr rc_bptr = &*other_cont.begin();
        other_cont.erase(rc_bptr);
        tie(left_rc_bptr, right_rc_bptr) = Read_Chunk::split(rc_bptr, c_brk, mut_left_cbptr, strict);
        if (left_rc_bptr and right_rc_bptr)
        {
            ASSERT(left_rc_bptr == rc_bptr);
            res[left_rc_bptr] = right_rc_bptr;
        }
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
    return res;
}

auto Read_Chunk_CE_Cont::erase_from_re_cont() const -> RC_Next_Map
{
    RC_Next_Map res;
    for (auto rc_cbptr : *this | referenced)
    {
        Read_Entry_BPtr re_bptr = rc_cbptr->re_bptr().unconst();
        if (re_bptr)
        {
            auto it = ++(re_bptr->chunk_cont().iterator_to(*rc_cbptr));
            res[rc_cbptr] = it;
            re_bptr->chunk_cont().erase(rc_cbptr);
        }
    }
    return res;
}

void Read_Chunk_CE_Cont::insert_into_re_cont(const RC_Next_Map& rc_next_map) const
{
    for (auto rc_cbptr : *this | referenced)
    {
        Read_Entry_BPtr re_bptr = rc_cbptr->re_bptr().unconst();
        auto it = rc_next_map.at(rc_cbptr);
        ASSERT(it == re_bptr->chunk_cont().end() or rc_cbptr->get_r_end() <= it->get_r_start());
        re_bptr->chunk_cont().insert_before(it, rc_cbptr.unconst());
    }
}

} // namespace MAC
