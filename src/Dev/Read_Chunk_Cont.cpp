//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Read_Chunk_Cont.hpp"
#include "Read_Entry.hpp"


namespace MAC
{

auto Read_Chunk_CE_Cont::splice(
    Read_Chunk_CE_Cont& other_cont,
    Size_Type c_brk, Mutation_CBPtr mut_left_cbptr, bool strict) -> RC_Map
{
    ASSERT(not mut_left_cbptr
           or (mut_left_cbptr->rf_start() == c_brk and mut_left_cbptr->is_ins()));
    RC_Map res;
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

auto Read_Chunk_CE_Cont::erase_from_re_cont() const -> RC_Map
{
    RC_Map res;
    // first save iterators to the next chunk
    for (auto rc_cbptr : *this | referenced)
    {
        Read_Entry_BPtr re_bptr = rc_cbptr->re_bptr().unconst();
        //ASSERT(re_bptr);
        if (re_bptr)
        {
            auto it = ++(re_bptr->chunk_cont().iterator_to(*rc_cbptr));
            res[rc_cbptr] = (it != re_bptr->chunk_cont().end()? &*it : nullptr);
        }
    }
    // then remove them
    for (auto rc_cbptr : *this | referenced)
    {
        Read_Entry_BPtr re_bptr = rc_cbptr->re_bptr().unconst();
        if (re_bptr)
        {
            re_bptr->chunk_cont().erase(rc_cbptr);
        }
    }
    return res;
}

void Read_Chunk_CE_Cont::insert_into_re_cont(RC_Map& rc_map)
{
    while (not rc_map.empty())
    {
        // find a chunk to insert whose next chunk is not in the map (i.e., it is already inserted)
        auto it = rc_map.begin();
        while (it != rc_map.end() and rc_map.count(it->second) > 0) ++it;
        ASSERT(it != rc_map.end());
        Read_Chunk_BPtr rc_bptr = it->first.unconst();
        Read_Chunk_CBPtr rc_next_cbptr = it->second;
        ASSERT(not rc_next_cbptr or rc_bptr->get_r_end() <= rc_next_cbptr->get_r_start());
        Read_Entry_BPtr re_bptr = rc_bptr->re_bptr().unconst();
        auto rc_next_it = (rc_next_cbptr?
                           re_bptr->chunk_cont().iterator_to(*rc_next_cbptr)
                           : re_bptr->chunk_cont().end());
        re_bptr->chunk_cont().insert_before(rc_next_it, rc_bptr);
        rc_map.erase(it);
    }
}

} // namespace MAC
