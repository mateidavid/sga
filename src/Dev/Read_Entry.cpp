//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Read_Entry.hpp"

#include "print_seq.hpp"

using namespace std;


namespace MAC
{
    const Read_Chunk* Read_Entry::get_chunk_with_pos(Size_Type r_pos) const
    {
        if (r_pos >= get_len())
            return NULL;
        Read_Chunk_Cont::const_iterator it = _chunk_cont.lower_bound(r_pos);
        if (it != _chunk_cont.end() and it->get_r_start() == r_pos)
            return &(*it);
        else
            return &(*(--it));
    }

    ostream& operator << (ostream& os, const Read_Entry& rhs)
    {
        os << "(name=" << rhs.get_name() << ",\nchunk_list=\n  ";
        print_seq(os, rhs._chunk_cont, "\n  ");
        os << "\n)\n";
        return os;
    }
}
