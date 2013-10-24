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
    ostream& operator << (ostream& os, const Read_Entry& rhs)
    {
        os << "(name=" << rhs.get_name() << ",\nchunk_list=\n  ";
        print_seq(os, rhs.chunk_cont, "\n  ");
        os << "\n)\n";
        return os;
    }
}
