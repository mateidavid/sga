//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Read_Entry.hpp"
#include "Read_Chunk.hpp"


namespace MAC
{
    std::ostream& operator << (std::ostream& os, const Read_Entry& rhs)
    {
        os << "(name=" << rhs.name << ",\nchunk_list=\n";
        std::for_each(rhs.chunk_list.begin(), rhs.chunk_list.end(), [&] (const Read_Chunk& c) {
            os << "  " << c << "\n";
        });
        os << ")\n";
        return os;
    }
}
