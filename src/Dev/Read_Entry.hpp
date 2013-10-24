//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __READ_ENTRY_HPP
#define __READ_ENTRY_HPP

#include <string>
#include <iostream>

#include "MAC_forward.hpp"
#include "Read_Chunk.hpp"


namespace MAC
{
    class Read_Entry
    {
    public:
        Read_Entry() {}

        const std::string& get_name() const { return name; }

        friend void add_read(const std::string&, const Seq_Type&, Read_Entry*, Contig_Entry*);

        friend std::ostream& operator << (std::ostream& os, const Read_Entry& rhs);

    private:
        std::string name;
        Read_Chunk_Cont chunk_cont;
    };
}


#endif
