//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __CONTIG_ENTRY_HPP
#define __CONTIG_ENTRY_HPP

#include "MAC_forward.hpp"


namespace MAC
{
    class Contig_Entry
    {
    public:
        Contig_Entry() {}

        friend void add_read(const std::string&, const Seq_Type&, Read_Entry_Desc, Contig_Entry_Desc);

        friend std::ostream& operator << (std::ostream& os, const Contig_Entry& rhs);

    private:
        Seq_Type seq;
        std::list<Mutation> mut_list;
        std::list<Read_Chunk_Desc> chunk_dlist;
    };
}


#endif
