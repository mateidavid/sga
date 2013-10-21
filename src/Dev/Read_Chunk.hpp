//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __READ_CHUNK_HPP
#define __READ_CHUNK_HPP

#include "MAC_forward.hpp"


namespace MAC
{
    class Read_Chunk
    {
    public:
        Read_Chunk()
        : r_start(0), r_len(0), c_start(0), c_len(0) {}

        Read_Chunk(Read_Entry_Desc _re_desc, Contig_Entry_Desc _ce_desc, Size_Type len)
        : re_desc(_re_desc), ce_desc(_ce_desc), r_start(0), r_len(len), c_start(0), c_len(len), rc(false) {}

        Read_Entry_Desc get_read_entry_desc() const { return re_desc;}
        Contig_Entry_Desc get_contig_entry_desc() const { return ce_desc; }

        friend std::ostream& operator << (std::ostream& os, const Read_Chunk& rhs);

    private:
        Read_Entry_Desc re_desc;
        Contig_Entry_Desc ce_desc;
        Size_Type r_start, r_len;
        Size_Type c_start, c_len;
        bool rc;
        std::list<Mutation_Desc> mut_dlist;
    };
}


#endif
