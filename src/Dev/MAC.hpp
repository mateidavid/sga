//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __MAC_HPP
#define __MAC_HPP

#include "MAC_forward.hpp"
#include "Mutation.hpp"
#include "Read_Chunk.hpp"
#include "Read_Entry.hpp"
#include "Contig_Entry.hpp"


namespace MAC
{
    /**
     * @brief Process a read: create Read_Entry and Contig_Entry objects.
     * */
    void add_read(const std::string* name_ptr, const Seq_Type* seq_ptr, Read_Entry_Cont& re_cont, Contig_Entry_Cont& ce_cont);
}


#endif
