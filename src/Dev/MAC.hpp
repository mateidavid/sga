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
    /** Add a read.
     *
     * Create basic Read_Entry, Read_Chunk, and Contig_Entry objects,
     * initialize them, and place them in their respective containers.
     *
     * @param name_ptr Pointer to string with read name. (Read container takes ownership.)
     * @param seq_ptr Pointer to container with read sequence. (Contig container takes ownership.)
     * @param re_cont Container for read entry objects.
     * @param ce_cont Container for contig entry objects.
     */
    void add_read(const std::string* name_ptr, const Seq_Type* seq_ptr, Read_Entry_Cont& re_cont, Contig_Entry_Cont& ce_cont);
}


#endif
