//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __MAC_FORWARD_HPP
#define __MAC_FORWARD_HPP

#include <string>


/** Multi-Allelic Contig namespace */
namespace MAC
{
    /** Type for absolute and relative offsets in read and contig sequences. */
    typedef size_t Size_Type;
    /** Type holding a read sequence */
    typedef typename std::string Seq_Type;
    /** Type holding a single base pair */
    typedef typename Seq_Type::value_type Symb_Type;

    class Mutation;
    class Read_Chunk;
    class Read_Entry;
    class Contig_Entry;
}


#endif
