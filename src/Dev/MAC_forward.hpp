//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __MAC_FORWARD_HPP
#define __MAC_FORWARD_HPP

#include <string>


namespace MAC
{
    typedef size_t Size_Type;
    typedef typename std::string Seq_Type;
    typedef typename Seq_Type::value_type Symb_Type;

    class Mutation;
    class Read_Chunk;
    class Read_Entry;
    class Contig_Entry;
}


#endif
