//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __MAC_FORWARD_HPP
#define __MAC_FORWARD_HPP

#include <iostream>
#include <string>
#include <list>
#include <algorithm>


namespace MAC
{
    typedef size_t Size_Type;
    typedef typename std::string Seq_Type;
    typedef typename Seq_Type::value_type Symb_Type;

    class Read_Entry;
    class Contig_Entry;
    class Read_Chunk;
    class Mutation;

    typedef Read_Entry* Read_Entry_Desc;
    typedef Contig_Entry* Contig_Entry_Desc;
    typedef Read_Chunk* Read_Chunk_Desc;
    typedef Mutation* Mutation_Desc;
}


#endif
