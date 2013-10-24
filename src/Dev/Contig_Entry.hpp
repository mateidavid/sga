//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __CONTIG_ENTRY_HPP
#define __CONTIG_ENTRY_HPP

#include <string>
#include <iostream>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/sequenced_index.hpp>

#include "MAC_forward.hpp"
#include "Mutation.hpp"
#include "Read_Chunk.hpp"


namespace MAC
{
    class Contig_Entry
    {
    public:
        Contig_Entry() {}

        friend void add_read(const std::string&, const Seq_Type&, Read_Entry*, Contig_Entry*);

        friend std::ostream& operator << (std::ostream& os, const Contig_Entry& rhs);

    private:
        Seq_Type seq;
        Mutation_Cont mut_cont;
        Read_Chunk_Ptr_Cont chunk_ptr_cont;
    };

    typedef boost::multi_index_container<
      Contig_Entry,
      boost::multi_index::indexed_by<
        boost::multi_index::sequenced<>
      >
    > Contig_Entry_Cont;
}


#endif
