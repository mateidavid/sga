//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __CONTIG_ENTRY_HPP
#define __CONTIG_ENTRY_HPP

#include <string>
#include <iostream>
#include <memory>
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
        Contig_Entry(const Seq_Type* _seq_ptr) : seq_ptr(_seq_ptr) {}

        const Seq_Type& get_seq() const { return *seq_ptr; }

        void add_chunk(Read_Chunk_CPtr chunk_cptr) { chunk_cptr_cont.insert(chunk_cptr); }

        friend std::ostream& operator << (std::ostream& os, const Contig_Entry& rhs);

    private:
        std::shared_ptr<const Seq_Type> seq_ptr;
        Mutation_Cont mut_cont;
        Read_Chunk_CPtr_Cont chunk_cptr_cont;
    };

    typedef boost::multi_index_container<
      Contig_Entry,
      boost::multi_index::indexed_by<
        boost::multi_index::sequenced<>
      >
    > Contig_Entry_Cont;
}


#endif
