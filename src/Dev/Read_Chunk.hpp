//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __READ_CHUNK_HPP
#define __READ_CHUNK_HPP

#include <iostream>
#include <vector>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>

#include "MAC_forward.hpp"


namespace MAC
{
    class Read_Chunk
    {
    public:
        Read_Chunk()
        : r_start(0), r_len(0), c_start(0), c_len(0) {}

        Read_Chunk(Read_Entry* _re_ptr, Contig_Entry* _ce_ptr, Size_Type len)
        : re_ptr(_re_ptr), ce_ptr(_ce_ptr), r_start(0), r_len(len), c_start(0), c_len(len), rc(false) {}

        Read_Entry* get_re_ptr() const { return re_ptr;}
        Contig_Entry* get_ce_ptr_desc() const { return ce_ptr; }

        typedef Size_Type key_type;
        key_type get_key() const { return r_start; }

        friend std::ostream& operator << (std::ostream& os, const Read_Chunk& rhs);

    private:
        Read_Entry* re_ptr;
        Contig_Entry* ce_ptr;
        Size_Type r_start;
        Size_Type r_len;
        Size_Type c_start;
        Size_Type c_len;
        bool rc;
        std::vector<Mutation*> mut_ptr_cont;
    };

    typedef const Read_Chunk* Read_Chunk_CPtr;

    struct Read_Chunk_Key
    {
        typedef Read_Chunk::key_type result_type;
        result_type operator () (const Read_Chunk& c) const { return c.get_key(); }
        result_type operator () (const Read_Chunk_CPtr& c) const { return c->get_key(); }
    };

    typedef boost::multi_index_container<
      Read_Chunk,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<
          Read_Chunk_Key
        >
      >
    > Read_Chunk_Cont;

    typedef boost::multi_index_container<
      Read_Chunk_CPtr,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_non_unique<
          Read_Chunk_Key
        >
      >
    > Read_Chunk_Ptr_Cont;
}


#endif
