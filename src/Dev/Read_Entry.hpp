//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __READ_ENTRY_HPP
#define __READ_ENTRY_HPP

#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <functional>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>

#include "MAC_forward.hpp"
#include "Read_Chunk.hpp"


namespace MAC
{
    class Read_Entry
    {
    public:
        Read_Entry(const std::string* _name_ptr, const Seq_Type* seq_ptr) : name_ptr(_name_ptr)
        {
            assert(_name_ptr != NULL);
            assert(seq_ptr != NULL);
            Read_Chunk chunk(this, seq_ptr->size());
            chunk_cont.insert(chunk);
        }

        const std::string& get_name() const { return *name_ptr; }
        Read_Chunk_CPtr get_cptr_first_chunk() const { return &(*chunk_cont.begin()); }

        typedef const Seq_Type& key_type;
        key_type get_key() const { return *name_ptr; }

        void modify_read_chunk(Read_Chunk_CPtr chunk_cptr, Read_Chunk::modifier_type f)
        {
            Read_Chunk_Cont::iterator it = chunk_cont.iterator_to(*chunk_cptr);
            chunk_cont.modify(it, f);
        }

        friend std::ostream& operator << (std::ostream& os, const Read_Entry& rhs);

    private:
        std::shared_ptr<const std::string> name_ptr;
        Read_Chunk_Cont chunk_cont;
    };

    struct Read_Entry_Key
    {
        typedef Read_Entry::key_type result_type;
        result_type operator () (const Read_Entry& c) const { return c.get_key(); }
    };

    //typedef std::vector<Read_Entry> Read_Entry_Cont;
    typedef boost::multi_index_container<
      Read_Entry,
      boost::multi_index::indexed_by<
        boost::multi_index::hashed_unique<
          Read_Entry_Key
        >
      >
    > Read_Entry_Cont;
}


#endif
