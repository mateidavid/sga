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
    /** Holds information about a read.
     *
     * The read sequence gets assigned to various contigs. The Read_Entry object holds the read name,
     * as well as the (order) sequence of read chunks which are mapped to contigs.
     */
    class Read_Entry
    {
    public:
        /** Type of key use to store Read_Entry objects. */
        typedef const Seq_Type& key_type;

        /** Constructor.
         * @param name_ptr Pointer to string containing read name.
         * @param len Length of the read.
         */
        Read_Entry(const std::string* name_ptr, Size_Type len) : _name_ptr(name_ptr)
        {
            assert(name_ptr != NULL);
            Read_Chunk chunk(this, len);
            _chunk_cont.insert(chunk);
        }

        /** @name Getters */
        /**@{*/
        const std::string& get_name() const { return *_name_ptr; }
        Size_Type get_len() const { return _chunk_cont.rend()->get_r_end(); }
        key_type get_key() const { return *_name_ptr; }
        const Read_Chunk_Cont& get_chunk_cont() const { return _chunk_cont; }
        /**@}*/

        /** Find chunk which contains given read position.
         * @param r_pos Read position, 0-based.
         * @return Pointer to Read Chunk object, or NULL if no chunk contains r_pos.
         */
        const Read_Chunk* get_chunk_with_pos(Size_Type r_pos) const;

        /** Member read chunk modifier.
         * @param chunk_cptr Pointer to read chunk to modify.
         * @param f Unary modifier function to apply.
         */
        void modify_read_chunk(Read_Chunk_CPtr chunk_cptr, Read_Chunk::modifier_type f)
        {
            Read_Chunk_Cont::iterator it = _chunk_cont.iterator_to(*chunk_cptr);
            _chunk_cont.modify(it, f);
        }

        friend std::ostream& operator << (std::ostream& os, const Read_Entry& rhs);

    private:
        std::shared_ptr<const std::string> _name_ptr;
        Read_Chunk_Cont _chunk_cont;
    };

    namespace detail
    {
        /** Key extractor struct for boost::multi_index_container. */
        struct Read_Entry_Key
        {
            typedef Read_Entry::key_type result_type;
            result_type operator () (const Read_Entry& c) const { return c.get_key(); }
        };
    }

    /** Container for Read_Entry objects. */
    typedef boost::multi_index_container<
      Read_Entry,
      boost::multi_index::indexed_by<
        boost::multi_index::hashed_unique<
          detail::Read_Entry_Key
        >
      >
    > Read_Entry_Cont;
}


#endif
