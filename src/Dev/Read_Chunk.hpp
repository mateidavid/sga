//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __READ_CHUNK_HPP
#define __READ_CHUNK_HPP

#include <iostream>
#include <vector>
#include <functional>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>

#include "MAC_forward.hpp"


namespace MAC
{
    /** Holds information about a read chunk.
     *
     * Reads will be partitioned into chunks of non-zero size as they are mapped to contigs.
     * Each read chunk object contains information about the read it comes from (start, span),
     * the contig it is mapped to (start, span, orientation), and the contig mutations it observes.
     */
    class Read_Chunk
    {
    public:
        /** Type of key use to store Read_Chunk objects. */
        typedef Size_Type key_type;

        /** Type for an external unary modifier. */
        typedef std::function<void(Read_Chunk&)> modifier_type;

        //Read_Chunk()
        //: r_start(0), r_len(0), c_start(0), c_len(0) {}

        /** Constructs a single read chunk from a read, and maps it to a new contig.
         * @param re_ptr Pointer to read entry where this chunk comes from.
         * @param len Length of the read.
         */
        Read_Chunk(const Read_Entry* re_ptr, Size_Type len)
        : _re_ptr(re_ptr), _ce_ptr(NULL), _r_start(0), _r_len(len), _c_start(0), _c_len(len), _rc(false) {}

        /** @name Getters */
        /**@{*/
        const Read_Entry* get_re_ptr() const { return _re_ptr;}
        const Contig_Entry* get_ce_ptr_desc() const { return _ce_ptr; }
        Size_Type get_r_start() const { return _r_start; }
        Size_Type get_r_len() const { return _r_len; }
        Size_Type get_r_end() const { return _r_start + _r_len; }
        Size_Type get_c_start() const { return _c_start; }
        Size_Type get_c_len() const { return _c_len; }
        Size_Type get_c_end() const { return _c_start + _c_len; }
        bool get_rc() const { return _rc; }
        key_type get_key() const { return _r_start; }
        /**@}*/

        void set_ce_ptr(const Contig_Entry* ce_ptr) { _ce_ptr = ce_ptr; }

        /** Assign read chunk to contig.
         * @param ce_ptr Pointer to contig entry object.
         * @param c_start Start of contig base sequence this chunk is mapped to.
         * @param c_len Length of contig base sequence this chunk is mapped to.
         * @param rc Orientation of the mapping.
         * @param mut_ptr_cont List of contig mutations observed by this chunk.
         */
        void assign_to_contig(const Contig_Entry* ce_ptr, Size_Type c_start, Size_Type c_len, bool rc, const std::vector<const Mutation*> mut_ptr_cont)
        {
            _ce_ptr = ce_ptr;
            _c_start = c_start;
            _c_len = c_len;
            _rc = rc;
            _mut_ptr_cont = mut_ptr_cont;
        }

        friend std::ostream& operator << (std::ostream&, const Read_Chunk&);

    private:
        const Read_Entry* _re_ptr;
        const Contig_Entry* _ce_ptr;
        Size_Type _r_start;
        Size_Type _r_len;
        Size_Type _c_start;
        Size_Type _c_len;
        bool _rc;
        std::vector<const Mutation*> _mut_ptr_cont;
    };

    /** Pointer to constant Read_Chunk. */
    typedef const Read_Chunk* Read_Chunk_CPtr;

    namespace detail
    {
        /** Key extractor struct for boost::multi_index_container. */
        struct Read_Chunk_Key
        {
            typedef Read_Chunk::key_type result_type;
            result_type operator () (const Read_Chunk& c) const { return c.get_key(); }
            result_type operator () (const Read_Chunk_CPtr& c) const { return c->get_key(); }
        };
    }

    /** Container for Read_Chunk objects. */
    typedef boost::multi_index_container<
      Read_Chunk,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<
          detail::Read_Chunk_Key
        >
      >
    > Read_Chunk_Cont;

    /** Container for pointers to Read_Chunk objects. */
    typedef boost::multi_index_container<
      Read_Chunk_CPtr,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_non_unique<
          detail::Read_Chunk_Key
        >
      >
    > Read_Chunk_CPtr_Cont;
}


#endif
