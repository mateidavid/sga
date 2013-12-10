//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __CONTIG_ENTRY_HPP
#define __CONTIG_ENTRY_HPP

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <memory>
#include <functional>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/sequenced_index.hpp>

#include "MAC_forward.hpp"
#include "Mutation.hpp"
#include "Read_Chunk.hpp"


namespace MAC
{
    /** Holds information about a contig.
     *
     * Each object holds a base sequence, a list of observed mutations, and a list of pointers
     * to read chunks which are mapped to this contig.
     */
    class Contig_Entry
    {
    public:
        /** Type for an external unary modifier. */
        typedef std::function<void(Contig_Entry&)> modifier_type;

        /** Constructor.
         * @param seq_ptr Pointer to read sequence.
         */
        Contig_Entry(Seq_Type* seq_ptr, Size_Type seq_offset = 0) : _seq_ptr(seq_ptr), _seq_offset(seq_offset) {}

        /** @name Getters */
        /**@{*/
        const Seq_Type& get_seq() const { return *_seq_ptr; }
        Size_Type get_seq_offset() const { return _seq_offset; }
        Seq_Type substr(Size_Type start, Size_Type len) const
        { assert(start >= _seq_offset and start + len <= _seq_offset + _seq_ptr->size()); return _seq_ptr->substr(start - _seq_offset, len); }
        Size_Type get_len() const { return _seq_ptr->size(); }
        const Mutation_Cont& get_mut_cont() const { return _mut_cont; }
        Mutation_Cont& mut_cont() { return _mut_cont; }
        const Read_Chunk_CPtr_Cont& get_chunk_cptr_cont() const { return _chunk_cptr_cont; }
        /**@}*/

        /** Add to the list of read chunks.
         * @param rc_cptr Pointer to read chunk to add.
         */
        void add_chunk(Read_Chunk_CPtr rc_cptr) { _chunk_cptr_cont.push_back(rc_cptr); }

        /** Remove read chunk.
         * @param rc_cptr Pointer to read chunk to remove.
         */
        void remove_chunk(Read_Chunk_CPtr rc_cptr);

        /** Cut mutation at given offsets.
         * @param mut_cptr Pointer to mutation to cut.
         * @param c_offset Contig offset inside mutation where to cut.
         * @param r_offset Read offset inside mutation where to cut.
         * @return Pointer to mutation object containing leftover part of the original mutation.
         */
        const Mutation* cut_mutation(const Mutation* mut_cptr, Size_Type c_offset, Size_Type r_offset);

        /** Get chunks which contain a given mutation.
         * @param mut_cptr Pointer to mutation.
         * @return Vector of pointers to chunks that have this mutation.
         */
        std::vector< Read_Chunk_CPtr > get_chunks_with_mutation(const Mutation* mut_cptr) const;

        /** Get mutations that completely span a given contig position.
         * @param c_pos Contig position, 0-based.
         * @return Vector of pointers to Mutation objects completely spanning c_pos.
         */
        std::vector< const Mutation* > get_mutations_spanning_pos(Size_Type c_pos) const;

        /** Get mutations from the second half of a given Contig Entry object that is being cut in 2.
         * @param ce_cptr Contig_Entry being cut.
         * @param c_brk Position of the cut.
         * @param mut_left_cptr Pointer to insertion at c_pos that must appear on the left of the cut, if any.
         * @return Map with (key=old mutation cptr; value=new mutation cptr).
         */
        std::map< const Mutation*, const Mutation* > acquire_second_half_mutations(
            const Contig_Entry* ce_cptr, Size_Type c_brk, const Mutation* mut_left_cptr);

        /** Drop mutations that appear in the map.
         * @param mut_cptr_map Map produced by acquire_second_half_mutations().
         */
        void drop_mutations(const std::map< const Mutation*, const Mutation* >& mut_cptr_map);

        /** Drop unused mutations. */
        void drop_unused_mutations();

        /** Drop base sequence suffix.
         * @param c_brk Prefix length to keep.
         */
        void drop_base_seq(Size_Type c_brk) { _seq_ptr->resize(c_brk); }

        /** Integrity check. */
        bool check() const;

        friend std::ostream& operator << (std::ostream& os, const Contig_Entry& rhs);

    private:
        std::shared_ptr< Seq_Type> _seq_ptr;
        Size_Type _seq_offset;
        Mutation_Cont _mut_cont;
        Read_Chunk_CPtr_Cont _chunk_cptr_cont;
    };

    /** Container for Contig_Entry objects. */
    typedef boost::multi_index_container<
      Contig_Entry,
      boost::multi_index::indexed_by<
        boost::multi_index::sequenced<>
      >
    > Contig_Entry_Cont;
}


#endif
