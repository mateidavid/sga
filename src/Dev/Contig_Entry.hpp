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
#include <memory>
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
        /** Constructor.
         * @param seq_ptr Pointer to read sequence.
         */
        Contig_Entry(const Seq_Type* seq_ptr) : _seq_ptr(seq_ptr) {}

        /** @name Getters */
        /**@{*/
        const Seq_Type& get_seq() const { return *_seq_ptr; }
        const Mutation_Cont& get_mutation_cont() const { return _mut_cont; }
        /**@}*/

        /** Add to the list of read chunks.
         * @param chunk_cptr Pointer to read chunk object.
         */
        void add_chunk(Read_Chunk_CPtr chunk_cptr) { _chunk_cptr_cont.insert(chunk_cptr); }

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

        friend std::ostream& operator << (std::ostream& os, const Contig_Entry& rhs);

    private:
        std::shared_ptr<const Seq_Type> _seq_ptr;
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
