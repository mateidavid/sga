//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __READ_CHUNK_HPP
#define __READ_CHUNK_HPP

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <memory>
#include <functional>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/tuple/tuple.hpp>

#include "MAC_forward.hpp"
#include "Mutation.hpp"
#include "Cigar.hpp"


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

        /** @name Getters */
        /**@{*/
        /** Empty constructor. */
        Read_Chunk()
        : _re_ptr(NULL), _ce_ptr(NULL), _r_start(0), _r_len(0), _c_start(0), _c_len(0), _rc(false) {}

        /** Constructs a single read chunk from a read, and maps it to a new contig.
         * @param re_ptr Pointer to read entry where this chunk comes from.
         * @param len Length of the read.
         */
        Read_Chunk(const Read_Entry* re_ptr, Size_Type len)
        : _re_ptr(re_ptr), _ce_ptr(NULL), _r_start(0), _r_len(len), _c_start(0), _c_len(len), _rc(false) {}
        /**@}*/

        /** @name Getters */
        /**@{*/
        const Read_Entry* get_re_ptr() const { return _re_ptr;}
        const Contig_Entry* get_ce_ptr() const { return _ce_ptr; }
        Size_Type get_r_start() const { return _r_start; }
        Size_Type get_r_len() const { return _r_len; }
        Size_Type get_r_end() const { return _r_start + _r_len; }
        Size_Type get_c_start() const { return _c_start; }
        Size_Type get_c_len() const { return _c_len; }
        Size_Type get_c_end() const { return _c_start + _c_len; }
        bool get_rc() const { return _rc; }
        const std::vector< const Mutation* > get_mut_ptr_cont() const { return _mut_ptr_cont; }
        key_type get_key() const { return _r_start; }
        /**@}*/

        /** Set read entry pointer. */
        void set_re_ptr(const Read_Entry* re_ptr) { _re_ptr = re_ptr; }

        /** Set contig entry pointer. */
        void set_ce_ptr(const Contig_Entry* ce_ptr) { _ce_ptr = ce_ptr; }

        /** Check if read chunk has given mutation.
         * @param mut_cptr Pointer to mutation to look for.
         */
        bool have_mutation(const Mutation* mut_cptr) const;

        /** If read chunk has old mutation, add new mutation right after it.
         * @param m_old_cptr Old mutation to look for.
         * @param m_new_cptr New mutation to add if old one exists.
         */
        void cond_add_mutation(const Mutation* m_old_cptr, const Mutation* m_new_cptr);

        /** Assign read chunk to contig.
         * @param ce_cptr Pointer to contig entry object.
         * @param c_start Start of contig base sequence this chunk is mapped to.
         * @param c_len Length of contig base sequence this chunk is mapped to.
         * @param rc Orientation of the mapping.
         * @param mut_ptr_cont List of contig mutations observed by this chunk.
         */
        void assign_to_contig(const Contig_Entry* ce_cptr, Size_Type c_start, Size_Type c_len, bool rc, const std::vector< const Mutation* > mut_ptr_cont)
        {
            _ce_ptr = ce_cptr;
            _c_start = c_start;
            _c_len = c_len;
            _rc = rc;
            _mut_ptr_cont = mut_ptr_cont;
        }

        /** Prepare data necessary to perform a split at the given position.
         * @param brk Breakpoint position, 0-based; equivalently, length before breakpoint.
         * @param is_contig_brk True if position is on contig, false if it is on read.
         * @return A triple (c_pos, r_pos, idx).
         * Guarantees:
         * - if is_contig_brk: c_pos <= brk
         * - else: if not _rc: r_pos <= brk; else brk <= r_pos
         * - [c_start,c_pos) matched to [r_start,r_pos) / [r_pos,r_end)
         * - if X_pos == brk: idx is the index of the first (smallest m_c_start) mutation with m_c_start >= c_pos.
         * - if X_pos != brk: idx is the index of a mutation spanning both positions brk-1 and brk
         */
         boost::tuple< Size_Type, Size_Type, size_t > get_cut_data(Size_Type r_brk, bool is_contig_brk) const;

         /** If given breakpoint is spanned by a mutation, return it along with cut coordinates.
          * @param brk Position where to cut, 0-based. Alternatively, length before the cut.
          * @param is_contig_brk True if position is on contig; false if it is on read.
          * @param c_brk_cont Set of preferred breakpoints. (Optional)
          * @return Pointer to mutation to cut (NULL if none), base and alternate offsets where to cut.
          */
         boost::tuple< const Mutation*, Size_Type, Size_Type > get_mutation_to_cut(
             Size_Type brk, bool is_contig_brk, const std::set< Size_Type >& c_brk_cont = std::set< Size_Type >()) const;

         /** Split read chunk based on contig position and mutation map.
          * @param c_brk Contig breakpoint.
          * @param mut_cptr_map Map of mutations which go to the right side.
          * @param ce_cptr Pointer to new contig entry object.
          * @return Flag whether to move read chunk to the right contig side; additional Read_Chunk object mapped to right side, if the chunk gets split
          */
         boost::tuple< bool, std::shared_ptr< Read_Chunk > > apply_contig_split(
             Size_Type c_brk, const std::map< const Mutation*, const Mutation* >& mut_cptr_map, const Contig_Entry* ce_cptr);

         /** Given a read chunk breakpoint, compute the range of possible contig breakpoints.
          * @param pos Breakpoint position in the read, 0-based; equivalently, length of read before breakpoint.
          * @return Leftmost and rightmost breakpoint positions in contig.
          */
         std::pair< Size_Type, Size_Type > get_contig_brk_range(Size_Type pos);

         /** Create 2 Read_Chunk objects (and corresponding Mutation containers) from a cigar string.
         * @param r1_start Overlap start in r1, 0-based. Equivalently, length of r1 before the overlap.
         * @param r2_start Overlap start in r2, 0-based. Equivalently, length of r2 before the overlap.
         * @param cigar Cigar string (r1:reference, r2:query).
         * @return A vector of size 2; index 0 corresponds to (base r1, query r2),
         * index 1 corresponds to (base r2, query r1).
         *
        static std::vector<std::pair<Read_Chunk,Mutation_Cont>> make_chunks_from_cigar(
            Size_Type r1_start, Size_Type r2_start, const Cigar& cigar);
        */

         /** Collapse mutations corresponding to 2 mappings.
          * @param lhs Read_Chunk object corresponding to contig->read1 mapping.
          * @param lhs_me Container of Mutation_Extra objects from lhs.
          * @param rhs Read_Chunk object corresponding to read1->read2 mapping.
          * @return vector of contig->read2 mutations; for each:
          * a bool specifying if it is a contig mutation (true) or a read1 mutation (false);
          * pointer to the mutation object;
          * and a mutation start&length.
          */
         static std::vector< boost::tuple< bool, Read_Chunk_CPtr, Size_Type, Size_Type > > collapse_mutations(
             const Read_Chunk& rc1, const Mutation_Extra_Cont& rc1_me_cont, const Read_Chunk& rc2);

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

        /** Construct mutation ptr container with read positions. */
        Mutation_Extra_Cont get_me_cont() const;
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
