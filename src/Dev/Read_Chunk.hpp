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

#include "MAC_forward.hpp"
#include "Mutation.hpp"
#include "Cigar.hpp"


namespace MAC
{
    /** Set of coordinates used to traverse a Read_Chunk mapping. */
    struct Read_Chunk_Pos
    {
        /** Contig position. */
        Size_Type c_pos;
        /** Read position. */
        Size_Type r_pos;
        /** Number of mutations passed. */
        size_t mut_idx;
        /** If non-zero, offset into mutation mut_idx. */
        size_t mut_offset;

        /** Constructor. */
        Read_Chunk_Pos(Size_Type _c_pos = 0, Size_Type _r_pos = 0, size_t _mut_idx = 0, size_t _mut_offset = 0)
        : c_pos(_c_pos), r_pos(_r_pos), mut_idx(_mut_idx), mut_offset(_mut_offset) {}

        /** Comparison operator. */
        bool operator == (const Read_Chunk_Pos& rhs)
        {
            return c_pos == rhs.c_pos and r_pos == rhs.r_pos and mut_idx == rhs.mut_idx and mut_offset == rhs.mut_offset;
        }
        bool operator != (const Read_Chunk_Pos& rhs) { return not (*this == rhs); }
    };

    std::ostream& operator << (std::ostream&, const Read_Chunk_Pos&);

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
        Seq_Type get_seq() const;
        Seq_Type substr(Size_Type start, Size_Type len) const;
        /**@}*/

        /** @name Traversal using Read_Chunk_Pos objects */
        /**@{*/
        /** Asserts that Read_Chunk_Pos object is valid. */
        inline bool check_pos(const Read_Chunk_Pos& pos) const;

        /** Get mapping start position. */
        Read_Chunk_Pos get_start_pos() const { return Read_Chunk_Pos(get_c_start(), (not _rc? get_r_start() : get_r_end()), 0, 0); }

        /** Get mapping end position. */
        Read_Chunk_Pos get_end_pos() const { return Read_Chunk_Pos(get_c_end(), (not _rc? get_r_end() : get_r_start()), _mut_ptr_cont.size(), 0); }

        /** Get length of unmutated mapped stretch starting at given position.
         * @param pos Read_Chunk_Pos object.
         * @param forward Bool; true if looking for match stretch forward, false if backward.
         * @return The length of the next stretch of unmutated bases; 0 if on the breakpoint of a mutation.
         */
        inline Size_Type get_match_len_from_pos(const Read_Chunk_Pos& pos, bool forward = true) const;

        /** Increment position object past the next mapping stretch, stopping at the given breakpoint.
         * @param pos Read_Chunk_Pos position object.
         * @param brk Breakpoint position.
         * @param on_contig Bool; true if breakpoint is on contig, false if breakpoint on read.
         */
        void increment_pos(Read_Chunk_Pos& pos, Size_Type brk = 0, bool on_contig = true) const;

        /** Decrement position object past the next mapping stretch, stopping at the given breakpoint.
         * @param pos Read_Chunk_Pos position object.
         * @param brk Breakpoint position.
         * @param on_contig Bool; true if breakpoint is on contig, false if breakpoint on read.
         */
        void decrement_pos(Read_Chunk_Pos& pos, Size_Type brk = 0, bool on_contig = true) const;

        /** Advance position object: wrapper for increment and decrement.
         * @param pos Read_Chunk_Pos position object.
         * @param forward Direction; true: increment; false: decrement.
         * @param brk Breakpoint position.
         * @param on_contig Bool; true if breakpoint is on contig, false if breakpoint on read.
         */
        void advance_pos(Read_Chunk_Pos& pos, bool forward, Size_Type brk = 0, bool on_contig = true) const
        { if (forward) increment_pos(pos, brk, on_contig); else decrement_pos(pos, brk, on_contig); }

        /** Advance position, but using a read Mutation breakpoint.
         * Repeated calls to this function will produce mapping stretches prior to the read Mutation,
         * including sliced contig mutations, followed by the stretch corresponding to the read Mutation.
         * @param pos Position to advance.
         * @param mut Read Mutation to use as breakpoint.
         * @param forward Direction; true: increment; false: decrement.
         * @return Bool; true if breakpoint reached (i.e., mapping stretch corresponds to read Mutation), false ow.
         */
        bool advance_pos_til_mut(Read_Chunk_Pos& pos, const Mutation& mut, bool forward = true) const;
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
        void assign_to_contig(const Contig_Entry* ce_cptr, Size_Type c_start, Size_Type c_len, bool rc, const std::vector< const Mutation* >& mut_ptr_cont)
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
         std::tuple< Size_Type, Size_Type, size_t > get_cut_data(Size_Type r_brk, bool is_contig_brk) const;

         /** If given breakpoint is spanned by a mutation, return it along with cut coordinates.
          * @param brk Position where to cut, 0-based. Alternatively, length before the cut.
          * @param is_contig_brk True if position is on contig; false if it is on read.
          * @param c_brk_cont Set of preferred breakpoints. (Optional)
          * @return Pointer to mutation to cut (NULL if none), base and alternate offsets where to cut.
          */
         std::tuple< const Mutation*, Size_Type, Size_Type > get_mutation_to_cut(
             Size_Type brk, bool is_contig_brk, const std::set< Size_Type >& c_brk_cont = std::set< Size_Type >()) const;

         /** Split read chunk based on contig position and mutation map.
          * @param c_brk Contig breakpoint.
          * @param mut_cptr_map Map of mutations which go to the right side.
          * @param ce_cptr Pointer to new contig entry object.
          * @return Flag whether to move read chunk to the right contig side; additional Read_Chunk object mapped to right side, if the chunk gets split
          */
         std::tuple< bool, std::shared_ptr< Read_Chunk > > apply_contig_split(
             Size_Type c_brk, const std::map< const Mutation*, const Mutation* >& mut_cptr_map, const Contig_Entry* ce_cptr);

         /** Given a read chunk breakpoint, compute the range of possible contig breakpoints.
          * @param pos Breakpoint position in the read, 0-based; equivalently, length of read before breakpoint.
          * @return Leftmost and rightmost breakpoint positions in contig.
          */
         std::pair< Size_Type, Size_Type > get_contig_brk_range(Size_Type pos);

         /** Create Read_Chunk object and corresponding Mutation container from a cigar string.
         * @param cigar Cigar string.
         * @param qr Query string, optional;
         * if not given (empty), Mutation objects hold alternate sequence _lengths_ only;
         * if given, it can be either the entire query string, or just the part mapped by cigar.
         * @return A read chunk corresponding to the query in the cigar object, and a container for Mutation objects.
         */
         static std::tuple< std::shared_ptr< Read_Chunk >, std::shared_ptr< Mutation_Cont > > make_chunk_from_cigar(
             const Cigar& cigar, const std::string& qr = std::string());

         /** Construct a Read_Chunk corresponding to the read->contig mapping.
          * @return A new Read_Chunk object; a container for reversed Mutations;
          * and a Mutation translation object from original to reversed Mutations.
          */
         std::tuple< std::shared_ptr< Read_Chunk >, std::shared_ptr< Contig_Entry >, std::shared_ptr< Mutation_Trans_Cont > >
         reverse() const;

         /** Translate read mutations into contig mutations.
          * @param r_mut_cont Container with read mutations.
          * @return A Mutation translation container, and a container for new mutations.
          */
         std::tuple< std::shared_ptr< Mutation_Trans_Cont >, std::shared_ptr< Mutation_Cont > > lift_read_mutations(
             const Mutation_Cont& r_mut_cont) const;

         /** Compute old contig mutations observed by a read mapping.
          * @param r_mut_cptr_cont List of read mutations observed by mapping.
          * @param mut_map Map of read mutations to new contig mutations.
          * @param new_mut_cont_sptr Container for new Mutations; this function might add to this container any merged insertions.
          * @return A vector of: (mutation ptr, start, end, bool);
          * bool: true iff this is a translated read mutation;
          * start: the start of the slice (=0 for read mutations);
          * end: end of the slice, 0 if slice goes to the end of the mutation (=0 for read mutations).
          */
         std::shared_ptr< std::vector< std::tuple< Mutation_CPtr, Size_Type, Size_Type, bool > > >
         get_mutations_under_mapping(const std::vector< Mutation_CPtr >& r_mut_cptr_cont,
                                     const Mutation_Trans_Cont& mut_map,
                                     std::shared_ptr< Mutation_Cont > new_mut_cont_sptr) const;

         /** Collapse mutations corresponding to 2 mappings.
          * @param lhs Read_Chunk object corresponding to contig->read1 mapping.
          * @param lhs_me Container of Mutation_Extra objects from lhs.
          * @param rhs Read_Chunk object corresponding to read1->read2 mapping.
          * @return vector of contig->read2 mutations; for each:
          * a bool specifying if it is a contig mutation (true) or a read1 mutation (false);
          * pointer to the mutation object;
          * and a mutation start&length.
          *
         static std::vector< std::tuple< bool, Read_Chunk_CPtr, Size_Type, Size_Type > > collapse_mutations(
             const Read_Chunk& rc1, const Mutation_Extra_Cont& rc1_me_cont, const Read_Chunk& rc2);
         */

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

        /** Compute read positions of the contig mutations in this object. */
        std::vector< Size_Type > get_mut_pos() const;
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
