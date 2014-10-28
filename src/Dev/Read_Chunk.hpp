//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __READ_CHUNK_HPP
#define __READ_CHUNK_HPP

#include "MAC_forward.hpp"
#include "Mutation.hpp"
#include "Mutation_Ptr_Cont.hpp"


namespace MAC
{

/** Struct that represents a Read_Chunk position. */
struct Read_Chunk_Pos
{
public:
    /** Contig position. */
    Size_Type c_pos;
    /** Read position. */
    Size_Type r_pos;
    /** Number of mutations passed. */
    //size_t mut_idx;
    /** Iterator with next mutation pointer. */
    Mutation_Ptr_Cont::const_iterator mca_cit;
    /** If non-zero, offset into current mutation. */
    Size_Type mut_offset;
private:
    /** Read_Chunk parent object. */
    const Read_Chunk* rc_cptr;

private:
    // Constructed only by Read_Chunk objects
    friend class Read_Chunk;
    /** Constructor. */
    Read_Chunk_Pos(Size_Type _c_pos = 0,
                   Size_Type _r_pos = 0,
                   Mutation_Ptr_Cont::const_iterator _mca_cit = Mutation_Ptr_Cont::const_iterator(),
                   Size_Type _mut_offset = 0,
                   const Read_Chunk* _rc_cptr = NULL)
        : c_pos(_c_pos),
          r_pos(_r_pos),
          mca_cit(_mca_cit),
          mut_offset(_mut_offset),
          rc_cptr(_rc_cptr) {}

public:
    // allow copy and move
    DEFAULT_COPY_CTOR(Read_Chunk_Pos)
    DEFAULT_MOVE_CTOR(Read_Chunk_Pos)
    DEFAULT_COPY_ASOP(Read_Chunk_Pos)
    DEFAULT_MOVE_ASOP(Read_Chunk_Pos)

public:
    /** Check if position is past the last mutation in the read chunk. */
    bool past_last_mut() const;

    /** Check if position is past first mutation in the read chunk. */
    bool past_first_mut() const;

    /** Get current or next mutation. */
    const Mutation& mut() const
    {
        ASSERT(rc_cptr);
        ASSERT(not past_last_mut());
        return *(mca_cit->mut_cbptr());
    }

    /** Get previous mutation. */
    const Mutation& prev_mut() const
    {
        ASSERT(rc_cptr);
        ASSERT(past_first_mut());
        auto tmp_cit = mca_cit;
        return *((--tmp_cit)->mut_cbptr());
    }

    /** Get length of unmutated mapped stretch starting at given position.
     * @param forward Bool; true if looking for match stretch forward, false if backward.
     * @return The length of the next stretch of unmutated bases; 0 if on the breakpoint of a mutation.
     */
    inline Size_Type get_match_len(bool forward = true) const;

    /** Increment position object past the next mapping stretch, stopping at the given breakpoint.
     * @param brk Breakpoint position.
     * @param on_contig Bool; true if breakpoint is on contig, false if breakpoint on read.
     */
    void increment(Size_Type brk = 0, bool on_contig = true);

    /** Decrement position object past the next mapping stretch, stopping at the given breakpoint.
     * @param brk Breakpoint position.
     * @param on_contig Bool; true if breakpoint is on contig, false if breakpoint on read.
     */
    void decrement(Size_Type brk = 0, bool on_contig = true);

    /** Advance position object: wrapper for increment and decrement.
     * @param forward Direction; true: increment; false: decrement.
     * @param brk Breakpoint position.
     * @param on_contig Bool; true if breakpoint is on contig, false if breakpoint on read.
     */
    void advance(bool forward, Size_Type brk = 0, bool on_contig = true)
    { if (forward) { increment(brk, on_contig); } else { decrement(brk, on_contig); } }

    /** Get position corresponding to given breakpoint.
     * @param brk Breakpoint location.
     * @param on_contig True: contig position; False: read position.
     */
    Read_Chunk_Pos& jump_to_brk(Size_Type brk, bool on_contig);

    /** Advance position, but using a read Mutation breakpoint.
     * Repeated calls to this function will produce mapping stretches prior to the read Mutation,
     * including sliced contig mutations, followed by the stretch corresponding to the read Mutation.
     * @param mut Read Mutation to use as breakpoint.
     * @param forward Direction; true: increment; false: decrement.
     * @return Bool; true if breakpoint reached (i.e., mapping stretch corresponds to read Mutation), false ow.
     */
    bool advance_til_mut(const Mutation& mut, bool forward = true);

    /** Advance position past any deletions.
     * @param forward Direction.
     */
    void advance_past_del(bool forward = true);

    /** Asserts that Pos object is valid. */
    void check() const;

    /** To stream. */
    //friend std::ostream& operator << (std::ostream&, const Read_Chunk_Pos&);
    boost::property_tree::ptree to_ptree() const;
};

bool operator == (const Read_Chunk_Pos&, const Read_Chunk_Pos&);
bool operator != (const Read_Chunk_Pos&, const Read_Chunk_Pos&);

namespace detail
{

struct Read_Chunk_ITree_Node_Traits;
struct Read_Chunk_Set_Node_Traits;

}

/** Holds information about a read chunk.
 *
 * Reads will be partitioned into chunks of non-zero size as they are mapped to contigs.
 * Each read chunk object contains information about the read it comes from (start, span),
 * the contig it is mapped to (start, span, orientation), and the contig mutations it observes.
 */
class Read_Chunk
{
public:
    typedef Read_Chunk_Pos Pos;
private:
    // Can only be constructed by Factory object
    friend class bounded::Factory< Read_Chunk >;

    /** @name Constructors */
    /**@{*/
    /** Empty constructor. */
    Read_Chunk()
        : _r_start(0),
          _r_len(0),
          _c_start(0),
          _c_len(0),
          _re_bptr(nullptr),
          _ce_bptr(nullptr),
          _mut_ptr_cont(),
          _bits(0)
    {
        _set_rc(false);
        _set_is_unmappable(false);
    }

    /** Constructs a single read chunk from a read, and maps it to a new contig.
     * @param re_ptr Pointer to read entry where this chunk comes from.
     * @param len Length of the read.
     */
    Read_Chunk(Read_Entry_BPtr re_bptr, Size_Type len)
        : _r_start(0),
          _r_len(len),
          _c_start(0),
          _c_len(len),
          _re_bptr(re_bptr),
          _ce_bptr(NULL),
          _mut_ptr_cont(),
          _bits(0)
    {
        _set_rc(false);
        _set_is_unmappable(false);
    }

    /** Construct chunk mapping full read to single contig.
     * Pre: Read and contig must be of the same length.
     * @param re_bptr Read entry.
     * @param ce_bptr Contig entry.
     */
    Read_Chunk(Read_Entry_BPtr re_bptr, Contig_Entry_BPtr ce_bptr);
    /**@}*/

    // allow move only when unlinked
    DELETE_COPY_CTOR(Read_Chunk)
    Read_Chunk(Read_Chunk&& rhs) { *this = std::move(rhs); }

    // when deallocating, mut_ptr_cont must have been cleared
    ~Read_Chunk()
    {
        ASSERT(mut_ptr_cont().empty());
        ASSERT(is_unlinked());
    }

public:
    DELETE_COPY_ASOP(Read_Chunk)
    Read_Chunk& operator = (Read_Chunk&& rhs)
    {
        if (this != &rhs)
        {
            ASSERT(is_unlinked() and rhs.is_unlinked());
            _re_bptr = std::move(rhs._re_bptr);
            _ce_bptr = std::move(rhs._ce_bptr);
            _mut_ptr_cont = std::move(rhs._mut_ptr_cont);
            _r_start = std::move(rhs._r_start);
            _r_len = std::move(rhs._r_len);
            _c_start = std::move(rhs._c_start);
            _c_len = std::move(rhs._c_len);
            _bits = std::move(rhs._bits);
        }
        return *this;
    }

    /** @name Getters */
    /**@{*/
    const Read_Entry_BPtr& re_bptr() const { return _re_bptr; }
    Read_Entry_BPtr& re_bptr() { return _re_bptr; }
    const Contig_Entry_BPtr& ce_bptr() const { return _ce_bptr; }
    Contig_Entry_BPtr& ce_bptr() { return _ce_bptr; }
    Size_Type get_r_start() const { return _r_start; }
    Size_Type get_r_len() const { return _r_len; }
    Size_Type get_r_end() const { return _r_start + _r_len; }
    Size_Type get_c_start() const { return _c_start; }
    Size_Type get_c_len() const { return _c_len; }
    Size_Type& c_len() { return _c_len; }
    Size_Type get_c_end() const { return _c_start + _c_len; }
    bool get_rc() const { return _get_rc(); }
    const Mutation_Ptr_Cont& mut_ptr_cont() const { return _mut_ptr_cont; }
    Mutation_Ptr_Cont& mut_ptr_cont() { return _mut_ptr_cont; }
    Seq_Type get_seq() const;
    Seq_Type substr(Size_Type start, Size_Type len) const;
    bool is_unmappable() const { return _get_is_unmappable(); }
    Size_Type get_read_len() const;
    /**@}*/

    /** Set of coordinates used to traverse a Read_Chunk mapping. */
    /** Get mapping start position. */
    Pos get_start_pos() const
    {
        return Pos(get_c_start(), (not get_rc()? get_r_start() : get_r_end()), _mut_ptr_cont.begin(), 0, this);
    }

    /** Get mapping end position. */
    Pos get_end_pos() const
    {
        return Pos(get_c_end(), (not get_rc()? get_r_end() : get_r_start()), _mut_ptr_cont.end(), 0, this);
    }

    /** Split Read_Chunk based on contig position.
     * Pre: No Mutation may span c_brk.
     * Pre: Read_Chunk must be unlinked from all containers.
     * Post: If a split is not required, one component is the original pointer, and the other is NULL.
     * If a split is required, a new Read_Chunk is allocated, and backpointers in Mutation_Ptr_Cont are updated.
     * @param rc_bptr Chunk to split.
     * @param c_brk Contig breakpoint.
     * @param mut_left_cbptr Insertion at c_pos to remain on the left of the cut, if any.
     * @return Tuple of Read_Chunk pointers that go on left&right of the cut.
     */
    static std::tuple< Read_Chunk_BPtr, Read_Chunk_BPtr >
    split(Read_Chunk_BPtr rc_bptr, Size_Type c_brk, Mutation_CBPtr mut_left_cbptr);

    /** Split read chunk based on contig position and mutation map.
     * @param c_brk Contig breakpoint.
     * @param mut_cptr_map Map of mutations which go to the right side.
     * @param ce_cptr Pointer to new contig entry object.
     * @return Flag whether to move read chunk to the right contig side; additional Read_Chunk object mapped to right side, if the chunk gets split
     */
    //std::tuple< bool, std::shared_ptr< Read_Chunk > > split(
    //    Size_Type c_brk, const std::map< const Mutation*, const Mutation* >& mut_cptr_map, const Contig_Entry* ce_cptr);

    /** Create Read_Chunk object from a cigar string.
     * Also creates a corresponding Contig_Entry object.
     * @param cigar Cigar string.
     * @param rf_ptr Reference string (either whole, or just mapped portion); Contig_Entry object takes ownership.
     * @param qr Query string (either whole, or just mapped portion).
     */
    static Read_Chunk_BPtr
    make_chunk_from_cigar(const Cigar& cigar, Seq_Type&& rf_ptr, const Seq_Type& qr);

    /** Create Read_Chunk object from a cigar string.
     * This version assumes a Contig_Entry object exists holding the cigar rf sequence.
     * @param cigar Cigar string.
     * @param qr Query string (either whole, or just mapped portion).
     * @param ce_bptr Contig_Entry object holding the cigar rf sequence.
     */
    static Read_Chunk_BPtr
    make_chunk_from_cigar(const Cigar& cigar, const Seq_Type& qr, Contig_Entry_BPtr ce_bptr);

    /** Create Read_Chunk and Contig_Entry objects from a cigar string and 2 existing Read_Chunk objects.
     * @param rc1_cbptr Read_Chunk corresponding to rf.
     * @param rc2_cbptr Read_Chunk corresponding to qr.
     * @param cigar Cigar string.
     * @return A read chunk corresponding to the query in the cigar object, and a Contig_Entry corresponding to the rf.
     */
    static Read_Chunk_BPtr
    make_relative_chunk(Read_Chunk_CBPtr rc1_cbptr, Read_Chunk_CBPtr rc2_cbptr,
                        const Cigar& cigar);

    /** Collapse mappings: from c1<-rc1 and rc1<-rc2, compute c1<-rc2.
     * Pre: Entire rc1 part where rc2 is mapped must in turn be mapped to c1.
     * Post: New Mutation objects are created and linked;
     * if equivalent Mutations exist in the container, they are used instead.
     * Post: If an rc2 endpoint is mapped to an rc1 endpoint,
     * no c1 deletions will be omitted when the former is mapped to c1.
     * @param c1rc1_cbptr Mapping of rc1 on c1.
     * @param rc1rc2_cbptr Mapping of rc2 on rc1.
     * @param mut_cont Container of c1 mutations.
     * @return New Read_Chunk corresponding to mapping of rc2 on c1.
     */
    static Read_Chunk_BPtr
    collapse_mapping(Read_Chunk_BPtr c1rc1_cbptr, Read_Chunk_BPtr rc1rc2_cbptr,
                     Mutation_Cont& mut_cont);

    /** Invert mapping; from contig<-chunk, compute chunk<-contig.
     * Post: New Read_Chunk and Contig_Entry objects are allocated.
     * Post: Contig_Entry only has the sequence in the original chunk;
     * several methods are unavailable in this case.
     * @param rc_bptr Chunk to invert.
     */
    static Read_Chunk_BPtr
    invert_mapping(Read_Chunk_CBPtr rc_cbptr);

    /** Reverse the contig mapping; Mutations and their container must be reversed externally. */
    void reverse();

    /** Merge read chunk with the next chunk of the same read along the contig.
     * Pre: Chunks must be mapped to the same contig, in the same orientation, in order, continuously.
     * Pre: Chunks must be unlinked from their RE&CE containers.
     * Post: If the chunks have touching mutations at the breakpoint,
     * a new merged Mutation is used; if there is an equivalent Mutation in the continer,
     * it is used; otherwise a new merged Mutation is created and added to the container;
     * if either of the old mutations are no longer used, they are removed and deallcated.
     * Post: rc_next_bptr is deallcated.
     * @param rc_bptr First chunk.
     * @param rc_next_bptr Next chunk in read direction.
     * @param mut_cont Mutation continer to alter.
     */
    static void cat_c_right(Read_Chunk_BPtr rc_bptr, Read_Chunk_BPtr rc_next_bptr, Mutation_Cont& mut_cont);

    /** Shift contig coordinates of this Read_Chunk object.
     * @param delta Signed integer value to add to contig start point.
     */
    void shift(int delta)
    {
        ASSERT(int(_c_start) + delta >= 0);
        _c_start = Size_Type(int(_c_start) + delta);
    }

    /** Compute mapped range.
     * Given a read/contig range, compute the corresponding contig/read range
     * under the mapping described in this object.
     * Pre: rg_start <= rg_end.
     * Pre: the given range must be included in this chunk.
     * NOTE: "Maximal" means the returned range is maximal at that endpoint;
     * e.g. if the chunk is mapped to the positive strand (not rc) of the contig,
     * and the original range is on the contig,
     * then rg_start_maximal == true means the left endpoint of the returned read range
     * is as _small_ as possible.
     * NOTE: It is possible to get a negative range returned (meaning start > end)
     * iff the start range is empty and rg_start_maximal and rg_end_maximal are both false.
     * @param rg_start Start of original range.
     * @param rg_end End of original range.
     * @param on_contig True iff original range is on contig.
     * @param rg_start_maximal True: use maximal mapping for start endpoint of original range.
     * @param rg_end_maximal True: use maximal mapping for end endpoint of original range.
     */
    std::tuple< Size_Type, Size_Type >
    mapped_range(Size_Type rg_start, Size_Type rg_end, bool on_contig,
                 bool rg_start_maximal, bool rg_end_maximal) const;

    /** Make chunk unmappable.
     * Creates a new contig entry containing only this chunk mapped to the positive strand.
     * Pre: Chunk should be removed from its previous contig entry container (call during clear and destroy).
     * NOTE: New contig entry is not inserted in the graph.
     * @param rc_bptr Chunk to make unmappable.
     */
    static void make_unmappable(Read_Chunk_BPtr rc_bptr);

    /** Integrity check. */
    void check() const;

    //friend std::ostream& operator << (std::ostream&, const Read_Chunk&);
    boost::property_tree::ptree to_ptree() const;
    static std::string to_string(Read_Chunk_CBPtr rc_cbptr, bool r_dir = true, bool forward = true);

private:
    /** Hooks for storage:
     * 1. in intrusive interval trees inside Contig_Entry objects.
     * 2. in intrusive trees inside Read_Entry objects.
     */
    friend struct detail::Read_Chunk_ITree_Node_Traits;
    friend struct detail::Read_Chunk_Set_Node_Traits;
    bool is_unlinked() const { return not(_ce_parent or _ce_l_child or _ce_r_child or _re_parent or _re_l_child or _re_r_child); }

    Size_Type _r_start;
    Size_Type _r_len;
    Size_Type _c_start;
    Size_Type _c_len;
    Size_Type _ce_max_end;

    Read_Entry_BPtr _re_bptr;
    Contig_Entry_BPtr _ce_bptr;
    Mutation_Ptr_Cont _mut_ptr_cont;
    Read_Chunk_BPtr _ce_parent;
    Read_Chunk_BPtr _ce_l_child;
    Read_Chunk_BPtr _ce_r_child;
    Read_Chunk_BPtr _re_parent;
    Read_Chunk_BPtr _re_l_child;
    Read_Chunk_BPtr _re_r_child;

    uint8_t _bits;

    static const uint8_t _rc_bit = 1u;
    static const uint8_t _is_unmappable_bit = 2u;
    static const uint8_t _ce_col_bit = 4u;
    static const uint8_t _re_col_bit = 8u;

    bool _get_rc() const { return bitmask::any(_bits, _rc_bit); }
    void _set_rc(bool v) { bitmask::assign(_bits, _rc_bit, v); }
    bool _get_is_unmappable() const { return bitmask::any(_bits, _is_unmappable_bit); }
    void _set_is_unmappable(bool v) { bitmask::assign(_bits, _is_unmappable_bit, v); }
    bool _get_ce_col() const { return bitmask::any(_bits, _ce_col_bit); }
    void _set_ce_col(bool v) { bitmask::assign(_bits, _ce_col_bit, v); }
    bool _get_re_col() const { return bitmask::any(_bits, _re_col_bit); }
    void _set_re_col(bool v) { bitmask::assign(_bits, _re_col_bit, v); }
}; // class Read_Chunk

} // namespace MAC


#endif
