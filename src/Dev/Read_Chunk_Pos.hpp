//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __READ_CHUNK_POS_HPP
#define __READ_CHUNK_POS_HPP

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
    /** Iterator with next mutation pointer. */
    Mutation_Cont::const_iterator mut_cit;
    /** Allele index for next mutation. */
    Allele_Idx_Cont::const_iterator allele_idx_cit;
    /** If non-zero, offset into current mutation. */
    Size_Type mut_offset;
private:
    /** Read_Chunk parent object. */
    const Read_Chunk* rc_cptr;
    Mutation_Cont::const_iterator mut_it_begin;
    Mutation_Cont::const_iterator mut_it_end;

private:
    // Constructed only by Read_Chunk objects
    friend class Read_Chunk;
    /** Constructor. */
    Read_Chunk_Pos(Size_Type _c_pos,
                   Size_Type _r_pos,
                   Mutation_Cont::const_iterator _mut_cit,
                   Allele_Idx_Cont::const_iterator _allele_idx_cit,
                   Size_Type _mut_offset,
                   const Read_Chunk* _rc_cptr)
        : c_pos(_c_pos),
          r_pos(_r_pos),
          mut_cit(_mut_cit),
          allele_idx_cit(_allele_idx_cit),
          mut_offset(_mut_offset),
          rc_cptr(_rc_cptr)
    {
        mut_it_begin = rc_cptr->mut_it_begin();
        mut_it_end = rc_cptr->mut_it_end();
    }

public:
    // allow copy and move
    DEFAULT_COPY_CTOR(Read_Chunk_Pos);
    DEFAULT_MOVE_CTOR(Read_Chunk_Pos);
    DEFAULT_COPY_ASOP(Read_Chunk_Pos);
    DEFAULT_MOVE_ASOP(Read_Chunk_Pos);

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
        return *mut_cit;
    }

    /** Get previous mutation. */
    const Mutation& prev_mut() const
    {
        ASSERT(rc_cptr);
        ASSERT(past_first_mut());
        auto tmp_cit = mut_cit;
        return *(--tmp_cit);
    }

    /** Get length of unmutated mapped stretch starting at given position.
     * @param forward Bool; true if looking for match stretch forward, false if backward.
     * @return The length of the next stretch of unmutated bases; 0 if on the breakpoint of a mutation.
     */
    Size_Type get_match_len(bool forward = true) const;

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

    friend bool operator == (const Read_Chunk_Pos& lhs, const Read_Chunk_Pos& rhs)
    {
        return (lhs.c_pos == rhs.c_pos
                and lhs.r_pos == rhs.r_pos
                and lhs.mut_cit == rhs.mut_cit
                and lhs.mut_offset == rhs.mut_offset);
    }
    friend bool operator != (const Read_Chunk_Pos& lhs, const Read_Chunk_Pos& rhs) { return !(lhs == rhs); }
}; // struct Read_Chunk_Pos

} // namespace MAC


#endif
