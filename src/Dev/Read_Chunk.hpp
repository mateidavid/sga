//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __READ_CHUNK_HPP
#define __READ_CHUNK_HPP

#include "MAC_forward.hpp"
#include "Mutation_Cont.hpp"
#include "Allele_Idx_Cont.hpp"
#include "Read_Sub_Chunk.hpp"


namespace MAC
{

/**
 * Holds information about a read chunk.
 * Reads will be partitioned into chunks of non-zero size as they are mapped to contigs.
 * Each read chunk object contains information about the read it comes from (start, span),
 * the contig it is mapped to (start, span, orientation), and the contig mutations it observes.
 */
class Read_Chunk
{
private:
    // Can only be constructed by Factory object
    friend class bounded::Factory< Read_Chunk >;

    /// Default constructor.
    Read_Chunk()
        : _r_start(0),
          _r_len(0),
          _c_start(0),
          _c_len(0),
          _bits(0)
    {}

    // disallow copy or move
    DELETE_COPY_CTOR(Read_Chunk);
    DELETE_MOVE_CTOR(Read_Chunk);
    DELETE_COPY_ASOP(Read_Chunk);
    DELETE_MOVE_ASOP(Read_Chunk);

    /**
     * Constructor that maps a full read to a single contig.
     * Pre: Read and contig must be of the same length.
     * @param re_bptr Read_Entry object.
     * @param ce_bptr Contig_Entry object.
     */
    Read_Chunk(Read_Entry_BPtr re_bptr, Contig_Entry_BPtr ce_bptr);

    /// Constructor from mapping positions.
    Read_Chunk(Size_Type r_start, Size_Type r_len,
               Size_Type c_start, Size_Type c_len,
               bool rc)
        : _r_start(r_start),
          _r_len(r_len),
          _c_start(c_start),
          _c_len(c_len),
          _bits(0)
    {
        _set_rc(rc);
    }

    ~Read_Chunk()
    {
        ASSERT(is_unlinked());
        _allele_idx_cont.clear_and_dispose();
    }

public:
    /// @name Getters
    /// @{
    GETTER(Read_Entry_BPtr, re_bptr, _re_bptr)
    GETTER(Contig_Entry_BPtr, ce_bptr, _ce_bptr)
    GETTER(Size_Type, r_start, _r_start)
    GETTER(Size_Type, r_len, _r_len)
    GETTER(Size_Type, c_len, _c_len)
    Size_Type get_r_start() const { return _r_start; }
    Size_Type get_r_len() const { return _r_len; }
    Size_Type get_r_end() const { return _r_start + _r_len; }
    Size_Type get_c_start() const { return _c_start; }
    Size_Type get_c_len() const { return _c_len; }
    Size_Type get_c_end() const { return _c_start + _c_len; }
    bool get_rc() const { return _get_rc(); }
    Seq_Type get_seq() const;
    bool is_unbreakable() const { return _get_is_unbreakable(); }
    Size_Type get_read_len() const;
    Mutation_Cont::const_iterator mut_it_begin() const;
    Mutation_Cont::const_iterator mut_it_end() const;
    /// @}

    /// Get first subchunk
    Read_Sub_Chunk get_first_sub_chunk() const
    {
        Read_Sub_Chunk res{get_c_start(), 0, get_r_start(), 0, 0, 0, get_c_start(), get_c_len(),
                mut_it_begin(), mut_it_begin(), mut_it_end(), _allele_idx_cont.begin(),
                0, 0, 0, 0,
                get_rc(), false};
        res.set_convenience_fields();
        return res;
    }

    /// Get last subchunk
    Read_Sub_Chunk get_last_sub_chunk() const
    {
        Read_Sub_Chunk res{get_c_end(), 0, get_r_end(), 0, 0, 0, get_c_start(), get_c_len(),
                mut_it_end(), mut_it_begin(), mut_it_end(), _allele_idx_cont.end(),
                0, 0, 0, 0,
                get_rc(), false};
        res.set_convenience_fields();
        return res;
    }

    /**
     * Split Read_Chunk based on contig position.
     * Pre: No Mutation may span c_brk.
     * Pre: Read_Chunk must be unlinked from all containers.
     * Post: If a split is not required, one component is the original pointer, and the other is NULL.
     * If a split is required, a new Read_Chunk is allocated, and backpointers in Mutation_Ptr_Cont are updated.
     * Whenever a split is required, the original Read_Chunk carries the LHS of the split.
     * @param rc_bptr Chunk to split.
     * @param c_brk Contig breakpoint.
     * @param mut_left_cbptr Insertion at c_pos to remain on the left of the cut, if any.
     * @param strict If true, every non-terminal Read_Chunk that spans the break gets split
     * into exactly 2 pieces (even if one is empty).
     * @return Tuple of Read_Chunk pointers that go on left&right of the cut.
     *
    static pair< Read_Chunk_BPtr, Read_Chunk_BPtr >
    split(Read_Chunk_BPtr rc_bptr, Size_Type c_brk, Mutation_CBPtr mut_left_cbptr, bool strict = false);
    */

    /**
     * Construct Read_Chunk object from a Cigar object.
     * Also creates a corresponding Contig_Entry object.
     * @param cigar Cigar string.
     * @param rf_ptr Reference string (either whole, or just mapped portion); Contig_Entry object takes ownership.
     * @param qr Query string (either whole, or just mapped portion).
     *
    static Read_Chunk_BPtr
    make_chunk_from_cigar(const Cigar& cigar, Seq_Type&& rf_ptr, const Seq_Proxy_Type& qr);
    */

    /**
     * Construct Read_Chunk object from a Cigar object, given existing Contig_Entry.
     * This version assumes a Contig_Entry object exists holding the cigar rf sequence.
     * @param cigar Cigar string.
     * @param qr Query string (either whole, or just mapped portion).
     * @param ce_bptr Contig_Entry object holding the cigar rf sequence.
     *
    static Read_Chunk_BPtr
    make_chunk_from_cigar(const Cigar& cigar, const Seq_Proxy_Type& qr, Contig_Entry_BPtr ce_bptr);
    */

    /**
     * Create Read_Chunk and Contig_Entry objects from a cigar string and 2 existing Read_Chunk objects.
     * @param rc1_cbptr Read_Chunk corresponding to rf.
     * @param rc2_cbptr Read_Chunk corresponding to qr.
     * @param cigar Cigar string.
     * @return A read chunk corresponding to the query in the cigar object, and a Contig_Entry corresponding to the rf.
     *
    static Read_Chunk_BPtr
    make_relative_chunk(Read_Chunk_CBPtr rc1_cbptr, Read_Chunk_CBPtr rc2_cbptr,
                        const Cigar& cigar);
    */

    /**
     * Collapse mappings: from c1<-rc1 and rc1<-rc2, compute c1<-rc2.
     * Pre: Entire rc1 part where rc2 is mapped must in turn be mapped to c1.
     * Post: New Mutation objects are created and linked;
     * if equivalent Mutations exist in the container, they are used instead.
     * Post: If an rc2 endpoint is mapped to an rc1 endpoint,
     * no c1 deletions will be omitted when the former is mapped to c1.
     * @param c1rc1_cbptr Mapping of rc1 on c1.
     * @param rc1rc2_cbptr Mapping of rc2 on rc1.
     * @param mut_cont Container of c1 mutations.
     * @return New Read_Chunk corresponding to mapping of rc2 on c1.
     *
    static Read_Chunk_BPtr
    collapse_mapping(Read_Chunk_BPtr c1rc1_cbptr, Read_Chunk_BPtr rc1rc2_cbptr,
                     Mutation_Cont& mut_cont);
    */

    /// Reverse the contig mapping; Mutations and their container must be reversed externally.
    //void reverse();

    /**
     * Merge read chunk with the next chunk of the same read along the contig.
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
    //static void cat_c_right(Read_Chunk_BPtr rc_bptr, Read_Chunk_BPtr rc_next_bptr, Mutation_Cont& mut_cont);

    /**
     * Shift contig coordinates of this Read_Chunk object.
     * @param delta Signed integer value to add to contig start point.
     */
    void shift(int delta)
    {
        ASSERT(int(_c_start) + delta >= 0);
        _c_start = Size_Type(int(_c_start) + delta);
    }

    /**
     * Compute mapped range.
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
    /*
    Range_Type mapped_range(Range_Type rg, bool on_contig,
                            bool rg_start_maximal, bool rg_end_maximal) const;
    */

    /**
     * Make chunk unmappable.
     * Creates a new contig entry containing only this chunk mapped to the positive strand.
     * Pre: Chunk should be removed from its previous contig entry container (call during clear and destroy).
     * NOTE: New contig entry is not inserted in the graph.
     * @param rc_bptr Chunk to make unmappable.
     */
    //static void make_unmappable(Read_Chunk_BPtr rc_bptr);

    /// Integrity check.
    void check() const;

    //friend ostream& operator << (ostream&, const Read_Chunk&);
    boost::property_tree::ptree to_ptree() const;
    static string to_string(Read_Chunk_CBPtr rc_cbptr, bool r_dir = true, bool forward = true);

private:
    friend struct detail::Read_Chunk_ITree_Node_Traits;
    friend struct detail::Read_Chunk_Set_Node_Traits;
    bool is_unlinked() const
    {
        return not(_ce_parent or _ce_l_child or _ce_r_child
                   or _re_parent or _re_l_child or _re_r_child);
    }

    Size_Type _r_start;
    Size_Type _r_len;
    Size_Type _c_start;
    Size_Type _c_len;
    Size_Type _ce_max_end;

    Read_Entry_BPtr _re_bptr;
    Contig_Entry_BPtr _ce_bptr;
    Allele_Idx_Cont _allele_idx_cont;

    Read_Chunk_BPtr _ce_parent;
    Read_Chunk_BPtr _ce_l_child;
    Read_Chunk_BPtr _ce_r_child;
    Read_Chunk_BPtr _re_parent;
    Read_Chunk_BPtr _re_l_child;
    Read_Chunk_BPtr _re_r_child;

    uint8_t _bits;

    static const uint8_t _rc_bit = 1u;
    static const uint8_t _is_unbreakable_bit = 2u;
    static const uint8_t _ce_col_bit = 4u;
    static const uint8_t _re_col_bit = 8u;

    bool _get_rc() const { return bitmask::any(_bits, _rc_bit); }
    void _set_rc(bool v) { bitmask::assign(_bits, _rc_bit, v); }
    bool _get_is_unbreakable() const { return bitmask::any(_bits, _is_unbreakable_bit); }
    void _set_is_unbreakable(bool v) { bitmask::assign(_bits, _is_unbreakable_bit, v); }
    bool _get_ce_col() const { return bitmask::any(_bits, _ce_col_bit); }
    void _set_ce_col(bool v) { bitmask::assign(_bits, _ce_col_bit, v); }
    bool _get_re_col() const { return bitmask::any(_bits, _re_col_bit); }
    void _set_re_col(bool v) { bitmask::assign(_bits, _re_col_bit, v); }
}; // class Read_Chunk

} // namespace MAC


#endif
