//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __CONTIG_ENTRY_HPP
#define __CONTIG_ENTRY_HPP

#include "MAC_forward.hpp"
#include "Mutation.hpp"
#include "Mutation_Cont.hpp"
#include "Read_Chunk.hpp"
#include "Read_Chunk_Cont.hpp"


namespace MAC
{

namespace detail
{

struct Contig_Entry_List_Node_Traits;

}

/** Holds information about a contig.
 *
 * Each object holds a base sequence, a list of observed mutations, and a list of pointers
 * to read chunks which are mapped to this contig.
 */
class Contig_Entry
{
private:
    // only constructed by Factory object
    friend class bounded::Factory< Contig_Entry >;

    /** Constructor.
     * @param seq RVR to read sequence (assume ownership).
     */
    Contig_Entry(Seq_Type&& seq, Size_Type seq_offset = 0) : _seq(std::move(seq)), _seq_offset(seq_offset), _tag(0) {}

    // allow move only when unlinked
    DEFAULT_DEF_CTOR(Contig_Entry)
    DELETE_COPY_CTOR(Contig_Entry)
    Contig_Entry(Contig_Entry&& rhs) { *this = std::move(rhs); }

    ~Contig_Entry()
    {
        ASSERT(chunk_cont().empty());
        ASSERT(mut_cont().empty());
        ASSERT(is_unlinked());
    }

public:
    DELETE_COPY_ASOP(Contig_Entry)
    Contig_Entry& operator = (Contig_Entry&& rhs)
    {
        if (this != &rhs)
        {
            ASSERT(is_unlinked());
            _seq = std::move(rhs._seq);
            _seq_offset = std::move(rhs._seq_offset);
            _tag = std::move(rhs._tag);
            _mut_cont = std::move(rhs._mut_cont);
            _chunk_cont = std::move(rhs._chunk_cont);
        }
        return *this;
    }

    /** @name Getters */
    /**@{*/
    GETTER(Seq_Type, seq, _seq)
    Size_Type get_seq_offset() const { return _seq_offset; }
    Seq_Type substr(Size_Type start, Size_Type len) const
    {
        ASSERT(start >= _seq_offset and start + len <= _seq_offset + _seq.size());
        return _seq.substr(start - _seq_offset, len);
    }
    Size_Type get_len() const { return _seq.size(); }
    GETTER(Mutation_Cont, mut_cont, _mut_cont)
    GETTER(Read_Chunk_CE_Cont, chunk_cont, _chunk_cont)
    GETTER(uint64_t, tag, _tag)
    bool old_is_unmappable() const { return _chunk_cont.size() == 1 and _chunk_cont.begin()->is_unmappable(); }
    bool is_unmappable() const { ASSERT(old_is_unmappable() == is_ambiguous()); return is_ambiguous(); }

    bool is_ambiguous() const { return bitmask::any(_tag, is_ambiguous_mask); }
    void set_ambiguous() { bitmask::set(_tag, is_ambiguous_mask); }
    void reset_ambiguous() { bitmask::reset(_tag, is_ambiguous_mask); }
    bool is_lowcomplex() const { return bitmask::any(_tag, is_lowcomplex_mask); }
    void set_lowcomplex() { bitmask::set(_tag, is_lowcomplex_mask); }
    void reset_lowcomplex() { bitmask::reset(_tag, is_lowcomplex_mask); }
    /**@}*/

    /** Cut mutation at given offsets.
     * Read chunks containing the original Mutation will get instead 2 adjacent mutations.
     * @param mut_bptr Mutation to cut.
     * @param c_offset Contig offset inside mutation where to cut.
     * @param r_offset Read offset inside mutation where to cut.
     */
    void cut_mutation(Mutation_BPtr mut_bptr, Size_Type c_offset, Size_Type r_offset);

    /** Reverse the contig. */
    void reverse();

    /** Retrieve chunks leaving contig.
     * For every chunk rc leaving the contig in the specified direction,
     * and whose sibling chunk rc' is unmappable, the following policies are available:
     * 0: skip such rc;
     * 1: ignore issue: each such rc will appear in a bucket on its own;
     * 2: collect all such rc in a special bucket with key (NULL, false);
     * 3: skip unmappable chunks: instead of rc', continue traversing the read
     * and use the next mappable chunk instead; discard rc if rc' is terminal and unmappable
     * 4: same as 3, but if rc' is terminal and unmappable, collect rc in special bucket (NULL, false)
     * @param c_right Bool; true: spanning past contig end, false: spanning past contig start.
     * @param unmappable_policy See above.
     * @param ignore_threshold Ignore contigs&orientations supported by that many
     * chunks or less (0 means nothing is ignored).
     * @return A map containing, for each Contig_Entry and orientation,
     * a list of Read_Chunk objects whose sibling is in that Contig_Entry and has the same orientation;
     * two chunks have the same orientation if they are both mapped to the positive or negative
     * strand of their respective contigs.
     */
    typedef map< tuple< Contig_Entry_CBPtr, bool >, vector< Read_Chunk_CBPtr > > out_chunks_dir_type;
    out_chunks_dir_type out_chunks_dir(bool c_right, int unmappable_policy, size_t ignore_threshold = 0) const;

    /** Get out-edges counts.
     * @param unmappable_policy See above.
     * @param ignore_threshold See above.
     * @return A tuple (cnt_left, uniq_left, cnt_right, uniq_right), where cnt is the number
     * of read chunks spanning that breakpoint, and uniq is the number of different contig entries
     * where following chunks are mapped.
     */
    std::tuple< size_t, size_t, size_t, size_t >
    get_out_degrees(int unmappable_policy, size_t ignore_threshold = 0) const;

    /** Compute the min and max unmappable regions neighbouring this contig.
     * Given a set of chunks leaving the contig in the given direction,
     * compute the minimum and maximum unmappable sibling chunks.
     * Pre: Chunks must exit the contig in the given direction.
     * Pre: After skipping unmappable chunks, each chunk must eventually
     * have a sibling mappable chunk before the end of the read.
     * NOTE: This method does not assume the chunks eventually lead
     * into the same neighbour contig.
     * NOTE: This method is only be useful when chunk_cont is one of the
     * vectors returned by out_chunks_dir() under unmappable policy 3.
     * @param c_right Bool; true: spanning past contig end, false: spanning past contig start.
     * @param chunk_cont Container of chunks leaving the contig in direction c_right.
     * @return (Minimum, maximum) unmappable regions.
     */
    std::tuple< Size_Type, Size_Type >
    unmappable_neighbour_range(bool c_right, const std::vector< Read_Chunk_CBPtr >& chunk_cont) const;

    /** Check if contig is catenable with a single other contig in the given direction.
     * @param c_right Bool; true: merge past contig end, false: merge past contig start.
     * @return A tuple (Contig_Entry, same_orientation, Read_Chunk list);
     * if the chunk is not catenable, the first coordinate is NULL.
     */
    std::tuple< Contig_Entry_CBPtr, bool, std::vector< Read_Chunk_CBPtr > >
    can_cat_dir(bool c_right) const;

    /** Merge two contigs in the first.
     * Pre: The contig must be in the same orientation as this one.
     * Pre: The second contig must be unlinked.
     * Post: The second contig is deallocated.
     * @param ce_cptr Initial contig.
     * @param ce_next_cptr Next contig.
     * @param rc_cbptr_cont Read chunks from first contig that span into the next contig.
     */
    static void cat_c_right(Contig_Entry_BPtr ce_bptr, Contig_Entry_BPtr ce_next_bptr,
                            std::vector< Read_Chunk_CBPtr >& rc_cbptr_cont);

    /** Get neighbour contigs.
     * @param dir Bool; true: past contig end, false: past contig start.
     * @param skip_unmappable Bool; true: include mappable contig after skipping an unmappable chunk;
     * false: include only mappable contigs that are direct neighbours.
     */
    //std::shared_ptr< std::map< std::tuple< const Contig_Entry*, bool >, std::tuple< unsigned int, Size_Type, Size_Type > > >
    //get_neighbours(bool dir, bool skip_unmappable = true, bool trim_results = true) const;

    /** Get heterozygous mutations isolated from other nearby mutations.
     * @param min_support_report Minimum read chunk support for mutation and for base.
     * @param min_separation Minimum distance to a nearby mutation or contig end.
     */
    //std::vector< Mutation_CPtr > get_separated_het_mutations(
    //    size_t min_support_report, Size_Type min_separation) const;

    //void print_separated_het_mutations(std::ostream& os,
    //                                   size_t min_support_report, Size_Type min_separation) const;

    static void dispose(Contig_Entry_BPtr ce_bptr)
    {
        // remove chunks from their RE cont
        ce_bptr->chunk_cont().erase_from_re_cont();
        // deallocate chunks and mca-s
        ce_bptr->chunk_cont().clear_and_dispose();
        // deallocate mutations
        ce_bptr->mut_cont().clear_and_dispose();
        // deallocate CE
        Contig_Entry_Fact::del_elem(ce_bptr);
    }

    /** Integrity check. */
    void check() const;

    //friend std::ostream& operator << (std::ostream& os, const Contig_Entry& rhs);
    boost::property_tree::ptree to_ptree() const;

    static const uint64_t is_ambiguous_mask = 1u << 31;
    static const uint64_t is_lowcomplex_mask = 1u << 30;

private:
    friend struct detail::Contig_Entry_List_Node_Traits;
    bool is_unlinked() const { return not(_previous or _next); }

    Seq_Type _seq;
    Size_Type _seq_offset;
    uint64_t _tag;
    Mutation_Cont _mut_cont;
    Read_Chunk_CE_Cont _chunk_cont;

    /** Hooks for storage in intrusive list in Graph object. */
    Contig_Entry_BPtr _previous;
    Contig_Entry_BPtr _next;
}; // class Contig_Entry

} // namespace MAC


#endif
