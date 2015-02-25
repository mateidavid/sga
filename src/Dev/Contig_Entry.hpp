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
    Contig_Entry(Seq_Type&& seq, Size_Type seq_offset = 0) : _seq(move(seq)), _seq_offset(seq_offset), _tag(0) {}

    // allow move only when unlinked
    Contig_Entry() : _seq_offset(0), _tag(0) {}
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
    Size_Type seq_offset() const { return _seq_offset; }
    Seq_Proxy_Type substr(Size_Type start, Size_Type len) const
    {
        ASSERT(start >= _seq_offset and start + len <= _seq_offset + _seq.size());
        return seq().substr(start - _seq_offset, len);
    }
    Size_Type len() const { return _seq.size(); }
    GETTER(Mutation_Cont, mut_cont, _mut_cont)
    GETTER(Read_Chunk_CE_Cont, chunk_cont, _chunk_cont)
    GETTER(uint64_t, tag, _tag)
    bool is_unmappable() const { return bitmask::any(_tag, is_unmappable_mask); }
    void set_unmappable() { bitmask::set(_tag, is_unmappable_mask); }
    void reset_unmappable() { bitmask::reset(_tag, is_unmappable_mask); }
    bool is_lowcomplex() const { return bitmask::any(_tag, is_lowcomplex_mask); }
    void set_lowcomplex() { bitmask::set(_tag, is_lowcomplex_mask); }
    void reset_lowcomplex() { bitmask::reset(_tag, is_lowcomplex_mask); }
    bool is_normal() const { return not bitmask::any(_tag, is_unmappable_mask | is_lowcomplex_mask); }
    typedef size_t category_type;
    category_type category() const { return _tag & (is_unmappable_mask | is_lowcomplex_mask); }
    /**@}*/

    /** Cut mutation at given offsets.
     * Read chunks containing the original Mutation will get instead 2 adjacent mutations.
     * @param mut_bptr Mutation to cut.
     * @param c_offset Contig offset inside mutation where to cut.
     * @param r_offset Read offset inside mutation where to cut.
     */
    void cut_mutation(Mutation_BPtr mut_bptr, Size_Type c_offset, Size_Type r_offset);

    /** Cut Contig_Entry object.
     * @param c_brk Contig position of the cut.
     * @param mut_left_cbptr If not NULL, single insertion at breakpoint
     * to be kept on the LHS (others are moved to RHS).
     * @return RHS contig entry if one is created, nullptr otherwise.
     * NOTE: If a new RHS contig entry is created, it is not added to the contig entry container.
     */
    Contig_Entry_BPtr cut(Size_Type c_brk, Mutation_CBPtr mut_left_cbptr);

    /** Get chunks supporting qr and rf alleles of a mutation.
     * @param mut_cbptr Mutation of interest.
     * @return A tuple of sets; first = chunks supporting qr allele;
     * second = chunks supporting rf allele and fully spanning it;
     * third = chunks supporting rf allele but only partially spanning it.
     * NOTE: For a chunk to support the rf allele and fully span it:
     * if the rf allele is non-empty, the chunk must start at or before the  mutation start
     * and end at or after the mutation end;
     * if the rf allele is empty (length 0, an insertion), we further require >=1bp
     * to be mapped on either side of the mutation.
     */
    tuple< set< Read_Chunk_CBPtr >, set< Read_Chunk_CBPtr >, set< Read_Chunk_CBPtr > >
    mut_support(Mutation_CBPtr mut_cbptr) const;

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
    typedef map< pair< Contig_Entry_CBPtr, bool >, set< Read_Chunk_CBPtr > > out_chunks_dir_type;
    out_chunks_dir_type out_chunks_dir(bool c_right, int unmappable_policy, size_t ignore_threshold = 0) const;

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
    pair< Size_Type, Size_Type >
    unmappable_neighbour_range(bool c_right, const set< Read_Chunk_CBPtr >& chunk_cont) const;

    /** Options for neighbour traversal. */
    class Neighbour_Options
    {
    public:
        typedef size_t option_type;
        static const option_type opt_none = 0;
        static const option_type opt_drop = 1;
        static const option_type opt_skip = 2;
        static const option_type opt_skip_not_last = 3;
        Neighbour_Options& none(category_type t) { _map[t] = opt_none; return *this; }
        Neighbour_Options& drop(category_type t) { _map[t] = opt_drop; return *this; }
        Neighbour_Options& skip(category_type t) { _map[t] = opt_skip; return *this; }
        Neighbour_Options& skip_not_last(category_type t) { _map[t] = opt_skip_not_last; return *this; }
        Neighbour_Options& clear() { _map.clear(); return *this; }
        Neighbour_Options& default_opt(option_type o) { _default_opt = o; return *this; }
        option_type get_opt(category_type t) const { return _map.count(t) > 0? _map.at(t) : _default_opt; }
    private:
        map< category_type, option_type > _map;
        option_type _default_opt = opt_none;
    }; // Neighbour_Options

    /** Get neighbours of this Contig Entry. */
    typedef map< tuple< Contig_Entry_CBPtr, bool >, vector< tuple< Read_Chunk_CBPtr, Read_Chunk_CBPtr > > > neighbours_type;
    neighbours_type neighbours(bool forward, Neighbour_Options to, size_t ignore_threshold = 0);

    /** Check if contig is catenable with a single other contig in the given direction.
     * For this to be true, the contigs must have out-degree 1 in each-other's direction.
     * NOTE: This function does not check the contig's types are catenatable.
     * @param c_right Bool; true: merge past contig end, false: merge past contig start.
     * @return A tuple (Contig_Entry, same_orientation, Read_Chunk list);
     * if the chunk is not catenable, the first coordinate is NULL.
     */
    tuple< Contig_Entry_CBPtr, bool, set< Read_Chunk_CBPtr > >
    can_cat_dir(bool c_right) const;

    /** Merge two contigs in the first.
     * Pre: The contig must be in the same orientation as this one.
     * Pre: The second contig must be unlinked.
     * Post: The second contig is deallocated.
     * @param ce_bptr Initial contig.
     * @param ce_next_bptr Next contig.
     * @param rc_cbptr_cont Read chunks from first contig that span into the next contig.
     */
    static void cat_c_right(Contig_Entry_BPtr ce_bptr, Contig_Entry_BPtr ce_next_bptr,
                            set< Read_Chunk_CBPtr >& rc_cbptr_cont);

    /** Get neighbour contigs.
     * @param dir Bool; true: past contig end, false: past contig start.
     * @param skip_unmappable Bool; true: include mappable contig after skipping an unmappable chunk;
     * false: include only mappable contigs that are direct neighbours.
     */
    //std::shared_ptr< std::map< std::tuple< const Contig_Entry*, bool >, std::tuple< unsigned int, Size_Type, Size_Type > > >
    //get_neighbours(bool dir, bool skip_unmappable = true, bool trim_results = true) const;

    /** Check the Mutation objects in this Contig_Entry are well separated.
     * @param min_separation Minimum separation required.
     * @param include_endpoints Also check separation between first&last mutations and the endpoints
     * @return Bool; true iff all obejcts are well separated.
     */
    bool separated_mutations(Size_Type min_separation, bool include_endpoints) const
    {
        if (not mut_cont().empty())
        {
            Mutation_Cont::const_iterator it_last = mut_cont().begin();
            Mutation_Cont::const_iterator it = next(it_last);
            while (it != mut_cont().end())
            {
                if (it->rf_start() < it_last->rf_end() + min_separation) return false;
                it_last = it;
                ++it;
            }
            if (include_endpoints)
            {
                return min_separation <= mut_cont().begin()->rf_start()
                    and mut_cont().rbegin()->rf_end() + min_separation <= len();
            }
        }
        else // mut_cont is empty
        {
            if (include_endpoints)
            {
                return min_separation <= len();
            }
        }
        return true;
    }

    /// Swap reference and query alleles at a mutation.
    /// @param mut_cbptr Mutation whose alleles are to be swapped.
    /// @return New Mutation object.
    static Mutation_CBPtr swap_mutation_alleles(Contig_Entry_BPtr ce_bptr, Mutation_CBPtr mut_cbptr);

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

    static const uint64_t is_unmappable_mask = 1u << 31;
    static const uint64_t is_lowcomplex_mask = 1u << 30;
    static uint64_t supercontig_endpoint_mask(int i)
    {
        static uint64_t _supercontig_endpoint_mask[2] = { 1u << 29, 1u << 28 };
        return _supercontig_endpoint_mask[i];
    }
    static uint64_t visit_strand_mask(int i)
    {
        static uint64_t _visit_strand_mask[2] = { 1u << 27, 1u << 26 };
        return _visit_strand_mask[i];
    }

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
