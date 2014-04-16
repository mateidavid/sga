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

#include "MAC_forward.hpp"
#include "Mutation.hpp"
#include "Mutation_Cont.hpp"
#include "Read_Chunk.hpp"
#include "Read_Chunk_Cont.hpp"


namespace MAC
{

struct Contig_Entry_List_Node_Traits;
/** Holds information about a contig.
 *
 * Each object holds a base sequence, a list of observed mutations, and a list of pointers
 * to read chunks which are mapped to this contig.
 */
class Contig_Entry
{
private:
    // only constructed by Factory object
    friend class Factory< Contig_Entry >;

    /** Constructor.
     * @param seq RVR to read sequence (assume ownership).
     */
    Contig_Entry(Seq_Type&& seq, Size_Type seq_offset = 0) : _seq(std::move(seq)), _seq_offset(seq_offset), _colour(0) {}

    // allow move only when unlinked
    DEFAULT_DEF_CTOR(Contig_Entry)
    DELETE_COPY_CTOR(Contig_Entry)
    Contig_Entry(Contig_Entry&& rhs) { *this = std::move(rhs); }
    ~Contig_Entry() { ASSERT(is_unlinked()); }
public:
    DELETE_COPY_ASOP(Contig_Entry)
    Contig_Entry& operator = (Contig_Entry&& rhs)
    {
        if (this != &rhs)
        {
            ASSERT(is_unlinked());
            _seq = std::move(rhs._seq);
            _mut_cont = std::move(rhs._mut_cont);
            _chunk_cont = std::move(rhs._chunk_cont);
            _seq_offset = std::move(rhs._seq_offset);
            _contig_id = std::move(rhs._contig_id);
            _colour = std::move(rhs._colour);
        }
        return *this;
    }

    /** @name Getters */
    /**@{*/
    const Seq_Type& seq() const { return _seq; }
    Seq_Type& seq() { return _seq; }
    Size_Type get_seq_offset() const { return _seq_offset; }
    Seq_Type substr(Size_Type start, Size_Type len) const
    {
        ASSERT(start >= _seq_offset and start + len <= _seq_offset + _seq.size());
        return _seq.substr(start - _seq_offset, len);
    }
    Size_Type get_len() const { return _seq.size(); }
    const Mutation_Cont& mut_cont() const { return _mut_cont; }
    Mutation_Cont& mut_cont() { return _mut_cont; }
    const Read_Chunk_CE_Cont& chunk_cont() const { return _chunk_cont; }
    Read_Chunk_CE_Cont& chunk_cont() { return _chunk_cont; }
    const int& colour() const { return _colour; }
    int& colour() { return _colour; }
    const size_t& contig_id() const { return _contig_id; }
    size_t& contig_id() { return _contig_id; }
    bool is_unmappable() const { return _chunk_cont.size() == 1 and _chunk_cont.begin()->is_unmappable(); }
    /**@}*/

    /** Find bounded pointer to this object.
     * Pre: Must be linked.
     */
    Contig_Entry_CBPtr bptr_to() const
    {
        ASSERT(not is_unlinked());
        return _next->_previous;
    }
    Contig_Entry_BPtr bptr_to()
    {
        return boost::intrusive::pointer_traits< Contig_Entry_BPtr >::const_cast_from(
            const_cast< const Contig_Entry* >(this)->bptr_to());
    }

    /** Cut mutation at given offsets.
     * Read chunks containing the original Mutation will get instead 2 adjacent mutations.
     * @param mut_bptr Mutation to cut.
     * @param c_offset Contig offset inside mutation where to cut.
     * @param r_offset Read offset inside mutation where to cut.
     */
    void cut_mutation(Mutation_BPtr mut_bptr, Size_Type c_offset, Size_Type r_offset);

    /** Add mutation to container if it doesn't already exists.
     * @param m Mutation to look for.
     * @return Pointer to mutation.
     */
    //Mutation_CPtr add_mutation(const Mutation& m) { return add_mut_to_cont(_mut_cont, m); }

    /** Get chunks which contain a given mutation.
     * @param mut_cptr Pointer to mutation.
     * @return Vector of pointers to chunks that have this mutation.
     */
    //std::vector< Read_Chunk_CPtr > get_chunks_with_mutation(const Mutation* mut_cptr) const;

    /** Reverse the contig. */
    void reverse();

    /** Get out-edges counts.
     * @return A tuple (cnt_left, uniq_left, cnt_right, uniq_right), where cnt is the number
     * of read chunks spanning that breakpoint, and uniq is the number of different contig entries
     * where following chunks are mapped.
     */
    //std::tuple< size_t, size_t, size_t, size_t > get_out_degrees() const;

    /** Retrieve read chunks completely spanning the given interval.
     * @param start Start of the interval, 0-based, closed.
     * @param end End of the interval, 0-based, open.
     */
    //std::shared_ptr< std::vector< Read_Chunk_CPtr > > get_chunks_spanning_pos(Size_Type start, Size_Type end) const;

    /** Retrieve chunks leaving contig.
     * @param c_right Bool; true: spanning past contig end, false: spanning past contig start.
     * @param skip_next_unmappable Bool; true: ignore chunks whose next chunk is unmappable.
     * @return A map containing, for each Contig_Entry and orientation,
     * a list of Read_Chunk objects whose sibling is in that Contig_Entry and has the same orientation;
     * two chunks have the same orientation if they are both mapped to the positive or negative
     * strand of their respective contigs.
     */
    std::map< std::tuple< Contig_Entry_CBPtr, bool >, std::vector< Read_Chunk_CBPtr > >
    out_chunks_dir(bool c_right, bool skip_next_unmappable = true) const;

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

    /** Integrity check. */
    bool check() const;
    //bool check_colour(bool dir) const;

    friend std::ostream& operator << (std::ostream& os, const Contig_Entry& rhs);

private:
    Seq_Type _seq;
    Mutation_Cont _mut_cont;
    Read_Chunk_CE_Cont _chunk_cont;
    Size_Type _seq_offset;
    size_t _contig_id;
    int _colour;

    /** Hooks for storage in intrusive list in Graph object. */
    friend struct Contig_Entry_List_Node_Traits;
    Contig_Entry_BPtr _previous;
    Contig_Entry_BPtr _next;
    bool is_unlinked() const { return not(_previous or _next); }
}; // class Contig_Entry

} // namespace MAC


#endif
