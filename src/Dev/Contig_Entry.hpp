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
     * @param seq_ptr Pointer to read sequence (assume ownership).
     */
    Contig_Entry(Seq_Type* seq_ptr, Size_Type seq_offset = 0) : _seq_ptr(seq_ptr), _seq_offset(seq_offset), _colour(0) {}

    /** Constructor out of a Read_Chunk. * /
    Contig_Entry(Read_Chunk_CPtr rc_cptr) : _seq_ptr(new Seq_Type(rc_cptr->get_seq())), _seq_offset(0), _colour(0)
    {
        _chunk_cptr_cont.push_back(rc_cptr);
    }
    */

    // allow move only when unlinked
    DELETE_COPY_CTOR(Contig_Entry)
    Contig_Entry(Contig_Entry&& rhs) { *this = std::move(rhs); }
public:
    DELETE_COPY_ASOP(Contig_Entry)
    Contig_Entry& operator = (Contig_Entry&& rhs)
    {
        if (this != &rhs)
        {
            ASSERT(is_unlinked());
            _seq_ptr = std::move(rhs._seq_ptr);
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
    const Seq_Type& get_seq() const { return *_seq_ptr; }
    Size_Type get_seq_offset() const { return _seq_offset; }
    Seq_Type substr(Size_Type start, Size_Type len) const
    { ASSERT(start >= _seq_offset and start + len <= _seq_offset + _seq_ptr->size()); return _seq_ptr->substr(start - _seq_offset, len); }
    Size_Type get_len() const { return _seq_ptr->size(); }
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
        return static_cast< Contig_Entry_BPtr >(const_cast< const Contig_Entry* >(this)->bptr_to());
    }

    /** Add to the list of read chunks.
     * @param rc_bptr Pointer to read chunk to add.
     */
    void add_chunk(Read_Chunk_BPtr rc_bptr) { _chunk_cont.insert(rc_bptr); }

    /** Remove read chunk.
     * @param rc_cptr Pointer to read chunk to remove.
     */
    //void remove_chunk(Read_Chunk_CPtr rc_cptr);

    /** Remove read chunks.
     * @param rc_cptr_set Set of read chunks to remove.
     */
    //void remove_chunks(const std::set< Read_Chunk_CPtr >& rc_cptr_set);

    /** Cut mutation at given offsets.
     * @param mut_cptr Pointer to mutation to cut.
     * @param c_offset Contig offset inside mutation where to cut.
     * @param r_offset Read offset inside mutation where to cut.
     * @return Pointer to mutation object containing leftover part of the original mutation.
     */
    //const Mutation* cut_mutation(const Mutation* mut_cptr, Size_Type c_offset, Size_Type r_offset);

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

    /** Get mutations that completely span a given contig position.
     * @param c_pos Contig position, 0-based.
     * @return Vector of pointers to Mutation objects completely spanning c_pos.
     */
    //std::vector< const Mutation* > get_mutations_spanning_pos(Size_Type c_pos) const;

    /** Get mutations from the second half of a given Contig Entry object that is being cut in 2.
     * @param ce_cptr Contig_Entry being cut.
     * @param c_brk Position of the cut.
     * @param mut_left_cptr Pointer to insertion at c_pos that must appear on the left of the cut, if any.
     * @return Map with (key=old mutation cptr; value=new mutation cptr).
     */
    //std::map< const Mutation*, const Mutation* > acquire_second_half_mutations(
     //   const Contig_Entry* ce_cptr, Size_Type c_brk, const Mutation* mut_left_cptr);

    /** Drop mutations that appear in the map.
     * @param mut_cptr_map Map produced by acquire_second_half_mutations().
     */
    //void drop_mutations(const std::map< const Mutation*, const Mutation* >& mut_cptr_map);

    /** Drop unused mutations. */
    //void drop_unused_mutations();

    /** Drop base sequence suffix.
     * @param c_brk Prefix length to keep.
     */
    //void drop_base_seq(Size_Type c_brk) { _seq_ptr->resize(c_brk); }

    /** Reverse the contig. */
    //void reverse(const Read_Chunk::ext_mod_type& rc_reverse_mod);

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

    /** Given a chunk mapped to this contig, retreive the next chunk in the same read
     * in the given contig direction.
     * @param dir Bool; true: past contig end; false: past contig start.
     * @param rc_cptr Chunk whose sibling to look for.
     * @return Next chunk in that direction, or NULL if none exists.
     */
    //Read_Chunk_CPtr get_next_chunk(bool dir, Read_Chunk_CPtr rc_cptr) const;

    /** Retrieve chunks leaving contig.
     * @param dir Bool; true: spanning right of c_end, false: spanning left of 0.
     * @param skip_next_unmappable Bool; true: ignore chunks whose next chunk is unmappable.
     */
    //std::shared_ptr< std::vector< Read_Chunk_CPtr > > get_chunks_out(bool dir, bool skip_next_unmappable = true) const;

    /** Check if contig has a unique neighbour in the given direction.
     * @param dir Bool; true: past contig end, false: past contig start.
     * @return NULL if not; otherwise, a list of read chunks spanning past contig end
     * into another unique contig.
     */
    //std::shared_ptr< std::vector< Read_Chunk_CPtr > > is_mergeable_one_way(bool forward) const;

    /** Check if contig is mergeable with a single other contig in the given direction.
     * @param dir Bool; true: merge past contig end, false: merge past contig start.
     * @return NULL if not mergeable; otherwise, a list of read chunks spanning past contig end
     * into another unique contig.
     */
    //std::shared_ptr< std::vector< Read_Chunk_CPtr > > is_mergeable(bool dir) const;

    /** Merge the given contig into this one.
     * @param ce_next_cptr Contig to be merged into this one.
     * @param rc_rebase_mod External modifier that allows rebasing chunks into this contig.
     */
    //void merge_forward(const Contig_Entry* ce_next_cptr, const Read_Chunk::ext_mod_with_map_type& rc_rebase_mod);

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
    std::shared_ptr< Seq_Type > _seq_ptr;
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
};

struct Contig_Entry_List_Node_Traits
{
    typedef Holder< Contig_Entry > node;
    typedef Contig_Entry_Fact fact_type;
    typedef fact_type::ptr_type node_ptr;
    typedef fact_type::const_ptr_type const_node_ptr;

    static node_ptr get_previous(const_node_ptr n) { return n->_previous; }
    static void set_previous(node_ptr n, node_ptr ptr) { n->_previous = ptr; }
    static node_ptr get_next(const_node_ptr n) { return n->_next; }
    static void set_next(node_ptr n, node_ptr ptr) { n->_next = ptr; }
};

struct Contig_Entry_List_Value_Traits
{
    typedef Contig_Entry value_type;
    typedef Contig_Entry_List_Node_Traits node_traits;
    typedef node_traits::node_ptr node_ptr;
    typedef node_traits::const_node_ptr const_node_ptr;
    typedef node_ptr pointer;
    typedef const_node_ptr const_pointer;
    typedef node_traits::fact_type::ref_type reference;
    typedef node_traits::fact_type::const_ref_type const_reference;

    static const boost::intrusive::link_mode_type link_mode = boost::intrusive::normal_link;

    static node_ptr to_node_ptr (reference value) { return &value; }
    static const_node_ptr to_node_ptr (const_reference value) { return &value; }
    static pointer to_value_ptr(node_ptr n) { return n; }
    static const_pointer to_value_ptr(const_node_ptr n) { return n; }
};

typedef boost::intrusive::list<
    Contig_Entry,
    boost::intrusive::value_traits< Contig_Entry_List_Value_Traits >
> Contig_Entry_Cont;

} // namespace MAC


#endif
