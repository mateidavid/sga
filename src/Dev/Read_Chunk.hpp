//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __READ_CHUNK_HPP
#define __READ_CHUNK_HPP

#include <iostream>
#include <tuple>
#include <vector>
#include <set>
#include <map>
#include <memory>
#include <functional>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>

#include "MAC_forward.hpp"
#include "Read_Chunk_forward.hpp"
#include "Mutation.hpp"
#include "Cigar.hpp"


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
    Mutation_Ptr_Cont::const_iterator mut_ptr_node_cit;
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
                   Mutation_Ptr_Cont::const_iterator _mut_ptr_node_cit = Mutation_Ptr_Cont::const_iterator(),
                   Size_Type _mut_offset = 0,
                   const Read_Chunk* _rc_cptr = NULL)
    : c_pos(_c_pos),
    r_pos(_r_pos),
    mut_ptr_node_cit(_mut_ptr_node_cit),
    mut_offset(_mut_offset),
    rc_cptr(_rc_cptr) {}

    // allow copy and move
    DEFAULT_COPY_CTOR(Read_Chunk_Pos)
    DEFAULT_MOVE_CTOR(Read_Chunk_Pos)
public:
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
        return *(mut_ptr_node_cit->get());
    }

    /** Get previous mutation. */
    const Mutation& prev_mut() const
    {
        ASSERT(rc_cptr);
        ASSERT(past_first_mut());
        auto tmp_it = mut_ptr_node_cit;
        return *((--tmp_it)->get());
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
    { if (forward) increment(brk, on_contig); else decrement(brk, on_contig); }

    /** Get position corresponding to given breakpoint.
     * @param brk Breakpoint location.
     * @param on_contig Bool; true iff this is a contig breakpoint.
     */
    void jump_to_brk(Size_Type brk, bool on_contig);

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
    bool check() const;

    /** Comparison operator. */
    friend bool operator == (const Read_Chunk_Pos&, const Read_Chunk_Pos&);
    /** To stream. */
    friend std::ostream& operator << (std::ostream&, const Read_Chunk_Pos&);
};

struct Read_Chunk_ITree_Node_Traits;
struct Read_Chunk_Set_Node_Traits;

/** Holds information about a read chunk.
 *
 * Reads will be partitioned into chunks of non-zero size as they are mapped to contigs.
 * Each read chunk object contains information about the read it comes from (start, span),
 * the contig it is mapped to (start, span, orientation), and the contig mutations it observes.
 */
class Read_Chunk
{
private:
    // Can only be constructed by Factory object
    friend class Factory< Read_Chunk >;

    /** @name Constructors */
    /**@{*/
    /** Empty constructor. */
    Read_Chunk()
    : _re_bptr(NULL), _ce_bptr(NULL), _mut_ptr_cont(), _r_start(0), _r_len(0), _c_start(0), _c_len(0), _rc(false), _is_unmappable(false) {}

    /** Constructs a single read chunk from a read, and maps it to a new contig.
     * @param re_ptr Pointer to read entry where this chunk comes from.
     * @param len Length of the read.
     */
    Read_Chunk(Read_Entry_BPtr re_bptr, Size_Type len)
    : _re_bptr(re_bptr), _ce_bptr(NULL), _mut_ptr_cont(), _r_start(0), _r_len(len), _c_start(0), _c_len(len), _rc(false), _is_unmappable(false) {}
    /**@}*/

    // allow move only when unlinked
    DELETE_COPY_CTOR(Read_Chunk)
    Read_Chunk(Read_Chunk&& rhs) { *this = std::move(rhs); }
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
            _rc = std::move(rhs._rc);
            _is_unmappable = std::move(rhs._is_unmappable);
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
    bool get_rc() const { return _rc; }
    const Mutation_Ptr_Cont& mut_ptr_cont() const { return _mut_ptr_cont; }
    Mutation_Ptr_Cont& mut_ptr_cont() { return _mut_ptr_cont; }
    Seq_Type get_seq() const;
    Seq_Type substr(Size_Type start, Size_Type len) const;
    bool is_unmappable() const { return _is_unmappable; }
    void set_unmappable() { _is_unmappable = true; }
    Size_Type get_read_len() const;
    /**@}*/

    /** Set of coordinates used to traverse a Read_Chunk mapping. */
    /** Get mapping start position. */
    Read_Chunk_Pos get_start_pos() const
    {
        return Read_Chunk_Pos(get_c_start(), (not _rc? get_r_start() : get_r_end()), _mut_ptr_cont.begin(), 0, this);
    }

    /** Get mapping end position. */
    Read_Chunk_Pos get_end_pos() const
    {
        return Read_Chunk_Pos(get_c_end(), (not _rc? get_r_end() : get_r_start()), --_mut_ptr_cont.end(), 0, this);
    }

    /** Find bounded pointer to this object.
     * Pre: Must be linked.
     */
    Read_Chunk_CBPtr bptr_to() const
    {
        ASSERT(not is_unlinked());
        if (_re_parent->_re_l_child.raw() == this)
        {
            return _re_parent->_re_l_child;
        }
        if (_re_parent->_re_r_child.raw() == this)
        {
            return _re_parent->_re_r_child;
        }
        if (_re_parent->_re_parent.raw() == this)
        {
            return _re_parent->_re_parent;
        }
        ASSERT(false);
        return nullptr;
    }
    Read_Chunk_BPtr bptr_to()
    {
        return static_cast< Read_Chunk_BPtr >(const_cast< const Read_Chunk* >(this)->bptr_to());
    }

    /** Assign read chunk to contig.
     * @param ce_cptr Pointer to contig entry object.
     * @param c_start Start of contig base sequence this chunk is mapped to.
     * @param c_len Length of contig base sequence this chunk is mapped to.
     * @param rc Orientation of the mapping.
     * @param mut_ptr_cont Contig mutations observed by this chunk; container is taken over.
     */
    void assign_to_contig(Contig_Entry_BPtr ce_bptr, Size_Type c_start, Size_Type c_len, bool rc, Mutation_Ptr_Cont&& mut_ptr_cont)
    {
        _ce_bptr = ce_bptr;
        _c_start = c_start;
        _c_len = c_len;
        _rc = rc;
        _mut_ptr_cont = std::move(mut_ptr_cont);
    }

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
     * @param rf_ptr Reference string pointer (either whole, or just mapped portion); Contig_Entry object takes ownership.
     * @param qr Query string (either whole, or just mapped portion).
     */
    Read_Chunk_BPtr make_chunk_from_cigar(const Cigar& cigar, Seq_Type* rf_ptr, const Seq_Type& qr = std::string());

    /** Create Read_Chunk object and corresponding Mutation container from a cigar string and 2 existing Read_Chunk objects.
     * @param cigar Cigar string.
     * @param rc1 Read_Chunk corresponding to rf.
     * @param rc2 Read_Chunk corresponding to qr.
     * @return A read chunk corresponding to the query in the cigar object, and a Contig_Entry corresponding to the rf.
     */
    //static std::tuple< std::shared_ptr< Read_Chunk >, std::shared_ptr< Contig_Entry > > make_chunk_from_cigar_and_chunks(
    //    const Cigar& cigar, const Read_Chunk& rc1, const Read_Chunk& rc2);

    /** Construct a Read_Chunk corresponding to the read->contig mapping.
     * @return New Read_Chunk object and Contig_Entry objects;
     * and a Mutation translation object from original to reversed Mutations.
     */
    //std::tuple< std::shared_ptr< Read_Chunk >, std::shared_ptr< Contig_Entry >, std::shared_ptr< Mutation_Trans_Cont > >
    //reverse() const;

    /** Collapse mutations corresponding to 2 mappings.
     * @param lhs Read_Chunk object corresponding to contig->read1 mapping.
     * @param lhs_me Container of Mutation_Extra objects from lhs.
     * @param rhs Read_Chunk object corresponding to read1->read2 mapping.
     * @return vector of contig->read2 mutations; for each:
     * a bool specifying if it is a contig mutation (true) or a read1 mutation (false);
     * pointer to the mutation object;
     * and a mutation start&length.
     *
     *         static std::vector< std::tuple< bool, Read_Chunk_CPtr, Size_Type, Size_Type > > collapse_mutations(
     *             const Read_Chunk& rc1, const Mutation_Extra_Cont& rc1_me_cont, const Read_Chunk& rc2);
     */
    //std::shared_ptr< Read_Chunk > collapse_mapping(const Read_Chunk& rc2, Mutation_Cont& extra_mut_cont) const;

    /** Reverse the contig mapping (assumes mutations are reverse by contig entry). */
    //void reverse();

    /** Merge this read chunk with the next chunk of the same read.
     * Pre: Chunks must be mapped to the same contig, in the same orientation, continuously.
     * @param rc_next_cptr CPtr to next chunk.
     * @param add_mut_mod Modifier that allows chunk object to add mutations to contig.
     */
    //void merge_next(const Read_Chunk* rc_next_cptr, Mutation::add_mut_mod_type add_mut_mod);

    /** Rebase this chunk into another contig.
     * @param ce_cptr New contig.
     * @param mut_map Mutation translation map.
     * @param prefix_len Length of prefix by which new contig is larger than the old one.
     */
    //void rebase(const Contig_Entry* ce_cptr, const Mutation_Trans_Cont& mut_map, Size_Type prefix_len);

    bool check() const;

    friend std::ostream& operator << (std::ostream&, const Read_Chunk&);

private:
    Read_Entry_BPtr _re_bptr;
    Contig_Entry_BPtr _ce_bptr;
    Mutation_Ptr_Cont _mut_ptr_cont;
    Size_Type _r_start;
    Size_Type _r_len;
    Size_Type _c_start;
    Size_Type _c_len;
    bool _rc;
    bool _is_unmappable;

    /** Hooks for storage:
     * 1. in intrusive interval trees inside Contig_Entry objects.
     * 2. in intrusive trees inside Read_Entry objects.
     */
    friend struct Read_Chunk_ITree_Node_Traits;
    friend struct Read_Chunk_Set_Node_Traits;
    Read_Chunk_BPtr _ce_parent;
    Read_Chunk_BPtr _ce_l_child;
    Read_Chunk_BPtr _ce_r_child;
    Read_Chunk_BPtr _re_parent;
    Read_Chunk_BPtr _re_l_child;
    Read_Chunk_BPtr _re_r_child;
    Size_Type _ce_max_end;
    bool _ce_col;
    bool _re_col;
    bool is_unlinked() const { return not(_ce_parent or _ce_l_child or _ce_r_child or _re_parent or _re_l_child or _re_r_child); }
}; // class Read_Chunk

struct Read_Chunk_ITree_Node_Traits
{
    typedef Holder< Read_Chunk > node;
    typedef Read_Chunk_Fact fact_type;
    typedef fact_type::ptr_type node_ptr;
    typedef fact_type::const_ptr_type const_node_ptr;
    typedef int color;
    typedef Size_Type key_type;

    static node_ptr get_parent(const_node_ptr n) { return n->_ce_parent; }
    static void set_parent(node_ptr n, node_ptr ptr) { n->_ce_parent = ptr; }
    static node_ptr get_left(const_node_ptr n) { return n->_ce_l_child; }
    static void set_left(node_ptr n, node_ptr ptr) { n->_ce_l_child = ptr; }
    static node_ptr get_right(const_node_ptr n) { return n->_ce_r_child; }
    static void set_right(node_ptr n, node_ptr ptr) { n->_ce_r_child = ptr; }
    static color get_color(const_node_ptr n) { return n->_ce_col; }
    static void set_color(node_ptr n, color c) { n->_ce_col = c ; }
    static color black() { return 0; }
    static color red() { return 1; }
    static key_type get_max_end(const_node_ptr n) { return n->_ce_max_end; }
    static void set_max_end(node_ptr n, key_type k) { n->_ce_max_end = k ; }
};

struct Read_Chunk_ITree_Value_Traits
{
    typedef Read_Chunk value_type;
    typedef Read_Chunk_ITree_Node_Traits node_traits;
    typedef node_traits::key_type key_type;
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
    static key_type get_start(const_pointer n) { return n->get_c_start(); }
    static key_type get_start(const value_type* n) { return n->get_c_start(); }
    static key_type get_end(const_pointer n) { return n->get_c_end(); }
    static key_type get_end(const value_type* n) { return n->get_c_end(); }
};

class Read_Chunk_CE_Cont
    : private boost::intrusive::itree< Read_Chunk_ITree_Value_Traits >
{
private:
    typedef boost::intrusive::itree< Read_Chunk_ITree_Value_Traits > Base;

public:
    // use base constructors
    //using Base::Base;

    // allow move only
    DEFAULT_DEF_CTOR(Read_Chunk_CE_Cont)
    DELETE_COPY_CTOR(Read_Chunk_CE_Cont)
    DEFAULT_MOVE_CTOR(Read_Chunk_CE_Cont)
    DELETE_COPY_ASOP(Read_Chunk_CE_Cont)
    DEFAULT_MOVE_ASOP(Read_Chunk_CE_Cont)

    // check it is empty when deallocating
    ~Read_Chunk_CE_Cont() { ASSERT(this->size() == 0); }

    using Base::size;
    using Base::begin;
    using Base::end;

    /** Insert read chunk in this container. */
    void insert(Read_Chunk_BPtr rc_bptr)
    {
        static_cast< Base* >(this)->insert(*rc_bptr);
    }
};

struct Read_Chunk_Set_Node_Traits
{
    typedef Holder< Read_Chunk > node;
    typedef Read_Chunk_Fact fact_type;
    typedef fact_type::ptr_type node_ptr;
    typedef fact_type::const_ptr_type const_node_ptr;
    typedef int color;

    static node_ptr get_parent(const_node_ptr n) { return n->_re_parent; }
    static void set_parent(node_ptr n, node_ptr ptr) { n->_re_parent = ptr; }
    static node_ptr get_left(const_node_ptr n) { return n->_re_l_child; }
    static void set_left(node_ptr n, node_ptr ptr) { n->_re_l_child = ptr; }
    static node_ptr get_right(const_node_ptr n) { return n->_re_r_child; }
    static void set_right(node_ptr n, node_ptr ptr) { n->_re_r_child = ptr; }
    static color get_color(const_node_ptr n) { return n->_re_col; }
    static void set_color(node_ptr n, color c) { n->_re_col = c ; }
    static color black() { return 0; }
    static color red() { return 1; }
};

struct Read_Chunk_Set_Value_Traits
{
    typedef Read_Chunk value_type;
    typedef Read_Chunk_Set_Node_Traits node_traits;
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

/** Comparator for storage in RE Cont. */
struct Read_Chunk_Set_Comparator
{
    bool operator () (const Read_Chunk& lhs, const Read_Chunk& rhs) const
    {
        return lhs.get_r_start() < rhs.get_r_start();
    }
};

class Read_Chunk_RE_Cont
    : boost::intrusive::set<
          Read_Chunk,
          boost::intrusive::compare< Read_Chunk_Set_Comparator >,
          boost::intrusive::value_traits< Read_Chunk_Set_Value_Traits >
          >
{
public:
    typedef boost::intrusive::set<
                Read_Chunk,
                boost::intrusive::compare< Read_Chunk_Set_Comparator >,
                boost::intrusive::value_traits< Read_Chunk_Set_Value_Traits >
                > Base;

    // use base constructors
    using Base::Base;

    // allow move only
    DEFAULT_DEF_CTOR(Read_Chunk_RE_Cont)
    DELETE_COPY_CTOR(Read_Chunk_RE_Cont)
    DEFAULT_MOVE_CTOR(Read_Chunk_RE_Cont)
    DELETE_COPY_ASOP(Read_Chunk_RE_Cont)
    DEFAULT_MOVE_ASOP(Read_Chunk_RE_Cont)

    // check it is empty when deallocating
    ~Read_Chunk_RE_Cont() { ASSERT(this->size() == 0); }

    /** Insert read chunk in this container. */
    void insert(Read_Chunk_BPtr rc_bptr)
    {
        Base::iterator it;
        bool success;
        std::tie(it, success) = static_cast< Base* >(this)->insert(*rc_bptr);
        ASSERT(success);
    }

    /** Find chunk which contains given read position.
     * @param r_pos Read position, 0-based.
     * @return Pointer to Read Chunk object, or NULL if no chunk contains r_pos.
     */
    Read_Chunk_CBPtr get_chunk_with_pos(Size_Type r_pos) const
    {
        ASSERT(this->size() > 0);
        ASSERT(this->begin()->re_bptr());
        if (r_pos >= this->begin()->get_read_len())
        {
            return nullptr;
        }
        Base::const_iterator cit = this->lower_bound(r_pos, Base::value_compare());
        if (cit != this->cend() and cit->get_r_start() == r_pos)
        {
            return &*cit;
        }
        else
        {
            ASSERT(cit != this->cbegin());
            return &*(--cit);
        }
    }

    Read_Chunk_BPtr get_chunk_with_pos(Size_Type r_pos)
    {
        return static_cast< Read_Chunk_BPtr >(const_cast< const Read_Chunk_RE_Cont* >(this)->get_chunk_with_pos(r_pos));
    }

    /** Get the sibling of the given read chunk.
     * @param next Bool: if true, get next chunk; if false, get previous chunk.
     * @return Pointer to sibling chunk, or NULL if no sibling exists.
     */
    Read_Chunk_CBPtr get_next(Read_Chunk_CBPtr rc_cbptr, bool next) const
    {
        Base::const_iterator rc_cit = this->iterator_to(*rc_cbptr);
        if (next)
        {
            ++rc_cit;
            if (rc_cit == this->cend())
            {
                return nullptr;
            }
        }
        else
        {
            if (rc_cit == this->cbegin())
            {
                return nullptr;
            }
            --rc_cit;
        }
        return &*rc_cit;
    }
    Read_Chunk_BPtr get_next(Read_Chunk_CBPtr rc_cbptr, bool next)
    {
        return static_cast< Read_Chunk_BPtr >(const_cast< const Read_Chunk_RE_Cont* >(this)->get_next(rc_cbptr, next));
    }

};

} // namespace MAC


#endif
