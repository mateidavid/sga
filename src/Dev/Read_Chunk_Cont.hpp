//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __READ_CHUNK_CONT_HPP
#define __READ_CHUNK_CONT_HPP

#include "Read_Chunk.hpp"


namespace MAC
{

namespace detail
{

struct Read_Chunk_ITree_Node_Traits
{
    typedef Read_Chunk node;
    typedef Read_Chunk_Fact fact_type;
    typedef fact_type::ptr_type node_ptr;
    typedef fact_type::const_ptr_type const_node_ptr;
    typedef bool color;
    typedef Size_Type key_type;

    static node_ptr get_parent(const_node_ptr n) { return n->_ce_parent; }
    static void set_parent(node_ptr n, node_ptr ptr) { n->_ce_parent = ptr; }
    static node_ptr get_left(const_node_ptr n) { return n->_ce_l_child; }
    static void set_left(node_ptr n, node_ptr ptr) { n->_ce_l_child = ptr; }
    static node_ptr get_right(const_node_ptr n) { return n->_ce_r_child; }
    static void set_right(node_ptr n, node_ptr ptr) { n->_ce_r_child = ptr; }
    static color get_color(const_node_ptr n) { return n->_get_ce_col(); }
    static void set_color(node_ptr n, color c) { n->_set_ce_col(c); }
    static color black() { return false; }
    static color red() { return true; }
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
    typedef Read_Chunk_ITree_Value_Traits* value_traits_ptr;

    static const bi::link_mode_type link_mode = bi::safe_link;

    static node_ptr to_node_ptr (reference value) { return &value; }
    static const_node_ptr to_node_ptr (const_reference value) { return &value; }
    static pointer to_value_ptr(node_ptr n) { return n; }
    static const_pointer to_value_ptr(const_node_ptr n) { return n; }
    static key_type get_start(const_pointer n) { return n->get_c_start(); }
    static key_type get_start(const value_type* n) { return n->get_c_start(); }
    static key_type get_end(const_pointer n) { return n->get_c_end(); }
    static key_type get_end(const value_type* n) { return n->get_c_end(); }
};

} // namespace detail

class Read_Chunk_CE_Cont
    : private bi::itree< Read_Chunk,
                         bi::value_traits< detail::Read_Chunk_ITree_Value_Traits >,
                         bi::constant_time_size< false >,
                         bi::header_holder_type< bounded::Pointer_Holder< Read_Chunk > >
                       >
{
private:
    typedef bi::itree< Read_Chunk,
                       bi::value_traits< detail::Read_Chunk_ITree_Value_Traits >,
                       bi::constant_time_size< false >,
                       bi::header_holder_type< bounded::Pointer_Holder< Read_Chunk > >
                     > Base;

public:
    // allow move only
    DEFAULT_DEF_CTOR(Read_Chunk_CE_Cont)
    DELETE_COPY_CTOR(Read_Chunk_CE_Cont)
    DEFAULT_MOVE_CTOR(Read_Chunk_CE_Cont)
    DELETE_COPY_ASOP(Read_Chunk_CE_Cont)
    DEFAULT_MOVE_ASOP(Read_Chunk_CE_Cont)

    USING_INTRUSIVE_CONT(Base)
    using typename Base::intersection_const_iterator;
    using typename Base::intersection_const_iterator_range;
    using Base::iintersect;
    using Base::max_end;
    using Base::clear_and_dispose;
    using Base::check;

    // check it is empty when deallocating
    ~Read_Chunk_CE_Cont() { ASSERT(empty()); }

    Base::size_type size() = delete;
    bool single_node() const
    {
        return not empty()
            and not Base::node_traits::get_left(Base::node_traits::get_parent(Base::header_ptr()))
            and not Base::node_traits::get_right(Base::node_traits::get_parent(Base::header_ptr()));
    }

    /** Insert Read_Chunk in this container. */
    void insert(Read_Chunk_BPtr rc_bptr) { Base::insert(*rc_bptr); }

    /** Erase Read_Chunk from container. */
    void erase(Read_Chunk_CBPtr rc_cbptr) { Base::erase(iterator_to(*rc_cbptr)); }

    /** Split container and its Read_Chunk objects at a contig breakpoint.
     * Pre: No Mutations may span c_pos.
     * Pre: Chunks must be unlinked from their RE container.
     * NOTE: Chunks that do not span the break stay unchanged;
     * For Chunks that get split, we also fix Chunk back pointers
     * in their Mutation ptr containers.
     * @param c_brk Contig position of the breakpoint.
     * @param mut_left_cbptr Insertion at c_pos to remain on the left of the cut, if any.
     * @return New Read_Chunk_CE_Cont containing rhs.
     */
    Read_Chunk_CE_Cont split(Size_Type c_brk, Mutation_CBPtr mut_left_cbptr);

    /** Erase all Chunks from their respective RE container. */
    void erase_from_re_cont() const;
    /** Insert all Chunks into their respective RE container. */
    void insert_into_re_cont() const;

    /** Shift contig coordinates of all Read_Chunk objects in this container.
     * @param delta Signed integer value to add to all contig start points.
     */
    void shift(int delta)
    {
        for (auto rc_bref : *this)
        {
            Read_Chunk_BPtr rc_bptr = &rc_bref;
            rc_bptr->shift(delta);
        }
        Base::implement_shift(delta);
    }

    /** Set ce_ptr of all chunks.
     * @param ce_ptr New ce pointer.
     */
    void set_ce_ptr(Contig_Entry_BPtr ce_bptr)
    {
        for (auto rc_bref : *this)
        {
            Read_Chunk_BPtr rc_bptr = &rc_bref;
            rc_bptr->ce_bptr() = ce_bptr;
        }
    }

    /** Copy chunks from given container into this one.
     * @param other_cont Container to clear.
     * @param ce_bptr Pointer to Contig_Entry object of this container.
     */
    void splice(Read_Chunk_CE_Cont& other_cont, Contig_Entry_BPtr ce_bptr)
    {
        static_cast< Base& >(other_cont).clear_and_dispose([&] (Read_Chunk_BPtr rc_bptr)
        {
            rc_bptr->ce_bptr() = ce_bptr;
            insert(rc_bptr);
        });
    }

    /** Clear container and deallocate all Read_Chunk objects.
     * Pre: Read_Chunk objects must be unlinked from their RE containers.
     * Post: MCAs are also deallocated, but not the Mutations themselves.
     */
    void clear_and_dispose()
    {
        Base::clear_and_dispose([] (Read_Chunk_BPtr rc_bptr)
        {
            rc_bptr->mut_ptr_cont().clear_and_dispose();
            Read_Chunk_Fact::del_elem(rc_bptr);
        });
    }

    /** Reverse all chunks in this container.
     */
    void reverse()
    {
        Read_Chunk_CE_Cont new_cont;
        Base::clear_and_dispose([&new_cont] (Read_Chunk_BPtr rc_bptr)
        {
            rc_bptr->mut_ptr_cont().reverse();
            rc_bptr->reverse();
            new_cont.insert(rc_bptr);
        });
        *this = std::move(new_cont);
    }

}; // class Read_Chunk_CE_Cont

namespace detail
{

struct Read_Chunk_Set_Node_Traits
{
    typedef Read_Chunk node;
    typedef Read_Chunk_Fact fact_type;
    typedef fact_type::ptr_type node_ptr;
    typedef fact_type::const_ptr_type const_node_ptr;
    typedef bool color;

    static node_ptr get_parent(const_node_ptr n) { return n->_re_parent; }
    static void set_parent(node_ptr n, node_ptr ptr) { n->_re_parent = ptr; }
    static node_ptr get_left(const_node_ptr n) { return n->_re_l_child; }
    static void set_left(node_ptr n, node_ptr ptr) { n->_re_l_child = ptr; }
    static node_ptr get_right(const_node_ptr n) { return n->_re_r_child; }
    static void set_right(node_ptr n, node_ptr ptr) { n->_re_r_child = ptr; }
    static color get_color(const_node_ptr n) { return n->_get_re_col(); }
    static void set_color(node_ptr n, color c) { n->_set_re_col(c); }
    static color black() { return false; }
    static color red() { return true; }
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
    typedef Read_Chunk_Set_Value_Traits* value_traits_ptr;

    static const bi::link_mode_type link_mode = bi::safe_link;

    static node_ptr to_node_ptr (reference value) { return &value; }
    static const_node_ptr to_node_ptr (const_reference value) { return &value; }
    static pointer to_value_ptr(node_ptr n) { return n; }
    static const_pointer to_value_ptr(const_node_ptr n) { return n; }
};

} // namespace detail

/** Comparator for storage in RE Cont. */
struct Read_Chunk_Set_Comparator
{
    bool operator () (const Read_Chunk& lhs, const Read_Chunk& rhs) const
    {
        return lhs.get_r_start() < rhs.get_r_start();
    }
    bool operator () (const Read_Chunk& lhs, size_t rhs_val) const
    {
        return lhs.get_r_start() < rhs_val;
    }
    bool operator () (size_t lhs_val, const Read_Chunk& rhs) const
    {
        return lhs_val < rhs.get_r_start();
    }
};

class Read_Chunk_RE_Cont
    : private bi::set< Read_Chunk,
                       bi::compare< Read_Chunk_Set_Comparator >,
                       bi::value_traits< detail::Read_Chunk_Set_Value_Traits >,
                       bi::header_holder_type< bounded::Pointer_Holder< Read_Chunk > >
                     >
{
private:
    typedef bi::set< Read_Chunk,
                     bi::compare< Read_Chunk_Set_Comparator >,
                     bi::value_traits< detail::Read_Chunk_Set_Value_Traits >,
                     bi::header_holder_type< bounded::Pointer_Holder< Read_Chunk > >
                   > Base;
public:
    // allow move only
    DEFAULT_DEF_CTOR(Read_Chunk_RE_Cont)
    DELETE_COPY_CTOR(Read_Chunk_RE_Cont)
    DEFAULT_MOVE_CTOR(Read_Chunk_RE_Cont)
    DELETE_COPY_ASOP(Read_Chunk_RE_Cont)
    DEFAULT_MOVE_ASOP(Read_Chunk_RE_Cont)

    USING_INTRUSIVE_CONT(Base)
    using Base::check;

    // check it is empty when deallocating
    ~Read_Chunk_RE_Cont() { ASSERT(empty()); }

    Base::size_type size() = delete;

    /** Insert read chunk in this container. */
    void insert(Read_Chunk_BPtr rc_bptr)
    {
        iterator it;
        bool success;
        std::tie(it, success) = static_cast< Base* >(this)->insert(*rc_bptr);
        ASSERT(success);
    }

    /** Erase Read_Chunk from container. */
    void erase(Read_Chunk_CBPtr rc_cbptr) { Base::erase(iterator_to(*rc_cbptr)); }

    /** Find chunk which contains given read position.
     * @param r_pos Read position, 0-based.
     * @return Pointer to Read Chunk object, or NULL if no chunk contains r_pos.
     */
    Read_Chunk_CBPtr get_chunk_with_pos(Size_Type r_pos) const
    {
        ASSERT(not empty());
        ASSERT(begin()->re_bptr());
        if (r_pos >= begin()->get_read_len())
        {
            return nullptr;
        }
        auto cit = Base::lower_bound(r_pos, Base::value_compare());
        if (cit == end() or cit->get_r_start() != r_pos)
        {
            ASSERT(cit != begin());
            --cit;
        }
        return &*cit;
    }

    /** Get the sibling of the given read chunk.
     * @param rc_cbptr Original Read_Chunk.
     * @param read Bool; true: right/left wrt to read; false: right/left wrt to contig.
     * @param right Bool; true: get chunk to the right; false: get chunk to the left.
     * @return Pointer to sibling chunk, or NULL if no sibling exists.
     */
    Read_Chunk_CBPtr get_sibling(Read_Chunk_CBPtr rc_cbptr, bool read, bool right) const
    {
        const_iterator rc_cit = iterator_to(*rc_cbptr);
        bool r_right = (read? right : right != rc_cbptr->get_rc());
        if (r_right)
        {
            ++rc_cit;
            if (rc_cit == this->end())
            {
                return nullptr;
            }
        }
        else
        {
            if (rc_cit == this->begin())
            {
                return nullptr;
            }
            --rc_cit;
        }
        return &*rc_cit;
    }
}; // class Read_Chunk_RE_Cont

} // namespace MAC

#endif
