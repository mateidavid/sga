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
    // allow move only
    DEFAULT_DEF_CTOR(Read_Chunk_CE_Cont)
    DELETE_COPY_CTOR(Read_Chunk_CE_Cont)
    DEFAULT_MOVE_CTOR(Read_Chunk_CE_Cont)
    DELETE_COPY_ASOP(Read_Chunk_CE_Cont)
    DEFAULT_MOVE_ASOP(Read_Chunk_CE_Cont)

    // check it is empty when deallocating
    ~Read_Chunk_CE_Cont() { ASSERT(size() == 0); }

    USING_ITERATORS(Base)

    /** Insert Read_Chunk in this container. */
    void insert(Read_Chunk_BPtr rc_bptr) { Base::insert(*rc_bptr); }

    /** Erase Read_Chunk from container. */
    void erase(Read_Chunk_CBPtr rc_cbptr) { Base::erase(*rc_cbptr); }

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
}; // class Read_Chunk_CE_Cont

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
    : private boost::intrusive::set< Read_Chunk,
                                     boost::intrusive::compare< Read_Chunk_Set_Comparator >,
                                     boost::intrusive::value_traits< Read_Chunk_Set_Value_Traits >
                                   >
{
private:
    typedef boost::intrusive::set< Read_Chunk,
                                   boost::intrusive::compare< Read_Chunk_Set_Comparator >,
                                   boost::intrusive::value_traits< Read_Chunk_Set_Value_Traits >
                                 > Base;
public:
    // allow move only
    DEFAULT_DEF_CTOR(Read_Chunk_RE_Cont)
    DELETE_COPY_CTOR(Read_Chunk_RE_Cont)
    DEFAULT_MOVE_CTOR(Read_Chunk_RE_Cont)
    DELETE_COPY_ASOP(Read_Chunk_RE_Cont)
    DEFAULT_MOVE_ASOP(Read_Chunk_RE_Cont)

    // check it is empty when deallocating
    ~Read_Chunk_RE_Cont() { ASSERT(size() == 0); }

    USING_ITERATORS(Base)

    /** Insert read chunk in this container. */
    void insert(Read_Chunk_BPtr rc_bptr)
    {
        iterator it;
        bool success;
        std::tie(it, success) = static_cast< Base* >(this)->insert(*rc_bptr);
        ASSERT(success);
    }

    /** Erase Read_Chunk from container. */
    void erase(Read_Chunk_CBPtr rc_cbptr) { Base::erase(*rc_cbptr); }

    /** Find chunk which contains given read position.
     * @param r_pos Read position, 0-based.
     * @return Pointer to Read Chunk object, or NULL if no chunk contains r_pos.
     */
    Read_Chunk_CBPtr get_chunk_with_pos(Size_Type r_pos) const
    {
        ASSERT(size() > 0);
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
     * @param next Bool: if true, get next chunk; if false, get previous chunk.
     * @return Pointer to sibling chunk, or NULL if no sibling exists.
     */
    Read_Chunk_CBPtr get_next(Read_Chunk_CBPtr rc_cbptr, bool next) const
    {
        const_iterator rc_cit = this->iterator_to(*rc_cbptr);
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
}; // class Read_Chunk_RE_Cont

} // namespace MAC

#endif
