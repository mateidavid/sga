//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __READ_ENTRY_CONT_HPP
#define __READ_ENTRY_CONT_HPP

#include <boost/intrusive/set.hpp>

#include "Read_Entry.hpp"


namespace MAC
{

namespace detail
{

struct Read_Entry_Set_Node_Traits
{
    typedef Read_Entry node;
    typedef Read_Entry_Fact fact_type;
    typedef fact_type::ptr_type node_ptr;
    typedef fact_type::const_ptr_type const_node_ptr;
    typedef int color;

    static node_ptr get_parent(const_node_ptr n) { return n->_parent; }
    static void set_parent(node_ptr n, node_ptr ptr) { n->_parent = ptr; }
    static node_ptr get_left(const_node_ptr n) { return n->_l_child; }
    static void set_left(node_ptr n, node_ptr ptr) { n->_l_child = ptr; }
    static node_ptr get_right(const_node_ptr n) { return n->_r_child; }
    static void set_right(node_ptr n, node_ptr ptr) { n->_r_child = ptr; }
    static color get_color(const_node_ptr n) { return n->_col; }
    static void set_color(node_ptr n, color c) { n->_col = c ; }
    static color black() { return 0; }
    static color red() { return 1; }
};

struct Read_Entry_Set_Value_Traits
{
    typedef Read_Entry value_type;
    typedef Read_Entry_Set_Node_Traits node_traits;
    typedef node_traits::node_ptr node_ptr;
    typedef node_traits::const_node_ptr const_node_ptr;
    typedef node_ptr pointer;
    typedef const_node_ptr const_pointer;
    typedef node_traits::fact_type::ref_type reference;
    typedef node_traits::fact_type::const_ref_type const_reference;
    typedef Read_Entry_Set_Value_Traits* value_traits_ptr;

    static const bi::link_mode_type link_mode = bi::safe_link;

    static node_ptr to_node_ptr (reference value) { return &value; }
    static const_node_ptr to_node_ptr (const_reference value) { return &value; }
    static pointer to_value_ptr(node_ptr n) { return n; }
    static const_pointer to_value_ptr(const_node_ptr n) { return n; }
};

} // namespace detail

/** Comparator for storage in tree. */
struct Read_Entry_Comparator
{
    bool operator () (const Read_Entry& lhs, const Read_Entry& rhs) const
    {
        return lhs.name() < rhs.name();
    }
    bool operator () (const Read_Entry& lhs, const std::string& rhs_name) const
    {
        return lhs.name() < rhs_name;
    }
    bool operator () (const std::string& lhs_name, const Read_Entry& rhs) const
    {
        return lhs_name < rhs.name();
    }
};

class Read_Entry_Cont
    : private bi::set< Read_Entry,
                       bi::compare< Read_Entry_Comparator >,
                       bi::value_traits< detail::Read_Entry_Set_Value_Traits >,
                       bi::header_holder_type< bounded::Pointer_Holder< Read_Entry > >
                     >
{
private:
    typedef bi::set< Read_Entry,
                     bi::compare< Read_Entry_Comparator >,
                     bi::value_traits< detail::Read_Entry_Set_Value_Traits >,
                     bi::header_holder_type< bounded::Pointer_Holder< Read_Entry > >
                   > Base;
public:
    DEFAULT_DEF_CTOR(Read_Entry_Cont)
    DELETE_COPY_CTOR(Read_Entry_Cont)
    DEFAULT_MOVE_CTOR(Read_Entry_Cont)
    DELETE_COPY_ASOP(Read_Entry_Cont)
    DEFAULT_MOVE_ASOP(Read_Entry_Cont)

    USING_INTRUSIVE_CONT(Base)
    friend class Graph;

    // check it is empty when deallocating
    ~Read_Entry_Cont() { ASSERT(empty()); }

    /** Insert Read_Entry into container.
     * @param re_bptr Pointer to Read_Entry object.
     */
    void insert(Read_Entry_BPtr re_bptr)
    {
        auto res = Base::insert(*re_bptr);
        ASSERT(res.second);
    }

    /** Search for Read_Entry.
     * @param name Name of read to look for.
     * @return Read_Entry pointer, or nullptr if not found.
     */
    Read_Entry_CBPtr find(const std::string& name) const
    {
        auto cit = Base::find(name, Base::value_compare());
        return cit != end()? &*cit : nullptr;
    }

    /** Clear container and deallocate RE objects.
     * Pre: Chunks must have been removed from their RE containers.
     */
    void clear_and_dispose()
    {
        Base::clear_and_dispose([] (Read_Entry_BPtr re_bptr)
        {
            Read_Entry_Fact::del_elem(re_bptr);
        });
    }
};

} // namespace MAC

#endif
