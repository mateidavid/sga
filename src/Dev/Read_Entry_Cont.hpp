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

struct Read_Entry_Set_Node_Traits
{
    typedef Holder< Read_Entry > node;
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

    static const boost::intrusive::link_mode_type link_mode = boost::intrusive::normal_link;

    static node_ptr to_node_ptr (reference value) { return &value; }
    static const_node_ptr to_node_ptr (const_reference value) { return &value; }
    static pointer to_value_ptr(node_ptr n) { return n; }
    static const_pointer to_value_ptr(const_node_ptr n) { return n; }
};

typedef boost::intrusive::set< Read_Entry,
                               boost::intrusive::compare< Read_Entry_Comparator >,
                               boost::intrusive::value_traits< Read_Entry_Set_Value_Traits >
                             > Read_Entry_Cont;

} // namespace MAC

#endif
