//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __CONTIG_ENTRY_CONT_HPP
#define __CONTIG_ENTRY_CONT_HPP

#include <boost/intrusive/list.hpp>

#include "MAC_forward.hpp"
#include "Contig_Entry.hpp"


namespace MAC
{

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

typedef boost::intrusive::list< Contig_Entry,
                                boost::intrusive::value_traits< Contig_Entry_List_Value_Traits >
                              > Contig_Entry_Cont;

} // namespace MAC


#endif
