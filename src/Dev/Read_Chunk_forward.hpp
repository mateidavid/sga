//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

// This file declares Read_Chunk-related containers used by other classes
// prior to declaring Read_Chunk itself.

#ifndef __READ_CHUNK_FORWARD_HPP
#define __READ_CHUNK_FORWARD_HPP

#include <boost/intrusive/list.hpp>

#include "MAC_forward.hpp"

namespace MAC
{

/*
//
// Read_Chunks are stored in intrusive interval trees inside Contig_Entry objects
//
struct Read_Chunk_Base
{
    static const Read_Chunk_Base* to_base(const Read_Chunk*);
    static Read_Chunk_Base* to_base(Read_Chunk* rc_ptr)
    { return const_cast< Read_Chunk_Base* >(to_base(const_cast< const Read_Chunk* >(rc_ptr))); }

    static const Read_Chunk* from_base(const Read_Chunk_Base*);
    static Read_Chunk* to_base(Read_Chunk_Base* rcb_ptr)
    { return const_cast< Read_Chunk* >(from_base(const_cast< const Read_Chunk_Base* >(rcb_ptr))); }

    Read_Chunk_BPtr _ce_parent;
    Read_Chunk_BPtr _ce_l_child;
    Read_Chunk_BPtr _ce_r_child;
    Read_Chunk_BPtr _re_parent;
    Read_Chunk_BPtr _re_l_child;
    Read_Chunk_BPtr _re_r_child;
    Size_Type _ce_max_end;
    bool _ce_col;
    bool _re_col;
    Size_Type get_c_start() const;
    Size_Type get_c_end() const;
};

struct Read_Chunk_ITree_Node_Traits
{
    typedef Holder< Read_Chunk > node;
    typedef Read_Chunk_Fact fact_type;
    typedef fact_type::ptr_type node_ptr;
    typedef fact_type::const_ptr_type const_node_ptr;
    typedef int color;
    typedef Size_Type key_type;

    static const Read_Chunk_Base* to_base(const_node_ptr n) { return Read_Chunk_Base::to_base(n.raw()); }
    static Read_Chunk_Base* to_base(node_ptr n) { return Read_Chunk_Base::to_base(n.raw()); }

    static node_ptr get_parent(const_node_ptr n) { return to_base(n)->_ce_parent; }
    static void set_parent(node_ptr n, node_ptr ptr) { to_base(n)->_ce_parent = ptr; }
    static node_ptr get_left(const_node_ptr n) { return to_base(n)->_ce_l_child; }
    static void set_left(node_ptr n, node_ptr ptr) { to_base(n)->_ce_l_child = ptr; }
    static node_ptr get_right(const_node_ptr n) { return to_base(n)->_ce_r_child; }
    static void set_right(node_ptr n, node_ptr ptr) { to_base(n)->_ce_r_child = ptr; }
    static color get_color(const_node_ptr n) { return to_base(n)->_ce_col; }
    static void set_color(node_ptr n, color c) { to_base(n)->_ce_col = c ; }
    static color black() { return 0; }
    static color red() { return 1; }
    static key_type get_max_end(const_node_ptr n) { return to_base(n)->_ce_max_end; }
    static void set_max_end(node_ptr n, key_type k) { to_base(n)->_ce_max_end = k ; }
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
    static key_type get_start(const_pointer n) { return node_traits::to_base(n)->get_c_start(); }
    static key_type get_start(const value_type* n) { return Read_Chunk_Base::to_base(n)->get_c_start(); }
    static key_type get_end(const_pointer n) { return node_traits::to_base(n)->get_c_end(); }
    static key_type get_end(const value_type* n) { return Read_Chunk_Base::to_base(n)->get_c_end(); }
};

//typedef boost::intrusive::itree< Read_Chunk_ITree_Value_Traits > Read_Chunk_CE_Cont;

//
// Read_Chunks are stored in intrusive sets inside Read_Entry objects
//
struct Read_Chunk_Set_Node_Traits
{
    typedef Holder< Read_Chunk > node;
    typedef Read_Chunk_Fact fact_type;
    typedef fact_type::ptr_type node_ptr;
    typedef fact_type::const_ptr_type const_node_ptr;
    typedef int color;

    static const Read_Chunk_Base* to_base(const_node_ptr n) { return Read_Chunk_Base::to_base(n.raw()); }
    static Read_Chunk_Base* to_base(node_ptr n) { return Read_Chunk_Base::to_base(n.raw()); }

    static node_ptr get_parent(const_node_ptr n) { return to_base(n)->_re_parent; }
    static void set_parent(node_ptr n, node_ptr ptr) { to_base(n)->_re_parent = ptr; }
    static node_ptr get_left(const_node_ptr n) { return to_base(n)->_re_l_child; }
    static void set_left(node_ptr n, node_ptr ptr) { to_base(n)->_re_l_child = ptr; }
    static node_ptr get_right(const_node_ptr n) { return to_base(n)->_re_r_child; }
    static void set_right(node_ptr n, node_ptr ptr) { to_base(n)->_re_r_child = ptr; }
    static color get_color(const_node_ptr n) { return to_base(n)->_re_col; }
    static void set_color(node_ptr n, color c) { to_base(n)->_re_col = c ; }
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

typedef boost::intrusive::set<
    Read_Chunk,
    boost::intrusive::value_traits< Read_Chunk_Set_Value_Traits >
> Read_Chunk_RE_Cont;
*/

//
// Read_Chunk pointers are stored in intrusive lists inside Mutation objects
//
struct Read_Chunk_Ptr_Node_List_Node_Traits;

class Read_Chunk_Ptr_Node
{
public:
    Read_Chunk_Ptr_Node(Read_Chunk_BPtr rc_bptr = 0) : _rc_bptr(rc_bptr) {}

    Read_Chunk_BPtr get() const { return _rc_bptr; }

private:
    const Read_Chunk_BPtr _rc_bptr;

    friend struct Read_Chunk_Ptr_Node_List_Node_Traits;
    Read_Chunk_Ptr_Node_BPtr _previous;
    Read_Chunk_Ptr_Node_BPtr _next;
};

struct Read_Chunk_Ptr_Node_List_Node_Traits
{
    typedef Holder< Read_Chunk_Ptr_Node > node;
    typedef Read_Chunk_Ptr_Node_Fact fact_type;
    typedef fact_type::ptr_type node_ptr;
    typedef fact_type::const_ptr_type const_node_ptr;

    static node_ptr get_previous(const_node_ptr n) { return n->_previous; }
    static void set_previous(node_ptr n, node_ptr ptr) { n->_previous = ptr; }
    static node_ptr get_next(const_node_ptr n) { return n->_next; }
    static void set_next(node_ptr n, node_ptr ptr) { n->_next = ptr; }
};

struct Read_Chunk_Ptr_Node_List_Value_Traits
{
    typedef Read_Chunk_Ptr_Node value_type;
    typedef Read_Chunk_Ptr_Node_List_Node_Traits node_traits;
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
    Read_Chunk_Ptr_Node,
    boost::intrusive::value_traits< Read_Chunk_Ptr_Node_List_Value_Traits >
> Read_Chunk_Ptr_Node_Cont;

} // namespace MAC


#endif
