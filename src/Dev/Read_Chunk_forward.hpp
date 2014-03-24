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

//
// Read_Chunk pointers are stored in intrusive lists inside Mutation objects
//
struct Read_Chunk_Ptr_Node_List_Node_Traits;

class Read_Chunk_Ptr_Node
{
private:
    friend class Factory< Read_Chunk_Ptr_Node >;

    Read_Chunk_Ptr_Node(Read_Chunk_BPtr rc_bptr = 0) : _rc_bptr(rc_bptr) {}

    // allow move only
    DELETE_COPY_CTOR(Read_Chunk_Ptr_Node)
    DEFAULT_MOVE_CTOR(Read_Chunk_Ptr_Node)
public:
    DELETE_COPY_ASOP(Read_Chunk_Ptr_Node)
    DEFAULT_MOVE_ASOP(Read_Chunk_Ptr_Node)

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
