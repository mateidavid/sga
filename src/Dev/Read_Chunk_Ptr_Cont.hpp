//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __READ_CHUNK_PTR_CONT_HPP
#define __READ_CHUNK_PTR_CONT_HPP

#include <boost/intrusive/list.hpp>

#include "Mutation_Chunk_Adapter.hpp"


namespace MAC
{

struct Read_Chunk_Ptr_List_Node_Traits
{
    typedef Holder< Mutation_Chunk_Adapter > node;
    typedef Mutation_Chunk_Adapter_Fact fact_type;
    typedef fact_type::ptr_type node_ptr;
    typedef fact_type::const_ptr_type const_node_ptr;

    static node_ptr get_previous(const_node_ptr n) { return n->_chunk_ptr_previous; }
    static void set_previous(node_ptr n, node_ptr ptr) { n->_chunk_ptr_previous = ptr; }
    static node_ptr get_next(const_node_ptr n) { return n->_chunk_ptr_next; }
    static void set_next(node_ptr n, node_ptr ptr) { n->_chunk_ptr_next = ptr; }
};

struct Read_Chunk_Ptr_List_Value_Traits
{
    typedef Mutation_Chunk_Adapter value_type;
    typedef Read_Chunk_Ptr_List_Node_Traits node_traits;
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

class Read_Chunk_Ptr_Cont
    : private boost::intrusive::list< Mutation_Chunk_Adapter,
                                      boost::intrusive::value_traits< Read_Chunk_Ptr_List_Value_Traits >
                                    >
{
private:
    typedef boost::intrusive::list< Mutation_Chunk_Adapter,
                                    boost::intrusive::value_traits< Read_Chunk_Ptr_List_Value_Traits >
                                  > Base;
public:
    // allow move only
    DEFAULT_DEF_CTOR(Read_Chunk_Ptr_Cont)
    DELETE_COPY_CTOR(Read_Chunk_Ptr_Cont)
    DEFAULT_MOVE_CTOR(Read_Chunk_Ptr_Cont)
    DELETE_COPY_ASOP(Read_Chunk_Ptr_Cont)
    DEFAULT_MOVE_ASOP(Read_Chunk_Ptr_Cont)

    USING_ITERATORS(Base)
    using Base::splice;

    // check it is empty when deallocating
    ~Read_Chunk_Ptr_Cont() { ASSERT(size() == 0); }

    /** Insert read chunk in this container. */
    void insert(Mutation_Chunk_Adapter_BPtr mca_bptr)
    {
        Base::push_back(*mca_bptr);
    }
};

} // namespace MAC

#endif
