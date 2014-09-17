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

namespace detail
{

/** MCA Cloner: copy Read_Chunk pointers, use new Mutation pointer. */
class MCA_Read_Chunk_Ptr_Cloner
{
public:
    MCA_Read_Chunk_Ptr_Cloner(Mutation_CBPtr mut_cbptr) : _mut_cbptr(mut_cbptr) {}

    Mutation_Chunk_Adapter_BPtr operator () (const Mutation_Chunk_Adapter& mca)
    {
        return Mutation_Chunk_Adapter_Fact::new_elem(_mut_cbptr, mca.chunk_cbptr());
    }

private:
    const Mutation_CBPtr _mut_cbptr;
};

struct Read_Chunk_Ptr_List_Node_Traits
{
    typedef Mutation_Chunk_Adapter node;
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
    typedef Read_Chunk_Ptr_List_Value_Traits* value_traits_ptr;

    static const bi::link_mode_type link_mode = bi::safe_link;

    static node_ptr to_node_ptr (reference value) { return &value; }
    static const_node_ptr to_node_ptr (const_reference value) { return &value; }
    static pointer to_value_ptr(node_ptr n) { return n; }
    static const_pointer to_value_ptr(const_node_ptr n) { return n; }
};

} // namespace detail

class Read_Chunk_Ptr_Cont
    : private bi::list< Mutation_Chunk_Adapter,
                        bi::value_traits< detail::Read_Chunk_Ptr_List_Value_Traits >,
                        bi::constant_time_size< false >,
                        bi::header_holder_type< bounded::Pointer_Holder< Mutation_Chunk_Adapter > >
                      >
{
private:
    typedef bi::list< Mutation_Chunk_Adapter,
                      bi::value_traits< detail::Read_Chunk_Ptr_List_Value_Traits >,
                      bi::constant_time_size< false >,
                      bi::header_holder_type< bounded::Pointer_Holder< Mutation_Chunk_Adapter > >
                    > Base;
    typedef detail::MCA_Read_Chunk_Ptr_Cloner Cloner;
public:
    // allow move only
    DEFAULT_DEF_CTOR(Read_Chunk_Ptr_Cont)
    DELETE_COPY_CTOR(Read_Chunk_Ptr_Cont)
    DEFAULT_MOVE_CTOR(Read_Chunk_Ptr_Cont)
    DELETE_COPY_ASOP(Read_Chunk_Ptr_Cont)
    DEFAULT_MOVE_ASOP(Read_Chunk_Ptr_Cont)

    USING_INTRUSIVE_CONT(Base)
    using Base::splice;
    using Base::check;

    // check it is empty when deallocating
    ~Read_Chunk_Ptr_Cont() { ASSERT(empty()); }

    Base::size_type size() = delete;
    Base::size_type nonconst_size() const { return Base::size(); }

    /** Insert read chunk in this container. */
    void insert(Mutation_Chunk_Adapter_BPtr mca_bptr)
    {
        Base::push_back(*mca_bptr);
    }

    /** Clone from another Read_Chunk_Ptr_Cont.
     * @param src Source container.
     * @param new_mut_cbptr New Mutation pointer to use.
     */
    void clone_from(const Read_Chunk_Ptr_Cont& src, Mutation_CBPtr new_mut_cbptr)
    {
        Base::clone_from(src, Read_Chunk_Ptr_Cont::Cloner(new_mut_cbptr), Mutation_Chunk_Adapter_Fact::disposer_type());
    }

    /** Erase MCA from container. */
    void erase(Mutation_Chunk_Adapter_CBPtr mca_cbptr)
    {
        Base::erase(iterator_to(*mca_cbptr));
    }
};

} // namespace MAC

#endif
