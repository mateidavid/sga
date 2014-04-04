//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __MUTATION_PTR_CONT_HPP
#define __MUTATION_PTR_CONT_HPP

#include <boost/intrusive/list.hpp>

#include "MAC_forward.hpp"
#include "Mutation.hpp"
#include "Mutation_Chunk_Adapter.hpp"
#include "Mutation_Cont.hpp"


namespace MAC
{

struct Mutation_Ptr_List_Node_Traits
{
    typedef Holder< Mutation_Chunk_Adapter > node;
    typedef Mutation_Chunk_Adapter_Fact fact_type;
    typedef fact_type::ptr_type node_ptr;
    typedef fact_type::const_ptr_type const_node_ptr;

    static node_ptr get_previous(const_node_ptr n) { return n->_mut_ptr_previous; }
    static void set_previous(node_ptr n, node_ptr ptr) { n->_mut_ptr_previous = ptr; }
    static node_ptr get_next(const_node_ptr n) { return n->_mut_ptr_next; }
    static void set_next(node_ptr n, node_ptr ptr) { n->_mut_ptr_next = ptr; }
};

struct Mutation_Ptr_List_Value_Traits
{
    typedef Mutation_Chunk_Adapter value_type;
    typedef Mutation_Ptr_List_Node_Traits node_traits;
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

class Mutation_Ptr_Cont
    : private boost::intrusive::list< Mutation_Chunk_Adapter,
                                      boost::intrusive::value_traits< Mutation_Ptr_List_Value_Traits >
                                    >
{
private:
    typedef boost::intrusive::list< Mutation_Chunk_Adapter,
                                    boost::intrusive::value_traits< Mutation_Ptr_List_Value_Traits >
                                  > Base;
public:
    // allow move only
    DEFAULT_DEF_CTOR(Mutation_Ptr_Cont)
    DELETE_COPY_CTOR(Mutation_Ptr_Cont)
    DEFAULT_MOVE_CTOR(Mutation_Ptr_Cont)
    DELETE_COPY_ASOP(Mutation_Ptr_Cont)
    DEFAULT_MOVE_ASOP(Mutation_Ptr_Cont)

    /*
    typedef boost::iterator::transform_iterator< std::mem_fn(&Mutation_Chunk_Adapter::mut_cbptr),
                                                 Base::iterator
                                               > iterator;
    typedef boost::iterator::transform_iterator< std::mem_fn(&Mutation_Chunk_Adapter::mut_cbptr),
                                                 Base::const_iterator
                                               > const_iterator;
    */

    USING_ITERATORS(Base)

    // check it is empty when deallocating
    ~Mutation_Ptr_Cont() { ASSERT(size() == 0); }

    /** Construct a Mutation_Ptr_Cont from a Mutation_Cont.
     * Pre: Mutations may not touch in reference.
     */
    Mutation_Ptr_Cont(const Mutation_Cont& mut_cont, Read_Chunk_CBPtr chunk_cbptr)
    {
        Mutation_CBPtr last_mut_cbptr = nullptr;
        for (const auto& mut_cbref : mut_cont)
        {
            Mutation_CBPtr mut_cbptr = &mut_cbref;
            if (last_mut_cbptr)
            {
                ASSERT(last_mut_cbptr->get_end() < mut_cbptr->get_start());
            }
            Mutation_Chunk_Adapter_BPtr mca_bptr = Mutation_Chunk_Adapter_Fact::new_elem(mut_cbptr, chunk_cbptr);
            Base::push_back(*mca_bptr);
        }
    }

    /** Insert MCA before element pointed to by iterator. */
    void insert_before(const_iterator cit, Mutation_Chunk_Adapter_BPtr mca_bptr)
    {
        ASSERT(cit != begin());
        Base::insert(cit, *mca_bptr);
    }
    void insert_after(const_iterator cit, Mutation_Chunk_Adapter_BPtr mca_bptr)
    {
        ASSERT(cit != end());
        Base::insert(++cit, *mca_bptr);
    }

    /** Find a given mutation pointer in this container. */
    const_iterator find(Mutation_CBPtr mut_cbptr) const
    {
        for (const auto& mca_cbref : *this)
        {
            if (mca_cbref.raw().mut_cbptr() == mut_cbptr)
            {
                return Base::iterator_to(mca_cbref);
            }
        }
        return end();
    }

    /** Add new MCA if old mutation exists. */
    void cond_insert_after(Mutation_CBPtr old_mut_cbptr, Mutation_Chunk_Adapter_BPtr new_mca_bptr)
    {
        auto cit = find(old_mut_cbptr);
        if (cit != end())
        {
            insert_after(cit, new_mca_bptr);
        }
    }
};

} // namespace MAC


#endif
