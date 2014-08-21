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

namespace detail
{

struct Contig_Entry_List_Node_Traits
{
    typedef Contig_Entry node;
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
    typedef Contig_Entry_List_Value_Traits* value_traits_ptr;

    static const bi::link_mode_type link_mode = bi::safe_link;

    static node_ptr to_node_ptr (reference value) { return &value; }
    static const_node_ptr to_node_ptr (const_reference value) { return &value; }
    static pointer to_value_ptr(node_ptr n) { return n; }
    static const_pointer to_value_ptr(const_node_ptr n) { return n; }
};

} // namespace detail

class Contig_Entry_Cont
    : private bi::list< Contig_Entry,
                        bi::value_traits< detail::Contig_Entry_List_Value_Traits >,
                        bi::header_holder_type< bounded::Pointer_Holder< Contig_Entry > >
                      >
{
private:
    typedef bi::list< Contig_Entry,
                      bi::value_traits< detail::Contig_Entry_List_Value_Traits >,
                      bi::header_holder_type< bounded::Pointer_Holder< Contig_Entry > >
                    > Base;
public:
    DEFAULT_DEF_CTOR(Contig_Entry_Cont)
    DELETE_COPY_CTOR(Contig_Entry_Cont)
    DEFAULT_MOVE_CTOR(Contig_Entry_Cont)
    DELETE_COPY_ASOP(Contig_Entry_Cont)
    DEFAULT_MOVE_ASOP(Contig_Entry_Cont)
    // check it is empty when deallocating
    ~Contig_Entry_Cont() { ASSERT(size() == 0); }

    USING_INTRUSIVE_CONT(Base)
    friend class Graph;

    /** Insert Contig_Entry into container.
     * @param ce_bptr Pointer to Contig_Entry object.
     */
    void insert(Contig_Entry_BPtr ce_bptr)
    {
        Base::push_back(*ce_bptr);
    }

    /** Erase Contig_Entry from container. */
    void erase(Contig_Entry_CBPtr ce_cbptr)
    {
        Base::erase(iterator_to(*ce_cbptr));
    }

    /** Clear container and deallocate CE objects.
     * For each CE, its chunks are first removed from their RE containers and
     * they are deallocated.
     */
    void clear_and_dispose()
    {
        Base::clear_and_dispose([] (Contig_Entry_BPtr ce_bptr)
        {
            // remove chunks from their RE cont
            ce_bptr->chunk_cont().erase_from_re_cont();
            // deallocate chunks and mca-s
            ce_bptr->chunk_cont().clear_and_dispose();
            // deallocate mutations
            ce_bptr->mut_cont().clear_and_dispose();
            // deallocate CE
            Contig_Entry_Fact::del_elem(ce_bptr);
        });
    }
};

} // namespace MAC


#endif
