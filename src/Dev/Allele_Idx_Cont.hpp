#ifndef __ALLELE_IDX_CONT_HPP
#define __ALLELE_IDX_CONT_HPP

#include "Allele_Idx.hpp"

namespace MAC
{

namespace detail
{

struct Allele_Idx_List_Node_Traits
{
    typedef Allele_Idx node;
    typedef Allele_Idx_Fact fact_type;
    typedef fact_type::ptr_type node_ptr;
    typedef fact_type::const_ptr_type const_node_ptr;

    static node_ptr get_previous(const_node_ptr n) { return n->_previous; }
    static void set_previous(node_ptr n, node_ptr ptr) { n->_previous = ptr; }
    static node_ptr get_next(const_node_ptr n) { return n->_next; }
    static void set_next(node_ptr n, node_ptr ptr) { n->_next = ptr; }
}; // struct Allele_Idx_List_Node_Traits

struct Allele_Idx_List_Value_Traits
{
    typedef Allele_Idx value_type;
    typedef Allele_Idx_List_Node_Traits node_traits;
    typedef node_traits::node_ptr node_ptr;
    typedef node_traits::const_node_ptr const_node_ptr;
    typedef node_ptr pointer;
    typedef const_node_ptr const_pointer;
    typedef node_traits::fact_type::ref_type reference;
    typedef node_traits::fact_type::const_ref_type const_reference;
    typedef Allele_Idx_List_Value_Traits* value_traits_ptr;

    static const bi::link_mode_type link_mode = bi::safe_link;

    static node_ptr to_node_ptr (reference value) { return &value; }
    static const_node_ptr to_node_ptr (const_reference value) { return &value; }
    static pointer to_value_ptr(node_ptr n) { return n; }
    static const_pointer to_value_ptr(const_node_ptr n) { return n; }
}; // struct Allele_Idx_List_Value_Traits

} // namespace detail

class Allele_Idx_Cont
    : private bi::list< Allele_Idx,
                        bi::value_traits< detail::Allele_Idx_List_Value_Traits >,
                        bi::constant_time_size< false >,
                        bi::header_holder_type< bounded::Pointer_Holder< Allele_Idx > >
                      >
{
private:
    typedef bi::list< Allele_Idx,
                      bi::value_traits< detail::Allele_Idx_List_Value_Traits >,
                      bi::constant_time_size< false >,
                      bi::header_holder_type< bounded::Pointer_Holder< Allele_Idx > >
                      > Base;
public:
    // allow move only
    DEFAULT_DEF_CTOR(Allele_Idx_Cont);
    DELETE_COPY_CTOR(Allele_Idx_Cont);
    DEFAULT_MOVE_CTOR(Allele_Idx_Cont);
    DELETE_COPY_ASOP(Allele_Idx_Cont);
    DEFAULT_MOVE_ASOP(Allele_Idx_Cont);

    USING_INTRUSIVE_CONT(Base)
    using Base::check;

    // check it is empty when deallocating
    ~Allele_Idx_Cont() { ASSERT(empty()); }

    Base::size_type size() = delete;
    Base::size_type nonconst_size() const { return Base::size(); }

}; // class Allele_Idx_Cont

}

#endif
