#ifndef __HAP_ENTRY_LIST_HPP
#define __HAP_ENTRY_LIST_HPP

#include "Hap_Entry.hpp"


namespace MAC
{

namespace detail
{

struct Hap_Entry_List_Node_Traits
{
    typedef Hap_Entry node;
    typedef Hap_Entry_Fact fact_type;
    typedef fact_type::ptr_type node_ptr;
    typedef fact_type::const_ptr_type const_node_ptr;

    static node_ptr get_previous(const_node_ptr n) { return n->_previous; }
    static void set_previous(node_ptr n, node_ptr ptr) { n->_previous = ptr; }
    static node_ptr get_next(const_node_ptr n) { return n->_next; }
    static void set_next(node_ptr n, node_ptr ptr) { n->_next = ptr; }
}; // struct Hap_Entry_List_Node_Traits

struct Hap_Entry_List_Value_Traits
{
    typedef Hap_Entry value_type;
    typedef Hap_Entry_List_Node_Traits node_traits;
    typedef node_traits::node_ptr node_ptr;
    typedef node_traits::const_node_ptr const_node_ptr;
    typedef node_ptr pointer;
    typedef const_node_ptr const_pointer;
    typedef node_traits::fact_type::ref_type reference;
    typedef node_traits::fact_type::const_ref_type const_reference;
    typedef Hap_Entry_List_Value_Traits* value_traits_ptr;

    static const bi::link_mode_type link_mode = bi::safe_link;

    static node_ptr to_node_ptr (reference value) { return &value; }
    static const_node_ptr to_node_ptr (const_reference value) { return &value; }
    static pointer to_value_ptr(node_ptr n) { return n; }
    static const_pointer to_value_ptr(const_node_ptr n) { return n; }
}; // struct Hap_Entry_List_Value_Traits

} // namespace detail

class Hap_Entry_Cont
    : private bi::list< Hap_Entry,
                        bi::value_traits< detail::Hap_Entry_List_Value_Traits >,
                        bi::header_holder_type< bounded::Pointer_Holder< Hap_Entry > >
                      >
{
private:
    typedef bi::list< Hap_Entry,
                      bi::value_traits< detail::Hap_Entry_List_Value_Traits >,
                      bi::header_holder_type< bounded::Pointer_Holder< Hap_Entry > >
                    > Base;
public:
    DEFAULT_DEF_CTOR(Hap_Entry_Cont);
    DELETE_COPY_CTOR(Hap_Entry_Cont);
    DEFAULT_MOVE_CTOR(Hap_Entry_Cont);
    DELETE_COPY_ASOP(Hap_Entry_Cont);
    DEFAULT_MOVE_ASOP(Hap_Entry_Cont);

    USING_INTRUSIVE_CONT(Base)
    using Base::check;
    using Base::clear_and_dispose;

    /// Insert Hap_Entry into container.
    /// @param ce_bptr Pointer to Contig_Entry object.
    void insert(Hap_Entry_BPtr he_bptr)
    {
        Base::push_back(*he_bptr);
    }

    // Erase Hap_Entry
    void erase(Hap_Entry_CBPtr he_cbptr)
    {
        Base::erase(iterator_to(*he_cbptr));
    }

}; // class Hap_Entry_Cont

} // namespace MAC


#endif
