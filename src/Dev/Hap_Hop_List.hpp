#ifndef __HAP_HOP_LIST_HPP
#define __HAP_HOP_LIST_HPP

#include "Hap_Hop.hpp"


namespace MAC
{

namespace detail
{

struct Hap_Hop_List_Node_Traits
{
    typedef Hap_Hop node;
    typedef Hap_Hop_Fact fact_type;
    typedef fact_type::ptr_type node_ptr;
    typedef fact_type::const_ptr_type const_node_ptr;

    static node_ptr get_previous(const_node_ptr n) { return n->_previous; }
    static void set_previous(node_ptr n, node_ptr ptr) { n->_previous = ptr; }
    static node_ptr get_next(const_node_ptr n) { return n->_next; }
    static void set_next(node_ptr n, node_ptr ptr) { n->_next = ptr; }
}; // struct Hap_Hop_List_Node_Traits

struct Hap_Hop_List_Value_Traits
{
    typedef Hap_Hop value_type;
    typedef Hap_Hop_List_Node_Traits node_traits;
    typedef node_traits::node_ptr node_ptr;
    typedef node_traits::const_node_ptr const_node_ptr;
    typedef node_ptr pointer;
    typedef const_node_ptr const_pointer;
    typedef node_traits::fact_type::ref_type reference;
    typedef node_traits::fact_type::const_ref_type const_reference;
    typedef Hap_Hop_List_Value_Traits* value_traits_ptr;

    static const bi::link_mode_type link_mode = bi::safe_link;

    static node_ptr to_node_ptr (reference value) { return &value; }
    static const_node_ptr to_node_ptr (const_reference value) { return &value; }
    static pointer to_value_ptr(node_ptr n) { return n; }
    static const_pointer to_value_ptr(const_node_ptr n) { return n; }
}; // struct Hap_Hop_List_Value_Traits

} // namespace detail

class Hap_Hop_List
    : private bi::list< Hap_Hop,
                        bi::value_traits< detail::Hap_Hop_List_Value_Traits >,
                        bi::header_holder_type< bounded::Pointer_Holder< Hap_Hop > >
                      >
{
private:
    typedef bi::list< Hap_Hop,
                      bi::value_traits< detail::Hap_Hop_List_Value_Traits >,
                      bi::header_holder_type< bounded::Pointer_Holder< Hap_Hop > >
                    > Base;
public:
    DEFAULT_DEF_CTOR(Hap_Hop_List)
    DELETE_COPY_CTOR(Hap_Hop_List)
    DEFAULT_MOVE_CTOR(Hap_Hop_List)
    DELETE_COPY_ASOP(Hap_Hop_List)
    DEFAULT_MOVE_ASOP(Hap_Hop_List)

    USING_INTRUSIVE_CONT(Base)
    using Base::clear_and_dispose;

    // Get front and back elements
    Hap_Hop_CBPtr front() const { Hap_Hop_CBRef hh_cbref = Base::front(); return &hh_cbref; }
    Hap_Hop_BPtr front() { return &Base::front(); }
    Hap_Hop_CBPtr back() const { return &Base::back(); }
    Hap_Hop_BPtr back() { return &Base::back(); }

    // Insert HH before element pointed to by iterator.
    void insert_before(const_iterator cit, Hap_Hop_BPtr hh_bptr)
    {
        Base::insert(cit, *hh_bptr);
    }
    void insert_before(Hap_Hop_CBPtr existing_hh_cbptr, Hap_Hop_BPtr hh_bptr)
    {
        insert_before(Base::iterator_to(*existing_hh_cbptr), hh_bptr);
    }
    // Insert MCA after element pointed to by iterator.
    void insert_after(const_iterator cit, Hap_Hop_BPtr hh_bptr)
    {
        ASSERT(cit != end());
        Base::insert(++cit, *hh_bptr);
    }
    void insert_after(Hap_Hop_CBPtr existing_hh_cbptr, Hap_Hop_BPtr hh_bptr)
    {
        insert_after(Base::iterator_to(*existing_hh_cbptr), hh_bptr);
    }

    // Insert at the front.
    void push_front(Hap_Hop_BPtr hh_bptr) { insert_before(begin(), hh_bptr); }

    // Insert at the back.
    void push_back(Hap_Hop_BPtr hh_bptr) { insert_before(end(), hh_bptr); }

    // Reverse the list of hops; also flip the c_direction flags
    void reverse()
    {
        for (auto hh_bptr : *this | referenced)
        {
            hh_bptr->c_direction() = not hh_bptr->c_direction();
        }
        Base::reverse();
    }

    // Copy given container into this one; also set the he_cbptr to the given one
    void splice_right(Hap_Hop_List& other_cont, Hap_Entry_CBPtr he_cbptr)
    {
        for (auto hh_bptr : other_cont | referenced)
        {
            hh_bptr->he_cbptr() = he_cbptr;
        }
        Base::splice(end(), static_cast< Base& >(other_cont));
    }

}; // class Hap_Hop_List

} // namespace MAC


#endif
