#ifndef __HAP_HOP_SET_HPP
#define __HAP_HOP_SET_HPP

#include "Hap_Hop.hpp"


namespace MAC
{

namespace detail
{

struct Hap_Hop_Set_Node_Traits
{
    typedef Hap_Hop node;
    typedef Hap_Hop_Fact fact_type;
    typedef fact_type::ptr_type node_ptr;
    typedef fact_type::const_ptr_type const_node_ptr;
    typedef bool color;

    static node_ptr get_parent(const_node_ptr n) { return n->_parent; }
    static void set_parent(node_ptr n, node_ptr ptr) { n->_parent = ptr; }
    static node_ptr get_left(const_node_ptr n) { return n->_l_child; }
    static void set_left(node_ptr n, node_ptr ptr) { n->_l_child = ptr; }
    static node_ptr get_right(const_node_ptr n) { return n->_r_child; }
    static void set_right(node_ptr n, node_ptr ptr) { n->_r_child = ptr; }
    static color get_color(const_node_ptr n) { return n->_col; }
    static void set_color(node_ptr n, color c) { n->_col = c ; }
    static bool black() { return 0; }
    static bool red() { return 1; }
}; // struct Hap_Hop_Set_Node_Traits

struct Hap_Hop_Set_Value_Traits
{
    typedef Hap_Hop value_type;
    typedef Hap_Hop_Set_Node_Traits node_traits;
    typedef node_traits::node_ptr node_ptr;
    typedef node_traits::const_node_ptr const_node_ptr;
    typedef node_ptr pointer;
    typedef const_node_ptr const_pointer;
    typedef node_traits::fact_type::ref_type reference;
    typedef node_traits::fact_type::const_ref_type const_reference;
    typedef Hap_Hop_Set_Value_Traits* value_traits_ptr;

    static const bi::link_mode_type link_mode = bi::safe_link;

    static node_ptr to_node_ptr (reference value) { return &value; }
    static const_node_ptr to_node_ptr (const_reference value) { return &value; }
    static pointer to_value_ptr(node_ptr n) { return n; }
    static const_pointer to_value_ptr(const_node_ptr n) { return n; }
}; // struct Hap_Hop_Set_Value_Traits

/** Comparator for storage in tree. */
struct Hap_Hop_Comparator
{
    bool operator () (const Hap_Hop& lhs, const Hap_Hop& rhs) const
    {
        return lhs.allele_anchor() < rhs.allele_anchor();
    }
    bool operator () (const Hap_Hop& lhs, const Allele_Anchor& rhs) const
    {
        return lhs.allele_anchor() < rhs;
    }
    bool operator () (const Allele_Anchor& lhs, const Hap_Hop& rhs) const
    {
        return lhs < rhs.allele_anchor();
    }
}; // struct Hap_Hop_Comparator

} // namespace detail

class Hap_Hop_Set
    : private bi::multiset< Hap_Hop,
                            bi::compare< detail::Hap_Hop_Comparator >,
                            bi::value_traits< detail::Hap_Hop_Set_Value_Traits >,
                            bi::header_holder_type< bounded::Pointer_Holder< Hap_Hop > >
                          >
{
private:
    typedef bi::multiset< Hap_Hop,
                          bi::compare< detail::Hap_Hop_Comparator >,
                          bi::value_traits< detail::Hap_Hop_Set_Value_Traits >,
                          bi::header_holder_type< bounded::Pointer_Holder< Hap_Hop > >
                        > Base;
public:
    DEFAULT_DEF_CTOR(Hap_Hop_Set);
    DELETE_COPY_CTOR(Hap_Hop_Set);
    DEFAULT_MOVE_CTOR(Hap_Hop_Set);
    DELETE_COPY_ASOP(Hap_Hop_Set);
    DEFAULT_MOVE_ASOP(Hap_Hop_Set);

    USING_INTRUSIVE_CONT(Base)
    using Base::check;

    // Insert Hap_Hop into this container
    void insert(Hap_Hop_BPtr hh_bptr)
    {
        Base::insert(*hh_bptr);
    }
    // Erase Hap_Hop from this container
    void erase(Hap_Hop_BPtr hh_bptr)
    {
        Base::erase(iterator_to(*hh_bptr));
    }

    // Get hops at given anchor.
    pair< const_iterator, const_iterator > equal_range(const Allele_Anchor& anchor) const
    {
        return Base::equal_range(anchor, Base::value_comp());
    }
    pair< iterator, iterator > equal_range(const Allele_Anchor& anchor)
    {
        return Base::equal_range(anchor, Base::value_comp());
    }

    // Get hops at anchor, grouped by allele
    typedef map< Allele_Specifier, set< pair< Hap_Hop_CBPtr, bool > > > find_anchor_type;
    find_anchor_type find_anchor(const Allele_Anchor& anchor, bool c_direction) const
    {
        find_anchor_type res;
        auto p = equal_range(anchor);
        for (auto it = p.first; it != p.second; ++it)
        {
            auto hh_bptr = &*it;
            ASSERT(hh_bptr->allele_anchor() == anchor);
            res[hh_bptr->allele_specifier()].insert(make_pair(hh_bptr, hh_bptr->c_direction() != c_direction));
        }
        return res;
    }
    // Get hops at allele
    typedef find_anchor_type::mapped_type find_allele_type;
    find_allele_type find_allele(const Allele_Anchor& anchor, const Allele_Specifier& allele, bool c_direction) const
    {
        return move(find_anchor(anchor, c_direction).at(allele));
    }

}; // class Hap_Hop_Set

} // namespace MAC


#endif
