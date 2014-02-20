#ifndef __ITREE_HPP
#define __ITREE_HPP

#include <boost/intrusive/set.hpp>


namespace detail
{

// Check for Node Traits members
struct has_mem_get_parent
{
    template <typename T, typename T::node_ptr (*)(typename T::const_node_ptr) = &T::get_parent> struct get {};
};
struct has_mem_set_parent
{
    template <typename T, void (*)(typename T::node_ptr, typename T::node_ptr) = &T::set_parent> struct get {};
};
struct has_mem_get_left
{
    template <typename T, typename T::node_ptr (*)(typename T::const_node_ptr) = &T::get_left> struct get {};
};
struct has_mem_set_left
{
    template <typename T, void (*)(typename T::node_ptr, typename T::node_ptr) = &T::set_left> struct get {};
};
struct has_mem_get_right
{
    template <typename T, typename T::node_ptr (*)(typename T::const_node_ptr) = &T::get_right> struct get {};
};
struct has_mem_set_right
{
    template <typename T, void (*)(typename T::node_ptr, typename T::node_ptr) = &T::set_right> struct get {};
};
struct has_mem_get_color
{
    template <typename T, typename T::color (*)(typename T::const_node_ptr) = &T::get_color> struct get {};
};
struct has_mem_set_color
{
    template <typename T, void (*)(typename T::node_ptr, typename T::color) = &T::set_color> struct get {};
};
struct has_mem_black
{
    template <typename T, typename T::color (*)() = &T::black> struct get {};
};
struct has_mem_red
{
    template <typename T, typename T::color (*)() = &T::red> struct get {};
};
struct has_mem_get_max_end
{
    template <typename T, typename T::key_type (*)(typename T::const_node_ptr) = &T::get_max_end> struct get {};
};
struct has_mem_set_max_end
{
    template <typename T, void (*)(typename T::node_ptr, typename T::key_type) = &T::set_max_end> struct get {};
};

// Check for Value_Traits members
struct has_mem_to_node_ptr
{
    template <typename T, typename T::node_ptr (*)(typename T::reference) = &T::to_node_ptr> struct get {};
};
struct has_mem_to_node_ptr_const
{
    template <typename T, typename T::const_node_ptr (*)(typename T::const_reference) = &T::to_node_ptr> struct get {};
};
struct has_mem_to_value_ptr
{
    template <typename T, typename T::pointer (*)(typename T::node_ptr) = &T::to_value_ptr> struct get {};
};
struct has_mem_to_value_ptr_const
{
    template <typename T, typename T::const_pointer (*)(typename T::const_node_ptr) = &T::to_value_ptr> struct get {};
};
struct has_mem_get_start_raw
{
    template <typename T, typename T::key_type (*)(const typename T::value_type*) = &T::get_start> struct get {};
};
struct has_mem_get_start_wrap
{
    template <typename T, typename T::key_type (*)(typename T::const_pointer) = &T::get_start> struct get {};
};
struct has_mem_get_end_raw
{
    template <typename T, typename T::key_type (*)(const typename T::value_type*) = &T::get_end> struct get {};
};
struct has_mem_get_end_wrap
{
    template <typename T, typename T::key_type (*)(typename T::const_pointer) = &T::get_end> struct get {};
};

using boost::intrusive::detail::has_member;

/** Node Traits adaptor for Interval Tree.
 *
 * This Traits class defines the node maintenance methods that hook into
 * the rbtree algorithms.
 */
template <typename Value_Traits>
struct ITree_Node_Traits : public Value_Traits::node_traits
{
private:
    typedef typename Value_Traits::node_traits Base;
    static_assert(has_member<Base, has_mem_get_parent>::value, "Node Traits missing get_parent()");
    static_assert(has_member<Base, has_mem_set_parent>::value, "Node Traits missing set_parent()");
    static_assert(has_member<Base, has_mem_get_left>::value, "Node Traits missing get_left()");
    static_assert(has_member<Base, has_mem_set_left>::value, "Node Traits missing set_left()");
    static_assert(has_member<Base, has_mem_get_right>::value, "Node Traits missing get_right()");
    static_assert(has_member<Base, has_mem_set_right>::value, "Node Traits missing set_right()");
    static_assert(has_member<Base, has_mem_get_color>::value, "Node Traits missing get_color()");
    static_assert(has_member<Base, has_mem_set_color>::value, "Node Traits missing set_color()");
    static_assert(has_member<Base, has_mem_black>::value, "Node Traits missing black()");
    static_assert(has_member<Base, has_mem_red>::value, "Node Traits missing red()");
    static_assert(has_member<Base, has_mem_get_max_end>::value, "Node Traits missing get_max_end()");
    static_assert(has_member<Base, has_mem_set_max_end>::value, "Node Traits missing set_max_end()");
public:
    typedef ITree_Node_Traits node_traits;
    using typename Base::node;
    using typename Base::node_ptr;
    using typename Base::const_node_ptr;
    using typename Base::color;
    using typename Base::key_type;
    using Base::get_parent;
    using Base::set_parent;
    using Base::get_left;
    using Base::set_left;
    using Base::get_right;
    using Base::set_right;
    using Base::get_color;
    using Base::set_color;
    using Base::black;
    using Base::red;
    using Base::get_max_end;
    using Base::set_max_end;

    static void init_data(node_ptr n)
    {
        set_max_end(n, Value_Traits::get_end(Value_Traits::to_value_ptr(n)));
    }
    static void recompute_data(node_ptr n)
    {
        init_data(n);
        key_type tmp = get_max_end(n);
        if (get_left(n))
        {
            tmp = std::max(tmp, get_max_end(get_left(n)));
        }
        if (get_right(n))
        {
            tmp = std::max(tmp, get_max_end(get_right(n)));
        }
        set_max_end(n, tmp);
    }
    static void copy_data(node_ptr dest, const_node_ptr src)
    {
        set_max_end(dest, get_max_end(src));
    }
};

/** Value Traits adaptor class for Interval Tree.
 *
 * The only function is to change the node_traits typedef.
 */
template <typename Value_Traits>
struct ITree_Value_Traits : public Value_Traits
{
private:
    typedef Value_Traits Base;
    static_assert(has_member<Base, has_mem_to_node_ptr>::value, "Node Traits missing to_node_ptr()");
    static_assert(has_member<Base, has_mem_to_node_ptr_const>::value, "Node Traits missing to_node_ptr() const");
    static_assert(has_member<Base, has_mem_to_value_ptr>::value, "Node Traits missing to_value_ptr()");
    static_assert(has_member<Base, has_mem_to_value_ptr_const>::value, "Node Traits missing to_value_ptr() const");
    static_assert(has_member<Base, has_mem_get_start_raw>::value, "Node Traits missing get_start() raw");
    static_assert(has_member<Base, has_mem_get_start_wrap>::value, "Node Traits missing get_start() wrap");
    static_assert(has_member<Base, has_mem_get_end_raw>::value, "Node Traits missing get_end() raw");
    static_assert(has_member<Base, has_mem_get_end_wrap>::value, "Node Traits missing get_end() wrap");
public:
    typedef ITree_Node_Traits< Value_Traits > node_traits;
    using typename Base::value_type;
    using typename Base::node_ptr;
    using typename Base::const_node_ptr;
    using typename Base::pointer;
    using typename Base::const_pointer;
    using typename Base::reference;
    using typename Base::const_reference;
    using Base::to_node_ptr;
    using Base::to_value_ptr;
    using Base::get_start;
    using Base::get_end;
};

/** Comparator for Interval Tree.
 *
 */
template <typename Value_Traits>
struct ITree_Compare
{
    typedef typename Value_Traits::value_type value_type;
    bool operator () (const value_type& lhs, const value_type& rhs) const
    {
        return Value_Traits::get_start(&lhs) < Value_Traits::get_start(&rhs);
    }
};


template <typename Value_Traits>
class itree
  : public boost::intrusive::multiset<
        typename Value_Traits::value_type,
        boost::intrusive::compare< ITree_Compare< Value_Traits > >,
        boost::intrusive::value_traits< ITree_Value_Traits< Value_Traits > >
    >
{
public:
    typedef boost::intrusive::multiset<
        typename Value_Traits::value_type,
        boost::intrusive::compare< ITree_Compare< Value_Traits > >,
        boost::intrusive::value_traits< ITree_Value_Traits< Value_Traits > >
    > Base;

    // inherit multiset constructors
    using Base::Base;
};

}

using detail::itree;


#endif
