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
    static_assert(boost::intrusive::detail::has_member<Base, has_mem_get_parent>::value, "Node Traits missing get_parent()");
    static_assert(boost::intrusive::detail::has_member<Base, has_mem_set_parent>::value, "Node Traits missing set_parent()");
    static_assert(boost::intrusive::detail::has_member<Base, has_mem_get_left>::value, "Node Traits missing get_left()");
    static_assert(boost::intrusive::detail::has_member<Base, has_mem_set_left>::value, "Node Traits missing set_left()");
    static_assert(boost::intrusive::detail::has_member<Base, has_mem_get_right>::value, "Node Traits missing get_right()");
    static_assert(boost::intrusive::detail::has_member<Base, has_mem_set_right>::value, "Node Traits missing set_right()");
    static_assert(boost::intrusive::detail::has_member<Base, has_mem_get_color>::value, "Node Traits missing get_color()");
    static_assert(boost::intrusive::detail::has_member<Base, has_mem_set_color>::value, "Node Traits missing set_color()");
    static_assert(boost::intrusive::detail::has_member<Base, has_mem_black>::value, "Node Traits missing black()");
    static_assert(boost::intrusive::detail::has_member<Base, has_mem_red>::value, "Node Traits missing red()");
    static_assert(boost::intrusive::detail::has_member<Base, has_mem_get_max_end>::value, "Node Traits missing get_max_end()");
    static_assert(boost::intrusive::detail::has_member<Base, has_mem_set_max_end>::value, "Node Traits missing set_max_end()");
public:
    typedef ITree_Node_Traits node_traits;
    typedef typename Base::node node;
    typedef typename Base::node_ptr node_ptr;
    typedef typename Base::const_node_ptr const_node_ptr;
    typedef typename Base::color color;
    typedef typename Base::key_type key_type;

    static void init_data(node_ptr n)
    {
        Base::set_max_end(n, Value_Traits::get_end(Value_Traits::to_value_ptr(n)));
    }
    static void recompute_data(node_ptr n)
    {
        init_data(n);
        key_type tmp = Base::get_max_end(n);
        if (Base::get_left(n))
        {
            tmp = std::max(tmp, Base::get_max_end(Base::get_left(n)));
        }
        if (Base::get_right(n))
        {
            tmp = std::max(tmp, Base::get_max_end(Base::get_right(n)));
        }
        Base::set_max_end(n, tmp);
    }
    static void copy_data(node_ptr dest, const_node_ptr src)
    {
        Base::set_max_end(dest, Base::get_max_end(src));
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
    static_assert(boost::intrusive::detail::has_member<Base, has_mem_to_node_ptr>::value, "Node Traits missing to_node_ptr()");
    static_assert(boost::intrusive::detail::has_member<Base, has_mem_to_node_ptr_const>::value, "Node Traits missing to_node_ptr() const");
    static_assert(boost::intrusive::detail::has_member<Base, has_mem_to_value_ptr>::value, "Node Traits missing to_value_ptr()");
    static_assert(boost::intrusive::detail::has_member<Base, has_mem_to_value_ptr_const>::value, "Node Traits missing to_value_ptr() const");
    static_assert(boost::intrusive::detail::has_member<Base, has_mem_get_start_raw>::value, "Node Traits missing get_start() raw");
    static_assert(boost::intrusive::detail::has_member<Base, has_mem_get_start_wrap>::value, "Node Traits missing get_start() wrap");
    static_assert(boost::intrusive::detail::has_member<Base, has_mem_get_end_raw>::value, "Node Traits missing get_end() raw");
    static_assert(boost::intrusive::detail::has_member<Base, has_mem_get_end_wrap>::value, "Node Traits missing get_end() wrap");
public:
    typedef ITree_Node_Traits< Value_Traits > node_traits;
    typedef typename Base::value_type value_type;
    typedef typename Base::node_ptr node_ptr;
    typedef typename Base::const_node_ptr const_node_ptr;
    typedef typename Base::pointer pointer;
    typedef typename Base::const_pointer const_pointer;
    typedef typename Base::reference reference;
    typedef typename Base::const_reference const_reference;
};

/** Comparator for Interval Tree.
 *
 */
template <typename Value_Traits>
struct ITree_Compare
{
    bool operator () (const typename Value_Traits::value_type& lhs, const typename Value_Traits::value_type& rhs) const
    {
        return Value_Traits::get_start(&lhs) < Value_Traits::get_start(&rhs);
    }
};

}


template <typename Value_Traits>
class itree
  : public boost::intrusive::multiset<
        typename Value_Traits::value_type,
        boost::intrusive::compare< detail::ITree_Compare< Value_Traits > >,
        boost::intrusive::value_traits< detail::ITree_Value_Traits< Value_Traits > >
    >
{
public:
    typedef boost::intrusive::multiset<
        typename Value_Traits::value_type,
        boost::intrusive::compare<
            detail::ITree_Compare< Value_Traits >
        >,
        boost::intrusive::value_traits<
            detail::ITree_Value_Traits< Value_Traits >
        >
    > Base;

    template <typename... Args>
    itree(Args&&... args) : Base(std::forward<Args>(args)...) {}
};


#endif
