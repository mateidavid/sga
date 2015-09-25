#ifndef __TYPE_TRAITS_EXTRA_HPP
#define __TYPE_TRAITS_EXTRA_HPP

#include <type_traits>

namespace std_extra
{
using namespace std;

/** boost::mpl::if_c */
template < bool c, typename T1, typename T2 >
using if_c = conditional< c, T1, T2 >;

/** boost::mpl::if_ */
template < typename C, typename T1, typename T2 >
struct if_ : if_c< C::value, T1, T2 > {};

/** boost::mpl::eval_if_c */
template < bool, typename T1, typename T2 >
struct eval_if_c : common_type< typename T1::type > {};

template < typename T1, typename T2 >
struct eval_if_c< false, T1, T2 > : common_type< typename T2::type > {};

/** boost::mpl::eval_if */
template < typename C, typename T1, typename T2 >
struct eval_if : eval_if_c< C::value, T1, T2 > {};

/** boost::mpl::not_ */
template < typename T >
struct not_ : integral_constant< bool, not T::value > {};

/** boost::mpl::and_ */
template < typename ... >
struct and_ : true_type {};

template < typename T >
struct and_< T >
    : conditional< T::type::value, true_type, false_type >::type {};

template < typename T, typename ...Args >
struct and_< T, Args... >
    : conditional< T::type::value, and_< Args... >, false_type >::type {};

/** boost::mpl::or_ */
template < typename ... >
struct or_ : false_type {};

template < typename T >
struct or_< T >
    : conditional< T::type::value, true_type, false_type >::type {};

template < typename T, typename ...Args >
struct or_< T, Args... >
    : conditional< T::type::value, true_type, or_< Args... > >::type {};

/**
 * If option.
 *
 * Given 2 types C, T, if_option produces a struct R where
 * R::value := C::value, and R::type := T.
 * When a single type is given, C defaults to std::true_type.
 * See also: if_list.
 */
template < typename ... >
struct if_option;

template < typename T >
struct if_option< T > : if_option< std::true_type, T > {};

template < bool c, typename T >
struct if_option_c
{
    static const bool value = c;
    typedef T type;
};

template < typename C, typename T >
struct if_option< C, T > : if_option_c< C::value, T > {};

/**
 * If list.
 *
 * For each template argument type T in order, if T::value is true,
 * the result of if_list is a type R, where R::type := T::type.
 * Crucially, for argument types T2 following T, the inner type T2::type
 * is never instantiated.
 *
 * Together with if_option, this template allows for the definition
 * of a target template type based on a set of conditions, using the following
 * syntax:
 *
 * template < typename T >
 * struct option_1 {...};
 * template < typename T >
 * struct option_2 {...};
 * template < typename T >
 * struct option_default {...};
 * template < typename T >
 * struct target
 *   : if_list<
 *       if_option< condition_1< T >, option_1< T > >,
 *       if_option< condition_2< T >, option_2< T > >,
 *       if_option< option_default< T >
 *     >::type {};
 *
 */
template < typename ... >
struct if_list {};

template < typename T, typename ...Args >
struct if_list< T, Args... >
    : eval_if_c< T::value, T, if_list< Args... > > {};

} // namespace std_extra

/**
 * Test if a type contains a member function with a specific signature.
 *
 * NOTE: As opposed to BOOST_TTI_HAS_MEMBER_FUNCTION, this test:
 * - allows for implicit conversions when matching arguments and return type;
 * - allows for the member function to be defined in a subclass.
 *
 * Template arguments:
 *   1  : Type to test for member function presence (required)
 *   2  : Member function return type (required)
 *   3+ : Member function argument types (optional)
 *
 * Usage:
 *   1. Define test struct using macro:
 *        LOOSE_HAS_MEM_FUN(resize)
 *   2. Use struct:
 *        if(has_mem_fun_resize< Cont, void, size_t >::value)
 */
#define LOOSE_HAS_MEM_FUN(name) \
    template < typename T, typename, typename ...Args > \
    struct has_mem_fun_ ## name ## _aux1 \
    { \
    template < typename ... > \
    static std::false_type test(...); \
    \
    template < typename T2, typename ...Args2 > \
    static auto test(std::nullptr_t) \
        -> decltype(((T2*)nullptr)->name(*(Args2*)nullptr...), std::true_type()); \
    \
    typedef decltype(test< T, Args... >(nullptr)) type; \
    }; \
    \
    template < typename T, typename Return, typename ...Args > \
    struct has_mem_fun_ ## name ## _aux2 \
        : std::is_convertible< decltype(((T*)nullptr)->name(*(Args*)nullptr...)), Return > {}; \
    \
    template < typename T, typename Return, typename ...Args > \
    struct has_mem_fun_ ## name \
        : std_extra::and_< has_mem_fun_ ## name ## _aux1< typename std::remove_reference< T >::type, Return, Args... >, \
                           has_mem_fun_ ## name ## _aux2< typename std::remove_reference< T >::type, Return, Args... > > {};

/**
 * Test if a type contains a member type.
 *
 * This test allows for the member type to be defined in a subclass.
 * NOTE: As opposed to BOOST_TTI_HAS_TYPE, this test:
 * - removes all references before the test.
 *
 * Template arguments:
 *   1  : Type to test for member function presence (required)
 *
 * Usage:
 *   1. Define test struct using macro:
 *        LOOSE_HAS_MEM_TYPE(iterator)
 *   2. Use struct:
 *        if(has_mem_type_iterator< Cont >::value)
 */
#define LOOSE_HAS_MEM_TYPE(name) \
    template < typename T > \
    struct has_mem_type_ ## name ## _aux \
    { \
    template < typename ... > \
    static std::false_type test(...); \
    \
    template < typename T2 > \
    static auto test(std::nullptr_t) \
        -> decltype(*(typename T2::name*)nullptr, std::true_type()); \
    \
    typedef decltype(test< T >(nullptr)) type; \
    }; \
    \
    template < typename T > \
    struct has_mem_type_ ## name \
        : has_mem_type_ ## name ## _aux< typename std::remove_reference< T >::type >::type {};

#endif
