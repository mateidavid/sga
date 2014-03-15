#ifndef __NONCONST_METHODS_HPP
#define __NONCONST_METHODS_HPP

#include <utility>


namespace detail
{
    /** Type function which adds constness to references and pointers. */
    template < typename T >
    struct const_referred_type {};
    template < typename T >
    struct const_referred_type< T& >
    {
        typedef const T& type;
    };
    template < typename T >
    struct const_referred_type< T* >
    {
        typedef const T* type;
    };
    template < typename T >
    struct const_referred_type< T&& >
    {
        typedef const T&& type;
    };

    /** Extract type from parantheses; to be used with macros. */
    template < typename T >
    struct argument_type;
    template < typename T >
    struct argument_type< void(T) >
    {
        typedef T type;
    };

    template <typename Return, class T, typename... Args>
    Return nonconst_method(
        const T* obj,
        typename const_referred_type< Return >::type (T::* member_fun)(Args...) const,
        Args... args)
    {
        return const_cast< Return >((obj->*member_fun)(args...));
    }

    template <typename Return, class T>
    Return nonconst_conversion(const T* obj)
    {
        return const_cast< Return >(obj->operator typename const_referred_type< Return >::type ());
    }
}

using detail::nonconst_method;
using detail::nonconst_conversion;

#define DEF_NONCONST_METHOD_T(_type, _function) \
    template < typename... Args > \
    typename detail::argument_type< void(_type) >::type _function(Args&&... args) \
    { \
        return const_cast< typename argument_type< void(_type) >::type >( \
            const_cast< typename detail::const_referred_type< decltype(this) >::type >(this) \
                ->_function(std::forward< Args >(args)...) ); \
    }

#define DEF_NONCONST_CONVERSION(_type) \
    operator typename detail::argument_type< void(_type) >::type () \
    { \
        return const_cast< typename argument_type< void(_type) >::type >( \
            const_cast< typename detail::const_referred_type< decltype(this) >::type >(this) \
                ->operator typename detail::const_referred_type< typename detail::argument_type< void(_type) >::type >::type ()); \
    }


#endif
