#ifndef __NONCONST_METHODS_HPP
#define __NONCONST_METHODS_HPP

#include <utility>


namespace nonconst
{
    /** Type function which adds constness to references and pointers. */
    template < typename T >
    struct const_referred_type;
    template < typename T >
    struct const_referred_type< T& > { typedef const T& type; };
    template < typename T >
    struct const_referred_type< T* > { typedef const T* type; };
    template < typename T >
    struct const_referred_type< T&& > { typedef const T&& type; };

    /** Extract type from parantheses: allows passing type as macro argument. */
    template < typename T >
    struct argument_type;
    template < typename T >
    struct argument_type< void(T) > { typedef T type; };

    /*
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

    template < typename T >
    const T* const_this(T* ptr)
    {
        return const_cast< const T* >(ptr);
    }
    */
}

/*
using nonconst::nonconst_method;
using nonconst::nonconst_conversion;
using nonconst::const_this;
*/

/* Macro that defines the non-const version of a const method.
 *
 * Arguments:
 *   _type : Return type of non-const version, possibly paranthesized.
 *           Note, as macro argument, the type must be paranthesized
 *           if it contains commas.
 *   _function : Name of non-const method.
 *
 * NOTE: This macro can only be used when the return type of non-const _function
 * can be obtained from the return type of const _function by const_cast<>.
 * In particular, this works if the return type of const _function is a
 * const-reference or const-pointer, and the return type of non-const _function
 * is the corresponding nonconst-reference or nonconst-pointer.
 *
 * Example:
 *   class A {
 *     const T<B,C>& at(size_t i) const { return _t[i]; }
 *     NONCONST_METHOD((T<B,C>&), at)
 *     // produces:
 *     //   template < typename... Args > T<B,C>& at(Args&&... args)
 *     //   {
 *     //       return const_cast< T<B,C>& >(
 *     //           const_cast< const A* >(this)
 *     //               ->at(std::forward< Args >(args)...));
 *     //   }
 *     const U& operator[] (size_t i) const { return _u[i]; }
 *     NONCONST_METHOD(U&, operator[])
 *     // produces:
 *     //   template < typename... Args > U& operator[](Args&&... args)
 *     //   {
 *     //       return const_cast< U& >(
 *     //           const_cast< const A* >(this)
 *     //               ->operator[](std::forward< Args >(args)...));
 *     //   }
 *   }
 */
#define NONCONST_METHOD(_type, _function) \
    template < typename... Args > \
    typename nonconst::argument_type< void(_type) >::type _function(Args&&... args) \
    { \
        return const_cast< typename nonconst::argument_type< void(_type) >::type >( \
            const_cast< typename nonconst::const_referred_type< decltype(this) >::type >(this) \
                ->_function(std::forward< Args >(args)...) ); \
    }

/* Macro that defines the non-const version of a conversion operator.
 *
 * Similar to NONCONST_METHOD. Only works when the return type of the const operator
 * is a const-reference or const-pointer, and the return type of the nonconst operator
 * is the corresponding nonconst-reference or nonconst-pointer.
 */
#define NONCONST_CONVERSION(_type) \
    operator typename nonconst::argument_type< void(_type) >::type () \
    { \
        return const_cast< typename nonconst::argument_type< void(_type) >::type >( \
            const_cast< typename nonconst::const_referred_type< decltype(this) >::type >(this) \
                ->operator typename nonconst::const_referred_type< typename nonconst::argument_type< void(_type) >::type >::type ()); \
    }


#endif
