#ifndef __TYPE_STR
#define __TYPE_STR


#include <type_traits>
#include <typeinfo>
#include <string>


template < typename T >
struct type_str
{
    operator std::string () const
    {
        std::string res;
        // referenced type
        if (std::is_lvalue_reference< T >::value or std::is_rvalue_reference< T >::value)
        {
            res = type_str< typename std::remove_reference< T >::type >();
            res += " ";
        }
        else if (std::is_pointer< T >::value)
        {
            res = type_str< typename std::remove_pointer< T >::type >();
            res += " ";
        }
        // cv-qualifiers
        if (std::is_const< T >::value)
        {
            res += "const ";
        }
        if (std::is_volatile< T >::value)
        {
            res += "volatile ";
        }
        // main type
        if (std::is_lvalue_reference< T >::value)
        {
            res += "&";
        }
        else if (std::is_rvalue_reference< T >::value)
        {
            res += "&&";
        }
        else if (std::is_pointer< T >::value)
        {
            res += "*";
        }
        else
        {
            res += typeid(typename std::remove_cv< T >::type).name();
        }
        return res;
    }

    friend std::ostream& operator << (std::ostream& os, type_str rhs)
    {
        os << static_cast< std::string >(rhs);
        return os;
    }
}; // struct type_str


#endif
