//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __PTREE_HPP
#define __PTREE_HPP

#include "shortcuts.hpp" // for is_convertible fix

#include <string>
#include <functional>
#include <type_traits>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "type_traits_extra.hpp"

/** Fix return value from put() and put_child() to allow chaining. */
class ptree
    : public boost::property_tree::ptree
{
public:
    typedef boost::property_tree::ptree base_ptree;
    ptree(const base_ptree& t = base_ptree())
        : ptree(std::string(), t)
    {}

    ptree(const std::string& root, const base_ptree& t = base_ptree())
      : prefix_(root.empty()? root : root + ".")
    {
        if (not root.empty())
        {
            base_ptree::put_child(root, t);
        }
        else
        {
            *static_cast< base_ptree* >(this) = t;
        }
    }

    template < typename T >
    ptree(const std::string& root, const T& t)
      : prefix_(root.empty()? root : root + ".")
    {
        base_ptree::put(root, t);
    }

    template < typename T >
    ptree& put(std::string key, const T& val);

    friend std::ostream& operator << (std::ostream& os, const ptree& rhs)
    {
        boost::property_tree::write_json(os, rhs, false);
        return os;
    }

private:
    std::string prefix_;
}; // class ptree

template < typename Range, typename Printer_Type >
ptree range_to_ptree(const Range& range, Printer_Type&& to_ptree);

template < typename Range >
ptree range_to_ptree(const Range& range);

namespace detail
{
using namespace std_extra;

LOOSE_HAS_MEM_TYPE(const_iterator)
LOOSE_HAS_MEM_FUN(to_ptree)
LOOSE_HAS_MEM_FUN(begin)
LOOSE_HAS_MEM_FUN(end)

template < typename T >
struct is_range_aux
    : and_< has_mem_fun_begin< T, typename T::const_iterator >,
            has_mem_fun_end< T, typename T::const_iterator > > {};

template < typename T >
struct is_range
    : and_< not_< is_convertible< T, string > >,
            has_mem_type_const_iterator< T >,
            is_range_aux< T > > {};

template < typename T >
struct to_ptree_has_implicit_converter
{
    ptree operator () (const T& t) const { return t; }
};

template < typename T >
struct to_ptree_has_explicit_converter
{
    ptree operator() (const T& e) const { return e.to_ptree(); }
};

template < typename T >
struct to_ptree_is_range
{
    ptree operator() (const T& rg) const { return range_to_ptree(rg); }
};

template < typename T >
struct to_ptree_default
{
    ptree operator() (const T& e) const { return ptree("", e); }
};

template < typename T >
struct to_ptree_aux :
        conditional_list< conditional_option< is_convertible< T, boost::property_tree::ptree >,
                                              to_ptree_has_implicit_converter< T > >,
                          conditional_option< has_mem_fun_to_ptree< const T, ptree >,
                                              to_ptree_has_explicit_converter< T > >,
                          conditional_option< is_range< T >,
                                              to_ptree_is_range< T > >,
                          conditional_option< to_ptree_default< T > >
                          >::type {};

/**
 * Functor that converts an object of type T to a ptree
 * Resolution order:
 * 1. If T is implicitly convertible to boost::property_tree::ptree, use this conversion.
 * 2. If T has a to_ptree() method returning a ptree, call it.
 * 3. If T has begin/end methods returning iterators, print it is a range.
 * 4. Else, try a call to put_value with argument T.
 */
template < typename T >
struct to_ptree : to_ptree_aux< typename remove_reference< T >::type > {};

} // namespace detail

template < typename T >
ptree& ptree::put(std::string key, const T& val)
{
    auto pt = detail::to_ptree< T >()(val);
    base_ptree::put_child(prefix_ + key, pt);
    return *this;
}

template < typename Range, typename Printer_Type >
ptree range_to_ptree(const Range& range, Printer_Type&& to_ptree)
{
    ptree pt;
    auto& array = pt.get_child("");
    for (const auto& e : range)
    {
        array.push_back(std::make_pair("", to_ptree(e)));
    }
    return pt;
}

template < typename Range >
ptree range_to_ptree(const Range& range)
{
    return range_to_ptree(range, detail::to_ptree< typename Range::value_type >());
}

#endif
