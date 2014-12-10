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
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/tti/has_member_function.hpp>


/** Fix return value from put() and put_child() to allow chaining. */
class ptree
    : public boost::property_tree::ptree
{
public:
    ptree(const std::string& root = std::string(), const boost::property_tree::ptree& t = boost::property_tree::ptree())
      : prefix_(root.empty()? root : root + ".")
    {
        if (not root.empty())
        {
            put_child(root, t);
        }
        else
        {
            *static_cast< boost::property_tree::ptree* >(this) = t;
        }
    }

    template < typename T >
    ptree& put(std::string key, const T& val)
    {
        boost::property_tree::ptree::put(prefix_ + key, val);
        return *this;
    }
    ptree& put(std::string key, const boost::property_tree::ptree& pt)
    {
        boost::property_tree::ptree::put_child(prefix_ + key, pt);
        return *this;
    }
    std::string prefix_;
};


template < bool = true >
std::ostream& operator << (std::ostream& os, const boost::property_tree::ptree& rhs)
{
    boost::property_tree::write_json(os, rhs, false);
    return os;
}

namespace detail
{

// if T has member function to_ptree, call that
template < typename T, bool has_to_ptree >
struct elem_to_ptree_impl
{
    boost::property_tree::ptree operator() (const T& e)
    {
        return e.to_ptree();
    }
};

// if T doesn't have member function to_ptree, try to construct ptree from T
template < typename T >
struct elem_to_ptree_impl< T, false >
{
    boost::property_tree::ptree operator() (const T& e)
    {
        return boost::property_tree::ptree(e);
    }
};

BOOST_TTI_HAS_MEMBER_FUNCTION(to_ptree)

// functor that converts an object of type T to boost::property_tree::ptree
template < typename T >
struct elem_to_ptree : elem_to_ptree_impl< T, has_member_function_to_ptree< const T, boost::property_tree::ptree >::value >
{};

} // namespace detail

template < typename Cont, typename Elem_Type = typename Cont::value_type >
boost::property_tree::ptree cont_to_ptree(const Cont& cont, std::function<boost::property_tree::ptree(const Elem_Type&)> to_ptree)
{
    boost::property_tree::ptree pt;
    auto& array = pt.get_child("");
    for (const typename Cont::value_type& e : cont)
    {
        array.push_back(std::make_pair("", to_ptree(e)));
    }
    return pt;
}

template < typename Cont >
boost::property_tree::ptree cont_to_ptree(const Cont& cont)
{
    typedef typename Cont::const_iterator::value_type elem_type;
    return cont_to_ptree< Cont, elem_type >(cont, detail::elem_to_ptree< elem_type >());
}

#endif
