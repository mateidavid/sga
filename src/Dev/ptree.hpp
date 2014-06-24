//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __PTREE_TOOLS_HPP
#define __PTREE_TOOLS_HPP

#include "shortcuts.hpp" // for is_convertible fix

#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>


/** Fix return value from put() and put_child() to allow chaining. */
class ptree
    : public boost::property_tree::ptree
{
public:
    template < typename T >
    ptree& put(std::string key, const T& val)
    {
        boost::property_tree::ptree::put(key, val);
        return *this;
    }
    ptree& put(std::string key, const boost::property_tree::ptree& pt)
    {
        boost::property_tree::ptree::put_child(key, pt);
        return *this;
    }
};

template < bool = true >
std::ostream& operator << (std::ostream& os, const boost::property_tree::ptree& rhs)
{
    boost::property_tree::write_json(os, rhs, false);
    return os;
}

template < typename Cont >
boost::property_tree::ptree cont_to_ptree(const Cont& cont)
{
    boost::property_tree::ptree pt;
    auto& array = pt.get_child("");
    for (const auto& e : cont)
    {
        array.push_back(std::make_pair("", e.to_ptree()));
    }
    return pt;
}

#endif
