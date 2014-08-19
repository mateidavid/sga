//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __PTREE_HPP
#define __PTREE_HPP

#include "shortcuts.hpp" // for is_convertible fix

#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>


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
