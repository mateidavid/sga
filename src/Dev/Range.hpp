#ifndef __RANGE_HPP
#define __RANGE_HPP

#include <tuple>
#include <utility>
#include "shortcuts.hpp"
#include "ptree.hpp"


namespace range
{

/** A Range is a tuple of coordinates. */
template < typename T >
class Range
    : public std::pair< T, T >
{
public:
    typedef std::pair< T, T > Base;
    using Base::Base;

    const T& begin() const { return this->first; }
    T& begin() { return this->first; }
    const T& end() const { return this->second; }
    T& end() { return this->second; }
    T len() const { return end() - begin(); }

    bool empty() const { return end() <= begin(); }
    bool touch(const Range& other) const
    {
        return (begin() <= other.begin() and other.begin() <= end())
            or (other.begin() <= begin() and begin() <= other.end());
    }
    Range& contract(const Range& other)
    {
        begin() = std::max(begin(), other.begin());
        end() = std::min(end(), other.end());
        return *this;
    }
    Range& extend(const Range& other)
    {
        begin() = std::min(begin(), other.begin());
        end() = std::max(end(), other.end());
        return *this;
    }

    ptree to_ptree() const
    {
        return ptree().put("begin", begin()).put("end", end());
    }
};

} // namespace Range


#endif
