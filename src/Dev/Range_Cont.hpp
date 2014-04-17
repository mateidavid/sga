//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __RANGE_CONT_HPP
#define __RANGE_CONT_HPP

#include <tuple>
#include <set>
#include "shortcuts.hpp"


/** A range is a tuple of coordinates. */
template < typename T >
using Range = std::tuple< T, T >;

/** A range container is a set that stores disjoint ranges. */
template < typename T >
class Range_Cont
    : private std::set< Range< T > >
{
private:
    typedef std::set< Range< T > > Base;
public:
    typedef T value_type;
    typedef Range< T > range_type;

    // allow move but not copy
    DEFAULT_DEF_CTOR(Range_Cont)
    DELETE_COPY_CTOR(Range_Cont)
    DEFAULT_MOVE_CTOR(Range_Cont)
    DELETE_COPY_ASOP(Range_Cont)
    DEFAULT_MOVE_ASOP(Range_Cont)

    USING_ITERATORS(Base)

    /** If the range overlaps or touches and existing range,
     * extend that range; else insert this one.
     */
    void insert(const range_type& range)
    {
        // ignore negative ranges
        if (std::get<1>(range) < std::get<0>(range))
        {
            return;
        }
        // find first range that could overlap this one
        auto it = Base::lower_bound(range);
        if (it != begin())
        {
            --it;
        }
        while (it != end() and std::get<0>(*it) <= std::get<1>(range))
        {
            if (std::get<0>(range) <= std::get<1>(*it))
            {
                // overlap:
                // either get<0>(range) <  get<0>(*it)   [ <= get<1>(range) ]
                // or     get<0>(*it)   <= get<0>(range) [ <= get<1>(*it)   ]
                range_type new_range = std::make_tuple(std::min(std::get<0>(*it), std::get<0>(range)),
                                                       std::max(std::get<0>(*it), std::get<0>(range)));
                Base::erase(it);
                Base::insert(new_range);
                return;
            }
            ++it;
        }
        // no overlaps were found; just insert new range
        Base::insert(range);
    }
};


#endif
