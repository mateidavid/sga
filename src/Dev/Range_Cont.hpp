//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __RANGE_CONT_HPP
#define __RANGE_CONT_HPP

#include <tuple>
#include <set>
#include "Range.hpp"
#include "shortcuts.hpp"


namespace range
{

/** A range container is a set that stores disjoint ranges. */
template < typename T >
class Range_Cont
    : private std::set< Range< T > >
{
private:
    typedef std::set< Range< T > > Base;
public:
    typedef Range< T > range_type;

    // allow move but not copy
    DEFAULT_DEF_CTOR(Range_Cont);
    DELETE_COPY_CTOR(Range_Cont);
    DEFAULT_MOVE_CTOR(Range_Cont);
    DELETE_COPY_ASOP(Range_Cont);
    DEFAULT_MOVE_ASOP(Range_Cont);

    USING_STD_CONT(Base)
    using Base::erase;

    /** If the range overlaps or touches and existing range,
     * extend that range; else insert this one.
     */
    void insert(range_type range)
    {
        // ignore negative ranges
        if (range.end() <= range.start())
        {
            return;
        }
        // find all ranges that could overlap this one, and merge them into it
        auto it = Base::lower_bound(range);
        if (it != begin())
        {
            --it;
        }
        while (it != end() and it->start() <= range.end())
        {
            auto it_copy = it++;
            if (range.start() <= it_copy->end())
            {
                // overlap:
                // either range.start <  it->start   <= range.end
                // or     it->start   <= range.start <= it->end
                range = range_type(std::min(it_copy->start(), range.start()),
                                   std::max(it_copy->end(), range.end()));
                Base::erase(it_copy);
            }
        }
        // done extending, insert new range
        Base::insert(range);
        check();
    }
private:
    void check()
    {
#ifndef DISABLE_ASSERTS
        auto it_last = begin();
        if (it_last == end()) return;
        auto it = std::next(it_last);
        while (it != end())
        {
            ASSERT(it_last->end() < it->start());
            it_last = it++;
        }
#endif
    }
};

} // namespace Range

#endif
