#ifndef __RANGE_HPP
#define __RANGE_HPP

#include <tuple>
#include <utility>
#include "shortcuts.hpp"


namespace range
{

/** A Range is a tuple of coordinates. */
template < typename T >
class Range
{
public:
    DEFAULT_DEF_CTOR(Range);

    Range(const T& start, const T& end) : _start(start), _end(end) {}
    Range(const std::tuple< T, T >& t) : _start(std::get<0>(t)), _end(std::get<1>(t)) {}

    operator std::tuple< T, T >() const { return std::make_tuple(_start, _end); }

    GETTER(T, start, _start)
    GETTER(T, end, _end)
    T len() const { return end() - start(); }

    friend bool operator == (const Range& lhs, const Range& rhs)
    {
        return lhs._start == rhs._start and lhs._end == rhs._end;
    }
    friend bool operator <  (const Range& lhs, const Range& rhs)
    {
        return lhs._start < rhs._start or (lhs._start == rhs._start and lhs._end < rhs._end);
    }
    friend bool operator != (const Range& lhs, const Range& rhs) { return !(lhs == rhs); }
    friend bool operator <= (const Range& lhs, const Range& rhs) { return lhs == rhs or lhs < rhs; }
    friend bool operator >  (const Range& lhs, const Range& rhs) { return !(lhs <= rhs); }
    friend bool operator >= (const Range& lhs, const Range& rhs) { return !(lhs < rhs); }

private:
    T _start;
    T _end;
};

} // namespace Range


#endif
