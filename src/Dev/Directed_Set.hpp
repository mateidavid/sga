#ifndef __DIRECTED_SET_HPP
#define __DIRECTED_SET_HPP

#include <array>
#include <set>
#include "filter_cont.hpp"

template < typename T, typename ...Args >
class Directed_Set
    : public std::array< std::set< T, Args... >, 2 >
{
public:
    typedef std::set< T, Args... > Set;
    typedef std::array< Set, 2 > Base;

    using Base::Base;
    using Base::at;
    Directed_Set() = default;

    template < typename Input_Iterator, typename Orient_Function >
    Directed_Set(Input_Iterator it, Input_Iterator it_end, Orient_Function&& f)
    {
        for ( ; it != it_end; ++it)
        {
            at(f(*it)).insert(*it);
        }
    }
    template < typename Input_Range, typename Orient_Function >
    Directed_Set(Input_Range&& r, Orient_Function&& f)
        : Directed_Set(r.begin(), r.end(), f) {}

    bool empty() const { return at(0).empty() and at(1).empty(); }
    size_t size() const { return at(0).size() + at(1).size(); }
    void reverse(bool do_reverse = true) { if (do_reverse) swap(at(0), at(1)); }
    Set collapse() const
    {
        Set s;
        s.insert(at(0).begin(), at(0).end());
        s.insert(at(1).begin(), at(1).end());
        return s;
    }

    template < typename Filter_Function >
    void filter(Filter_Function&& f)
    {
        for (int d = 0; d < 2; ++d)
        {
            filter_cont(at(d), [&] (const T& v) { return f(v, bool(d)); });
        }
    }

    void intersect(const Directed_Set& other, bool rel_direction)
    {
        filter([&] (const T& v, bool d) { return other.at((d + rel_direction) % 2).count(v) > 0; });
    }

    void subtract(const Directed_Set& other, bool rel_direction)
    {
        filter([&] (const T& v, bool d) { return not other.at((d + rel_direction) % 2).count(v); });
    }

    void add(const Directed_Set& other, bool rel_direction)
    {
        for (int d = 0; d < 2; ++d)
        {
            at(d).insert(other.at((d + rel_direction) % 2).begin(), other.at((d + rel_direction) % 2).end());
        }
    }

    static std::pair< Directed_Set, Directed_Set >
    separate(const Directed_Set& s1, const Directed_Set& s2, bool rel_direction)
    {
        std::pair< Directed_Set, Directed_Set > res;
        for (int d = 0; d < 2; ++d)
        {
            for (auto re_cbptr : s1.at(d))
            {
                if (s2.at((d + rel_direction) % 2).count(re_cbptr))
                {
                    res.first.at(d).insert(re_cbptr);
                }
                else
                {
                    res.second.at(d).insert(re_cbptr);
                }
            }
        }
        return res;
    }

    static Directed_Set
    intersect(const Directed_Set& s1, const Directed_Set& s2, bool rel_direction)
    { return separate(s1, s2, rel_direction).first; }

    static Directed_Set
    subtract(const Directed_Set& s1, const Directed_Set& s2, bool rel_direction)
    { return separate(s1, s2, rel_direction).second; }

}; // class Directed_Set

#endif
