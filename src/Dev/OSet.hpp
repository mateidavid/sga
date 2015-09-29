#ifndef __OSET_HPP
#define __OSET_HPP

#include <array>
#include <set>

template < typename T, typename ...Args >
class OSet
    : public std::array< std::set< T, Args... >, 2 >
{
public:
    typedef std::set< T, Args... > Set;
    typedef std::array< Set, 2 > Base;

    using Base::Base;
    using Base::at;
    OSet() = default;

    template < typename Orient_Function >
    OSet(const Set& s, Orient_Function&& f)
    {
        for (auto it = s.begin(); it != s.end(); ++it)
        {
            at(f(*it)).insert(*it);
        }
    }
    template < typename Forward_Iterator, typename Orient_Function >
    OSet(Forward_Iterator it, Forward_Iterator it_end, Orient_Function&& f)
    {
        for ( ; it != it_end; ++it)
        {
            at(f(*it)).insert(*it);
        }
    }

    bool empty() const { return at(0).empty() and at(1).empty(); }
    size_t size() const { return at(0).size() + at(1).size(); }
    void flip(bool do_flip = true) { if (do_flip) swap(at(0), at(1)); }

}; // class OSet

#endif
