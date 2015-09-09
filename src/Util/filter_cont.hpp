#ifndef __FILTER_SET_HPP
#define __FILTER_SET_HPP

#include <functional>

template < typename Cont_Type >
void filter_cont(Cont_Type& s, std::function< bool(typename Cont_Type::const_iterator) > f)
{
    for (auto it = s.begin(); it != s.end(); )
    {
        auto it_copy = it++;
        if (not f(it_copy))
        {
            s.erase(it_copy);
        }
    }
} // filter_cont

template < typename Cont_Type >
void filter_cont(Cont_Type& s, std::function< bool(const typename Cont_Type::value_type&) > f)
{
    filter_cont(s, [&f] (typename Cont_Type::const_iterator cit) { return f(*cit); });
} // filter_cont

#endif
