#ifndef __FILTER_SET_HPP
#define __FILTER_SET_HPP

#include <functional>

/**
 * Filter container in place; iterator argument version.
 * For those elements where f returns false, call erase_f on an iterator pointing at them.
 * By default, erase_f calls s.erase(it).
 * NOTE: For this function to work, iterators must not be invalidated following a call to erase()!
 * This is suitable for lists/trees, but not for vectors/deques.
 */
template < typename Cont_Type >
void
filter_cont(Cont_Type& s,
            std::function< bool(typename Cont_Type::const_iterator) > f,
            std::function< void(typename Cont_Type::const_iterator) > erase_f = nullptr)
{
    if (not erase_f)
    {
        erase_f = [&] (typename Cont_Type::const_iterator it) { s.erase(it); };
    }
    for (auto it = s.begin(); it != s.end(); )
    {
        auto it_copy = it++;
        if (not f(it_copy))
        {
            erase_f(it_copy);
        }
    }
} // filter_cont

/**
 * Filter container in place; value argument version.
 * For those elements where f returns false, call erase_f on an iterator pointing at them.
 * By default, erase_f calls s.erase(it).
 * NOTE: For this function to work, iterators must not be invalidated following a call to erase()!
 * This is suitable for lists and trees, but not vectors or deques.
 */
template < typename Cont_Type >
void
filter_cont(Cont_Type& s,
            std::function< bool(const typename Cont_Type::value_type&) > f,
            std::function< void(typename Cont_Type::const_iterator) > erase_f = nullptr)
{
    filter_cont(s,
                [&f] (typename Cont_Type::const_iterator cit) { return f(*cit); },
                erase_f);
} // filter_cont

#endif
