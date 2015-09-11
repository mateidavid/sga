#ifndef __MIN_MAX_OF_HPP
#define __MIN_MAX_OF_HPP

#include <algorithm>
#include <type_traits>

template < class ForwardIterator, class UnaryFunction >
ForwardIterator min_of(ForwardIterator first, ForwardIterator last, UnaryFunction f)
{
    if (first == last) return last;
    auto min_it = first++;
    for (auto it = first; it != last; ++it)
    {
        if (f(*it) < f(*min_it))
        {
            min_it = it;
        }
    }
    return min_it;
}

template < class ForwardRange, class UnaryFunction >
typename ForwardRange::iterator min_of(const ForwardRange& rg, UnaryFunction f)
{
    return min_of(rg.begin(), rg.end(), f);
}

template < class ForwardIterator, class UnaryFunction >
typename std::result_of< UnaryFunction(typename ForwardIterator::value_type) >::type
min_value_of(ForwardIterator first, ForwardIterator last, UnaryFunction f)
{
    auto it = min_of(first, last, f);
    return it != last
        ? f(*it)
        : std::result_of< UnaryFunction(ForwardIterator) >::type();
}

template < class ForwardRange, class UnaryFunction >
typename std::result_of< UnaryFunction(typename ForwardRange::iterator::value_type) >::type
min_value_of(const ForwardRange& rg, UnaryFunction f)
{
    return min_value_of(rg.begin(), rg.end(), f);
}

template < class ForwardIterator, class UnaryFunction >
ForwardIterator max_of(ForwardIterator first, ForwardIterator last, UnaryFunction f)
{
    if (first == last) return last;
    auto max_it = first++;
    for (auto it = first; it != last; ++it)
    {
        if (f(*it) > f(*max_it))
        {
            max_it = it;
        }
    }
    return max_it;
}

template < class ForwardRange, class UnaryFunction >
typename ForwardRange::iterator max_of(const ForwardRange& rg, UnaryFunction f)
{
    return max_of(rg.begin(), rg.end(), f);
}

template < class ForwardIterator, class UnaryFunction >
typename std::result_of< UnaryFunction(typename ForwardIterator::value_type) >::type
max_value_of(ForwardIterator first, ForwardIterator last, UnaryFunction f)
{
    auto it = max_of(first, last, f);
    return it != last
        ? f(*it)
        : std::result_of< UnaryFunction(ForwardIterator) >::type();
}

template < class ForwardRange, class UnaryFunction >
typename std::result_of< UnaryFunction(typename ForwardRange::iterator::value_type) >::type
max_value_of(const ForwardRange& rg, UnaryFunction f)
{
    return max_value_of(rg.begin(), rg.end(), f);
}

template < class ForwardIterator, class UnaryFunction >
std::pair< ForwardIterator, ForwardIterator > minmax_of(ForwardIterator first, ForwardIterator last, UnaryFunction f)
{
    if (first == last) return make_pair(last, last);
    auto min_it = first++;
    auto max_it = min_it;
    for (auto it = first; it != last; ++it)
    {
        if (f(*it) < f(*min_it))
        {
            min_it = it;
        }
        if (f(*it) > f(*max_it))
        {
            max_it = it;
        }
    }
    return make_pair(min_it, max_it);
}

template < class ForwardRange, class UnaryFunction >
std::pair< typename ForwardRange::iterator, typename ForwardRange::iterator >
minmax_of(const ForwardRange& rg, UnaryFunction f)
{
    return minmax_of(rg.begin(), rg.end(), f);
}

template < class ForwardIterator, class UnaryFunction >
std::pair< typename std::result_of< UnaryFunction(typename ForwardIterator::value_type) >::type,
           typename std::result_of< UnaryFunction(typename ForwardIterator::value_type) >::type >
minmax_value_of(ForwardIterator first, ForwardIterator last, UnaryFunction f)
{
    auto p = minmax_of(first, last, f);
    return p.first != last
        ? make_pair(f(*first), f(*last))
        : make_pair(typename std::result_of< UnaryFunction(typename ForwardIterator::value_type) >::type(),
                    typename std::result_of< UnaryFunction(typename ForwardIterator::value_type) >::type());
}

template < class ForwardRange, class UnaryFunction >
std::pair< typename std::result_of< UnaryFunction(typename ForwardIterator::value_type) >::type,
           typename std::result_of< UnaryFunction(typename ForwardIterator::value_type) >::type >
minmax_value_of(const ForwardRange& rg, UnaryFunction f)
{
    return min_value_of(rg.begin(), rg.end(), f);
}

#endif
