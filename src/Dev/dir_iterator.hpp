#ifndef __DIR_ITERATOR_HPP
#define __DIR_ITERATOR_HPP

#include <iterator>
#include <boost/iterator/iterator_adaptor.hpp>
#include "shortcuts.hpp"


template < typename Iterator >
class dir_iterator
    : public boost::iterator_adaptor< dir_iterator< Iterator >, Iterator >
{
private:
    typedef boost::iterator_adaptor< dir_iterator< Iterator >, Iterator > Base;
public:
    dir_iterator()
        : dir_iterator::iterator_adaptor_(0), _reversed(false) {}
    dir_iterator(const Iterator& other, bool reversed = false)
        : dir_iterator::iterator_adaptor_(other), _reversed(reversed) {}

    template < class Other_Iterator,
               T_ENABLE_IF((boost::is_convertible< Other_Iterator, Iterator >::value)) >
    dir_iterator(const dir_iterator< Other_Iterator >& other)
        : dir_iterator::iterator_adaptor_(other.base()), _reversed(other._reversed) {}

private:
    friend class boost::iterator_core_access;

    typename Iterator::reference dereference() const
    {
        if (not _reversed)
        {
            return *(this->base());
        }
        else
        {
            Iterator tmp = this->base();
            return *(--tmp);
        }
    }
    void increment() { increment_dir(true); }
    void decrement() { increment_dir(false); }
    void increment_dir(bool fwd)
    {
        if (fwd == not _reversed)
        {
            ++this->base_reference();
        }
        else
        {
            --this->base_reference();
        }
    }

    bool _reversed;
};


#endif
