#ifndef __REF_ITERATOR_HPP
#define __REF_ITERATOR_HPP

#include <boost/iterator/transform_iterator.hpp>
#include "shortcuts.hpp"


/** Iterator adaptor class that applies operator& to produced elements.
 */
template < typename Iterator >
class ref_iterator
    : public boost::iterator_adaptor< ref_iterator< Iterator >,
                                      Iterator,
                                      typename Iterator::pointer,
                                      boost::use_default,
                                      typename Iterator::pointer >
{
public:
    ref_iterator()
        : ref_iterator::iterator_adaptor_(0) {}
    ref_iterator(const Iterator& other)
        : ref_iterator::iterator_adaptor_(other) {}

    template < class Other_Iterator,
               T_ENABLE_IF((boost::is_convertible< Other_Iterator, Iterator >::value)) >
    ref_iterator(const ref_iterator< Other_Iterator >& other)
        : ref_iterator::iterator_adaptor_(other.base()) {}

private:
    friend class boost::iterator_core_access;
    typename Iterator::pointer dereference() const { return &*(this->base()); }
}; // class ref_iterator


#endif
