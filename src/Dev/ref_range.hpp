#ifndef __REF_RANGE_HPP
#define __REF_RANGE_HPP

#include <boost/range/iterator_range.hpp>
#include "ref_iterator.hpp"


/** Given a range r, construct an equivalent referenced range.
 *  By iterating over ref_range(r), operator& is applied to every value produced.
 */
template < typename Range >
class ref_range
    : public boost::iterator_range< ref_iterator< typename boost::range_iterator< Range >::type > >
{
private:
    typedef boost::iterator_range< ref_iterator< typename boost::range_iterator< Range >::type > > base;
public:
    using typename base::iterator;
    ref_range(Range& r) : base(iterator(std::begin(r)), iterator(std::end(r))) {}
}; // class ref_range

/** Construct a ref_range without the need of template arguments. */
template < typename Range >
ref_range< Range > make_ref_range(Range& r) { return ref_range< Range >(r); }

namespace detail
{
struct ref_range_forwarder {};
}
namespace
{
const detail::ref_range_forwarder referenced = detail::ref_range_forwarder();
}

template < typename Range >
inline ref_range< Range > operator | (Range& r, detail::ref_range_forwarder)
{
    return ref_range< Range >(r);
}


#endif
