#ifndef __RE_DSET_HPP
#define __RE_DSET_HPP

#include "MAC_forward.hpp"
#include "Directed_Set.hpp"

namespace MAC
{

class RC_DSet;

class RE_DSet
    : public Directed_Set< Read_Entry_CBPtr >
{
public:
    typedef Directed_Set< Read_Entry_CBPtr > Base;
    using Base::Base;

    std::unique_ptr< RC_DSet > chunks(Contig_Entry_CBPtr ce_cbptr) const;
}; // class RE_DSet

} // namespace MAC

#endif
