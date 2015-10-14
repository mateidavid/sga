#ifndef __RC_DSET_HPP
#define __RC_DSET_HPP

#include "MAC_forward.hpp"
#include "Directed_Set.hpp"
#include "RE_DSet.hpp"

namespace MAC
{

class RE_DSet;

class RC_DSet
    : public Directed_Set< Read_Chunk_CBPtr >
{
public:
    typedef Directed_Set< Read_Chunk_CBPtr > Base;
    using Base::Base;

    RE_DSet reads() const;
}; // class RC_DSet

} // namespace MAC

#endif
