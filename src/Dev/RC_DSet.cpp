#include "RC_DSet.hpp"

#include "Read_Chunk.hpp"

namespace MAC
{

RE_DSet RC_DSet::reads() const
{
    RE_DSet res;
    for (int d = 0; d < 2; ++d)
    {
        for (auto rc_cbptr: this->at(d))
        {
            res.at(d).insert(rc_cbptr->re_bptr());
        }
    }
    return res;
} // RC_DSet::reads

} // namespace MAC
