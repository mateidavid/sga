#include "RC_DSet.hpp"
#include "RE_DSet.hpp"

#include "Read_Chunk.hpp"

namespace MAC
{

unique_ptr< RE_DSet > RC_DSet::reads() const
{
    unique_ptr< RE_DSet > res(new RE_DSet);
    for (int d = 0; d < 2; ++d)
    {
        for (auto rc_cbptr: this->at(d))
        {
            res->at(d).insert(rc_cbptr->re_bptr());
        }
    }
    return res;
} // RC_DSet::reads

} // namespace MAC
