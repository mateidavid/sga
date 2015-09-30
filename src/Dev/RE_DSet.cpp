#include "RE_DSet.hpp"
#include "RC_DSet.hpp"

#include "Contig_Entry.hpp"

namespace MAC
{

std::unique_ptr< RC_DSet > RE_DSet::chunks(Contig_Entry_CBPtr ce_cbptr) const
{
    std::unique_ptr< RC_DSet > res(new RC_DSet());
    for (int d = 0; d < 2; ++d)
    {
        for (auto re_cbptr : at(d))
        {
            auto rc_cbptr = ce_cbptr->chunk_cont().search_read(re_cbptr);
            if (rc_cbptr)
            {
                res->at(d).insert(rc_cbptr);
            }
        }
    }
    return res;
} // RE_DSet::chunks

} // namespace MAC
