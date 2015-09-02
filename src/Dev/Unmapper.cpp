#include "Unmapper.hpp"

namespace MAC
{

void
Unmapper::_unmap_loop(set< Contig_Entry_BPtr >& unmap_ce_set,
                      set< pair< Read_Entry_BPtr, Range_Type > >& unmap_re_set)
{
    set< pair< Read_Entry_BPtr, Range_Type > > extend_re_set;
    while (not unmap_ce_set.empty()
           or not unmap_re_set.empty()
           or not extend_re_set.empty())
    {
        // first priority: unmap ce-s
        if (not unmap_ce_set.empty())
        {
            Contig_Entry_BPtr ce_bptr = *unmap_ce_set.begin();
            unmap_ce_set.erase(ce_bptr);
            _unmap_contig(ce_bptr, extend_re_set);
        }
        // second priority: unmap re-s
        else if (not unmap_re_set.empty())
        {
            auto p = *unmap_re_set.begin();
            unmap_re_set.erase(p);
            _unmap_region(p, unmap_ce_set, unmap_re_set);
        }
        // if nothing to unmap, extend unmapped regions
        else
        {
            auto p = *extend_re_set.begin();
            extend_re_set.erase(p);
            _extend_region(p, unmap_ce_set);
        }
    }
}

void
Unmapper::_unmap_ce(Contig_Entry_BPtr ce_bptr,
                    set< pair< Read_Entry_BPtr, Range_Type > >& extend_re_set)
{
}

void
Unmapper::_unmap_re_region(const pair< Read_Entry_BPtr, Range_Type >& p,
                           set< Contig_Entry_BPtr >& unmap_ce_set,
                           set< pair< Read_Entry_BPtr, Range_Type > > unmap_re_set)
{
}

void
Unmapper::_extend_unmappable_re_region(const pair< Read_Entry_BPtr, Range_Type >& p,
                                       set< Contig_Entry_BPtr >& unmap_ce_set)
{
}

} // namespace MAC
