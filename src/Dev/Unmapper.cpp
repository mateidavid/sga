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
            _unmap_ce(ce_bptr, extend_re_set);
        }
        // second priority: unmap re-s
        else if (not unmap_re_set.empty())
        {
            auto p = *unmap_re_set.begin();
            unmap_re_set.erase(p);
            _unmap_re_region(p, unmap_ce_set, unmap_re_set);
        }
        // if nothing to unmap, extend unmapped regions
        else
        {
            auto p = *extend_re_set.begin();
            extend_re_set.erase(p);
            _extend_unmappable_re_region(p, unmap_re_set);
        }
    }
}

void
Unmapper::_unmap_ce(Contig_Entry_BPtr ce_bptr,
                    set< pair< Read_Entry_BPtr, Range_Type > >& extend_re_set)
{
    // add all chunks to be unmapped to the extend set
    ce_bptr->chunk_cont().clear_and_dispose(
        [&] (Read_Chunk_BPtr rc_bptr) {
            Range_Type rg(rc_bptr->get_r_start(), rc_bptr->get_r_end());
            extend_re_set.insert(make_pair(rc_bptr->re_bptr(), rg));
            Read_Chunk::make_unmappable(rc_bptr);
            _g.ce_cont().insert(rc_bptr->ce_bptr());
        });
    // done with old Contig_Entry
    ce_bptr->mut_cont().clear_and_dispose();
    _g.ce_cont().erase(ce_bptr);
    Contig_Entry_Fact::del_elem(ce_bptr);
} // Unmapper::_unmap_ce

void
Unmapper::_unmap_re_region(const pair< Read_Entry_BPtr, Range_Type >& p,
                           set< Contig_Entry_BPtr >& unmap_ce_set,
                           set< pair< Read_Entry_BPtr, Range_Type > > unmap_re_set)
{
    Read_Entry_BPtr re_bptr = p.first;
    Size_Type r_start = max(p.second.start(), re_bptr->start());
    Size_Type r_end = min(p.second.end(), re_bptr->end());
    if (r_end <= r_start) return;
    _g.cut_read_entry(re_bptr, r_start);
    _g.cut_read_entry(re_bptr, r_end);
    while (r_start < r_end)
    {
        Read_Chunk_BPtr rc_bptr = re_bptr->chunk_cont().get_chunk_with_pos(r_start).unconst();
        ASSERT(rc_bptr);
        ASSERT(rc_bptr->get_c_start() == 0 and rc_bptr->get_c_end() == rc_bptr->ce_bptr()->len());
        ASSERT(rc_bptr->get_r_end() <= r_end);
        if (not rc_bptr->ce_bptr()->is_unmappable())
        {
            unmap_ce_set.insert(rc_bptr->ce_bptr());
            if (rc_bptr->get_r_end() < r_end)
            {
                // there is a leftover range to unmap beyond this chunk
                Range_Type rg(rc_bptr->get_r_end(), r_end);
                unmap_re_set.insert(make_pair(re_bptr, rg));
            }
            return;
        }
        r_start = rc_bptr->get_r_end();
    }
} // Unmapper::_unmap_re_region

void
Unmapper::_extend_unmappable_re_region(const pair< Read_Entry_BPtr, Range_Type >& p,
                                       set< pair< Read_Entry_BPtr, Range_Type > > unmap_re_set)
{
    Read_Entry_BPtr re_bptr = p.first;
    Read_Chunk_BPtr rc_bptr = re_bptr->chunk_cont().get_chunk_with_pos(p.second.start()).unconst();
    ASSERT(rc_bptr or p.second.end() <= re_bptr->start() or re_bptr->end() <= p.second.start());
    if (not rc_bptr)
    {
        return;
    }
    ASSERT(rc_bptr->get_r_start() <= p.second.start() and p.second.end() <= rc_bptr->get_r_end());
    ASSERT(rc_bptr->ce_bptr()->is_unmappable());
    for (int dir = 0; dir < 2; ++dir)
    {
        Read_Chunk_BPtr next_rc_bptr;
        while (true)
        {
            next_rc_bptr = re_bptr->chunk_cont().get_sibling(rc_bptr, true, dir).unconst();
            if (not next_rc_bptr or not next_rc_bptr->ce_bptr()->is_unmappable())
            {
                break;
            }
            // merge rc and next_rc
            ASSERT(not rc_bptr->get_rc());
            ASSERT(not next_rc_bptr->get_rc());
            bool success = _g.cat_contigs(not dir? rc_bptr->ce_bptr() : next_rc_bptr->ce_bptr(), true);
            static_cast< void >(success);
            ASSERT(success);
            // recompute rc
            rc_bptr = re_bptr->chunk_cont().get_chunk_with_pos(p.second.start()).unconst();
            ASSERT(rc_bptr);
            ASSERT(rc_bptr->get_r_start() <= p.second.start() and p.second.end() <= rc_bptr->get_r_end());
            ASSERT(rc_bptr->ce_bptr()->is_unmappable());
        }
        if (next_rc_bptr)
        {
            // check if there is enough separation to the next unmappable region or read end
            Size_Type separation_len = next_rc_bptr->get_r_len();
            while (separation_len <= _g.unmap_trigger_len())
            {
                next_rc_bptr = re_bptr->chunk_cont().get_sibling(next_rc_bptr, true, dir).unconst();
                if (not next_rc_bptr or next_rc_bptr->ce_bptr()->is_unmappable())
                {
                    break;
                }
                separation_len += next_rc_bptr->get_r_len();
            }
            if (separation_len <= _g.unmap_trigger_len())
            {
                auto rg = (not dir
                           ? Range_Type(rc_bptr->get_r_end(), rc_bptr->get_r_end() + separation_len)
                           : Range_Type(rc_bptr->get_r_start() - separation_len, rc_bptr->get_r_start()));
                unmap_re_set.insert(make_pair(re_bptr, rg));
            }
        }
    } // for dir

    // check if run of unmappable regions extended to the end of the read
    // in this case, trim the read
    ASSERT(rc_bptr);
    ASSERT(rc_bptr->ce_bptr()->is_unmappable());
    if (rc_bptr->get_r_end() == re_bptr->end()
        or rc_bptr->get_r_start() == re_bptr->start())
    {
        _g.trim_terminal_unmappable_chunk(rc_bptr);
    }
} // Unmapper::_extend_unmappable_re_region

} // namespace MAC
