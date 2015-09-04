#include "Unmapper.hpp"

namespace MAC
{

void
Unmapper::_unmap_loop(ce_set_type&& unmap_ce_set, re_set_type&& unmap_re_set)
{
    re_set_type extend_re_set;
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
            auto& p = *unmap_re_set.begin();
            Read_Entry_BPtr re_bptr = p.first;
            ASSERT(not unmap_re_set.at(re_bptr).empty());
            Range_Type rg = *unmap_re_set.at(re_bptr).begin();
            unmap_re_set.at(re_bptr).erase(rg);
            if (unmap_re_set.at(re_bptr).empty())
            {
                unmap_re_set.erase(re_bptr);
            }
            _unmap_re_region(re_bptr, rg, unmap_ce_set, unmap_re_set);
        }
        // if nothing to unmap, extend unmapped regions
        else
        {
            auto& p = *extend_re_set.begin();
            Read_Entry_BPtr re_bptr = p.first;
            ASSERT(not extend_re_set.at(re_bptr).empty());
            Range_Type rg = *extend_re_set.at(re_bptr).begin();
            extend_re_set.at(re_bptr).erase(rg);
            if (extend_re_set.at(re_bptr).empty())
            {
                extend_re_set.erase(re_bptr);
            }
            _extend_unmappable_re_region(re_bptr, rg, unmap_re_set);
        }
    }
}

void
Unmapper::_unmap_ce(Contig_Entry_BPtr ce_bptr, re_set_type& extend_re_set)
{
    LOG("Unmapper", debug) << ptree("_unmap_ce")
        .put("ce_bptr", ce_bptr.to_int());
    // add all chunks to be unmapped to the extend set
    set< Read_Entry_CBPtr > re_to_check;
    ce_bptr->chunk_cont().clear_and_dispose(
        [&] (Read_Chunk_BPtr rc_bptr) {
            Range_Type rg(rc_bptr->get_r_start(), rc_bptr->get_r_end());
            extend_re_set[rc_bptr->re_bptr()].insert(rg);
            re_to_check.insert(rc_bptr->re_bptr());
            Read_Chunk::make_unmappable(rc_bptr);
            _g.ce_cont().insert(rc_bptr->ce_bptr());
        });
    // done with old Contig_Entry
    ce_bptr->mut_cont().clear_and_dispose();
    _g.ce_cont().erase(ce_bptr);
    Contig_Entry_Fact::del_elem(ce_bptr);
    _g.check(re_to_check);
} // Unmapper::_unmap_ce

void
Unmapper::_unmap_re_region(Read_Entry_BPtr re_bptr, const Range_Type& rg,
                           ce_set_type& unmap_ce_set, re_set_type&)
{
    LOG("Unmapper", debug) << ptree("_unmap_re_region")
        .put("re_bptr", re_bptr.to_int())
        .put("rg_start", rg.start())
        .put("rg_end", rg.end());
    Size_Type r_start = max(rg.start(), re_bptr->start());
    Size_Type r_end = min(rg.end(), re_bptr->end());
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
        }
        r_start = rc_bptr->get_r_end();
    }
} // Unmapper::_unmap_re_region

void
Unmapper::_extend_unmappable_re_region(Read_Entry_BPtr re_bptr, const Range_Type& rg,
                                       re_set_type& unmap_re_set)
{
    LOG("Unmapper", debug) << ptree("_extend_unmappable_re_region")
        .put("re_bptr", re_bptr.to_int())
        .put("rg_start", rg.start())
        .put("rg_end", rg.end());
    Read_Chunk_BPtr rc_bptr = re_bptr->chunk_cont().get_chunk_with_pos(rg.start()).unconst();
    ASSERT(rc_bptr or rg.end() <= re_bptr->start() or re_bptr->end() <= rg.start());
    if (not rc_bptr)
    {
        return;
    }
    ASSERT(rc_bptr->ce_bptr()->is_unmappable());
    for (int dir = 0; dir < 2; ++dir)
    {
        Read_Chunk_BPtr next_rc_bptr;
        while (true)
        {
            next_rc_bptr = re_bptr->chunk_cont().get_sibling(rc_bptr, true, not dir).unconst();
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
            rc_bptr = re_bptr->chunk_cont().get_chunk_with_pos(rg.start()).unconst();
            ASSERT(rc_bptr);
            ASSERT(rc_bptr->ce_bptr()->is_unmappable());
        }
        if (next_rc_bptr)
        {
            // check if there is enough separation to the next unmappable region or read end
            Size_Type separation_len = next_rc_bptr->get_r_len();
            while (separation_len <= _g.unmap_trigger_len())
            {
                next_rc_bptr = re_bptr->chunk_cont().get_sibling(next_rc_bptr, true, not dir).unconst();
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
                unmap_re_set[re_bptr].insert(rg);
            }
        }
    } // for dir

    // check if run of unmappable regions extended to the end of the read
    // in this case, trim the read
    if (_g.trim_during_unmapping())
    {
        ASSERT(rc_bptr);
        ASSERT(rc_bptr->ce_bptr()->is_unmappable());
        if (rc_bptr->get_r_end() == re_bptr->end()
            or rc_bptr->get_r_start() == re_bptr->start())
        {
            _g.trim_terminal_unmappable_chunk(rc_bptr);
        }
    }
    _g.check({re_bptr});
} // Unmapper::_extend_unmappable_re_region

} // namespace MAC
