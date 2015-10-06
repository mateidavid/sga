#include "Unmapper.hpp"
#include "Graph.hpp"

namespace MAC
{

void
Unmapper::_unmap_loop(ce_set_type&& unmap_ce_set, re_set_type&& unmap_re_set) const
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
Unmapper::_unmap_ce(Contig_Entry_BPtr ce_bptr, re_set_type& extend_re_set) const
{
    LOG("Unmapper", debug) << ptree("begin")
        .put("ce_bptr", ce_bptr.to_int());
    // add all chunks to be unmapped to the extend set
    set< Read_Entry_CBPtr > re_to_check;
    ce_bptr->chunk_cont().clear_and_dispose(
        [&] (Read_Chunk_BPtr rc_bptr) {
            Range_Type rg(rc_bptr->get_r_start(), rc_bptr->get_r_end());
            extend_re_set[rc_bptr->re_bptr()].insert(rg);
            re_to_check.insert(rc_bptr->re_bptr());
            Read_Chunk::make_unmappable(rc_bptr);
            _g_p->ce_cont().insert(rc_bptr->ce_bptr());
        });
    // done with old Contig_Entry
    ce_bptr->mut_cont().clear_and_dispose();
    _g_p->ce_cont().erase(ce_bptr);
    Contig_Entry_Fact::del_elem(ce_bptr);
    _g_p->check(re_to_check);
    LOG("Unmapper", debug) << ptree("end");
} // Unmapper::_unmap_ce

void
Unmapper::_unmap_re_region(Read_Entry_BPtr re_bptr, Range_Type rg,
                           ce_set_type& unmap_ce_set, re_set_type&) const
{
    LOG("Unmapper", debug) << ptree("begin")
        .put("re_bptr", re_bptr.to_int())
        .put("rg", rg);
    rg.contract(Range_Type(re_bptr->start(), re_bptr->end()));
    if (rg.empty()) return;
    _g_p->cut_read_entry(re_bptr, rg.begin());
    _g_p->cut_read_entry(re_bptr, rg.end());
    while (not rg.empty())
    {
        Read_Chunk_BPtr rc_bptr = re_bptr->chunk_cont().get_chunk_with_pos(rg.begin()).unconst();
        ASSERT(rc_bptr);
        ASSERT(rc_bptr->get_c_start() == 0 and rc_bptr->get_c_end() == rc_bptr->ce_bptr()->len());
        ASSERT(rc_bptr->get_r_end() <= rg.end());
        if (not rc_bptr->ce_bptr()->is_unmappable())
        {
            unmap_ce_set.insert(rc_bptr->ce_bptr());
        }
        rg.begin() = rc_bptr->get_r_end();
    }
    LOG("Unmapper", debug) << ptree("end");
} // Unmapper::_unmap_re_region

void
Unmapper::_extend_unmappable_re_region(Read_Entry_BPtr re_bptr, Range_Type rg,
                                       re_set_type& unmap_re_set) const
{
    LOG("Unmapper", debug) << ptree("begin")
        .put("re_bptr", re_bptr.to_int())
        .put("rg", rg);
    rg.contract(Range_Type(re_bptr->start(), re_bptr->end()));
    if (rg.empty()) return;
    Read_Chunk_BPtr rc_bptr = re_bptr->chunk_cont().get_chunk_with_pos(rg.begin()).unconst();
    ASSERT(rc_bptr);
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
            bool success = _g_p->cat_contigs(not dir? rc_bptr->ce_bptr() : next_rc_bptr->ce_bptr(), true);
            static_cast< void >(success);
            ASSERT(success);
            // recompute rc
            rc_bptr = re_bptr->chunk_cont().get_chunk_with_pos(rg.begin()).unconst();
            ASSERT(rc_bptr);
            ASSERT(rc_bptr->ce_bptr()->is_unmappable());
        }
        if (next_rc_bptr)
        {
            // check if there is enough separation to the next unmappable region or read end
            Size_Type separation_len = next_rc_bptr->get_r_len();
            while (separation_len <= _g_p->unmap_trigger_len())
            {
                next_rc_bptr = re_bptr->chunk_cont().get_sibling(next_rc_bptr, true, not dir).unconst();
                if (not next_rc_bptr or next_rc_bptr->ce_bptr()->is_unmappable())
                {
                    break;
                }
                separation_len += next_rc_bptr->get_r_len();
            }
            if (separation_len <= _g_p->unmap_trigger_len())
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
    if (_g_p->trim_tuc_step())
    {
        ASSERT(rc_bptr);
        ASSERT(rc_bptr->ce_bptr()->is_unmappable());
        if (rc_bptr->get_r_end() == re_bptr->end()
            or rc_bptr->get_r_start() == re_bptr->start())
        {
            _g_p->trim_tuc(rc_bptr);
        }
    }
    _g_p->check({re_bptr});
    LOG("Unmapper", debug) << ptree("end");
} // Unmapper::_extend_unmappable_re_region

void Unmapper::unmap_single_chunks() const
{
    LOG("Unmapper", info) << ptree("begin");
    ce_set_type unmap_ce_set;
    for (auto re_bptr : _g_p->re_cont() | referenced)
    {
        for (auto rc_bptr : re_bptr->chunk_cont() | referenced)
        {
            auto ce_bptr = rc_bptr->ce_bptr();
            if (not ce_bptr->is_unmappable() and size_one(ce_bptr->chunk_cont()))
            {
                LOG("Unmapper", debug) << ptree("single_chunk")
                    .put("re_bptr", re_bptr.to_int())
                    .put("rc_bptr", rc_bptr.to_int())
                    .put("ce_bptr", ce_bptr.to_int());
                unmap_ce_set.insert(ce_bptr);
            }
        }
    }
    _unmap_loop(move(unmap_ce_set), re_set_type());
    _g_p->check_all();
    LOG("Unmapper", info) << ptree("end");
} // Unmapper::unmap_single_chunks

void Unmapper::unmap_single_terminal_regions() const
{
    LOG("Unmapper", info) << ptree("begin");
    ce_set_type ce_set;
    re_set_type re_set;
    // Unmap terminal read chunk if it is mapped to a contig end, with no other supporting chunk.
    auto unmap_single_terminal_region = [&] (Read_Entry_BPtr re_bptr, bool r_dir) {
        if (re_bptr->chunk_cont().empty())
        {
            return;
        }
        Read_Chunk_BPtr rc_bptr = not r_dir? &*re_bptr->chunk_cont().rbegin() : &*re_bptr->chunk_cont().begin();
        Contig_Entry_BPtr ce_bptr = rc_bptr->ce_bptr();
        bool c_dir = (r_dir != rc_bptr->get_rc());
        if (ce_bptr->is_unmappable())
        {
            return;
        }
        if (size_one(ce_bptr->chunk_cont()))
        {
            ce_set.insert(ce_bptr);
            return;
        }
        if (c_dir)
        {
            if (// chunk does not span c_start
                rc_bptr->get_c_start() > 0
                // or it's not the first in the chunk container
                or &*ce_bptr->chunk_cont().begin() != rc_bptr
                // or a second chunk exists and starts at c_start
                or (not size_one(ce_bptr->chunk_cont()) and next(ce_bptr->chunk_cont().begin())->get_c_start() == 0))
            {
                return;
            }
            Size_Type c_brk = next(ce_bptr->chunk_cont().begin())->get_c_start();
            auto r_rg = rc_bptr->mapped_range(Range_Type(0, c_brk), true, true, true);        
            re_set[re_bptr].insert(r_rg);
        }
        else // not c_dir
        {
            if (// chunk does not span c_end
                rc_bptr->get_c_end() < ce_bptr->len())
            {
                return;
            }
            // this is trickier than the c_dir==true case
            // because chunks are not in the order of their end pos
            auto rg = ce_bptr->chunk_cont().iintersect(ce_bptr->len(), ce_bptr->len());
            // at least rc_bptr must span c_end
            ASSERT(rg.begin() != rg.end());
            if (not size_one(rg))
            {
                // more than 2 chunks span c_end
                return;
            }
            ASSERT(&*rg.begin() == rc_bptr);
            Size_Type c_brk = ce_bptr->chunk_cont().max_end(ce_bptr->len() - 1);
            // max_end smaller than contig end must exist because chunk_cont.size() >= 2 and rg.size() == 1
            ASSERT(c_brk < ce_bptr->len());
            // there can be no mutations in rc_bptr past c_brk
            ASSERT(rc_bptr->mut_ptr_cont().empty()
                   or rc_bptr->mut_ptr_cont().rbegin()->mut_cbptr()->rf_end() <= c_brk);
            auto r_rg = rc_bptr->mapped_range(Range_Type(c_brk, ce_bptr->len()), true, true, true);        
            re_set[re_bptr].insert(r_rg);
        }
    }; // unmap_single_terminal_region()
    for (auto re_bptr : _g_p->re_cont() | referenced)
    {
        unmap_single_terminal_region(re_bptr, false);
        unmap_single_terminal_region(re_bptr, true);
    }
    _unmap_loop(move(ce_set), move(re_set));
    _g_p->check_all();
    LOG("Unmapper", info) << ptree("end");
} // Unmapper::unmap_single_terminal_chunks

void Unmapper::unmap_short_contigs(unsigned min_len, unsigned max_deg) const
{
    LOG("Unmapper", info) << ptree("begin");
    ce_set_type ce_set;
    for (auto ce_bptr : _g_p->ce_cont() | referenced)
    {
        if (not ce_bptr->is_normal()) continue;
        if (ce_bptr->len() >= min_len) continue;
        auto oc_left = ce_bptr->out_chunks_dir(false, 3, 1);
        auto oc_right = ce_bptr->out_chunks_dir(true, 3, 1);
        if (oc_left.size() <= max_deg and oc_right.size() <= max_deg) continue;
        // unmap contig entry
        ce_set.insert(ce_bptr);
    }
    unmap(move(ce_set));
    _g_p->check_all();
    LOG("Unmapper", info) << ptree("end");
} // Unmapper::unmap_short_contigs

} // namespace MAC
