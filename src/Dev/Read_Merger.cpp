#include "Read_Merger.hpp"
#include "filter_cont.hpp"

namespace MAC
{

void Read_Merger::operator () ()
{
    LOG("graph", info) << ptree("begin");
    for (auto ce_bptr : _g.ce_cont() | referenced)
    {
        if (not ce_bptr->is_normal()) continue;
        for (auto mut_bptr : ce_bptr->mut_cont() | referenced)
        {
            Allele_Anchor anchor(mut_bptr);
            auto anchor_support = get_anchor_read_support(anchor, false);
            for (int al = 0; al < 2; ++al)
            {
                if (mut_bptr->get_copy_num(al) != 1) continue;
                // found a mutation allele with copy number 1
                // greb all reads supporting it
                Allele_Specifier allele(al);
                auto allele_support = get_allele_read_support(move(anchor_support), allele);
                extend_haploid_support(anchor, allele, false, allele_support);
                extend_haploid_support(anchor, allele, true, allele_support);
                merge_reads(allele_support);
            }
        }
    }
    LOG("graph", info) << ptree("end");
}

void Read_Merger::extend_haploid_support(
    const Allele_Anchor& _anchor, const Allele_Specifier& _allele, bool _c_direction,
    Allele_Read_Support& allele_support)
{
    auto crt_anchor = _anchor;
    auto crt_allele = _allele;
    auto crt_c_direction = _c_direction;
    while (true)
    {
        // compute next anchor
        Allele_Anchor next_anchor;
        Allele_Specifier next_allele;
        if (crt_anchor.is_endpoint() and crt_anchor.c_right() != crt_c_direction)
        {
            // this is the last anchor in the current contig direction
            // we jump to the associated anchor at the end of the contig edge
            if (not crt_allele.ce_next_cbptr())
            {
                // out-degree 0
                ASSERT(crt_allele.same_orientation());
                return;
            }
            bool next_c_direction = crt_c_direction == crt_allele.same_orientation();
            next_anchor = Allele_Anchor(crt_allele.ce_next_cbptr(), next_c_direction);
            next_allele = Allele_Specifier(crt_anchor.ce_cbptr(), crt_allele.same_orientation());
            swap(crt_c_direction, next_c_direction);
            swap(crt_anchor, next_anchor);
            swap(crt_allele, next_allele);
            continue;
        }
        next_anchor = crt_anchor.get_sibling(crt_c_direction);
        // compute support at next anchor
        auto next_anchor_support = get_anchor_read_support(next_anchor, crt_c_direction);
        // remove all reads not currently begin tracked
        for (auto& p : next_anchor_support)
        {
            for (int dir = 0; dir < 2; ++dir)
            {
                filter_cont(
                    p.second[dir],
                    [&] (Read_Entry_CBPtr re_cbptr) {
                        return allele_support[(dir + _c_direction) % 2].count(re_cbptr) > 0;
                    });
            }
        }
        filter_cont(
            next_anchor_support,
            [] (Anchor_Read_Support::const_iterator cit) {
                return not (cit->second[0].empty() and cit->second[1].empty());
            });
        if (next_anchor_support.empty())
        {
            break;
        }
        else if (next_anchor_support.size() > 1)
        {
            // check if there is a unique allele with support larger than the discordant threshold
            bool found_next_allele = false;
            for (const auto& p : next_anchor_support)
            {
                unsigned allele_support_size = p.second[0].size() + p.second[1].size();
                ASSERT(allele_support_size > 0);
                if (allele_support_size > _max_discordant_support)
                {
                    if (found_next_allele)
                    {
                        // this is a second allele with large support; we give up
                        LOG("Read_Merger", info) << ptree("discordant_support")
                            .put("anchor", _anchor)
                            .put("allele", _allele)
                            .put("crt_anchor", crt_anchor)
                            .put("crt_allele", crt_allele)
                            .put("next_anchor", next_anchor)
                            .put("next_allele_1", next_allele)
                            .put("next_allele_2", p.first);
                        allele_support[0].clear();
                        allele_support[1].clear();
                        return;
                    }
                    else
                    {
                        next_allele = p.first;
                        found_next_allele = true;
                    }
                }
            }
            if (not found_next_allele)
            {
                LOG("Read_Merger", info) << ptree("weak_ambiguous_support")
                    .put("anchor", _anchor)
                    .put("allele", _allele)
                    .put("crt_anchor", crt_anchor)
                    .put("crt_allele", crt_allele)
                    .put("next_anchor", next_anchor);
                allele_support[0].clear();
                allele_support[1].clear();
                return;
            }
            // support at next_anchor is ambiguous,
            // but a single allele has support greater than the discordant threshold
            for (auto& p : next_anchor_support)
            {
                unsigned allele_support_size = p.second[0].size() + p.second[1].size();
                ASSERT(allele_support_size > 0);
                if (allele_support_size > _max_discordant_support) continue;
                // we split the remaining reads in between anchor and next_anchor
                for (int dir = 0; dir < 2; ++dir)
                {
                    for (auto re_cbptr : p.second[dir])
                    {
                        ASSERT(allele_support[(dir + _c_direction) % 2].count(re_cbptr));
                        allele_support[(dir + _c_direction) % 2].erase(re_cbptr);
                        auto rp = split_read(
                            re_cbptr,
                            not crt_c_direction? crt_anchor : next_anchor,
                            not crt_c_direction? next_anchor : crt_anchor);
                        // we keep the side of the read that spans crt_anchor
                        ASSERT(get_allele_read_support(crt_anchor, crt_allele, crt_c_direction)[dir].count(rp.first) > 0);
                        allele_support[(dir + _c_direction) % 2].insert(rp.first);
                    }
                    p.second[dir].clear();
                }
            }
            filter_cont(
                next_anchor_support,
                [] (Anchor_Read_Support::const_iterator cit) {
                    return not (cit->second[0].empty() and cit->second[1].empty());
                });
        } // if (next_anchor_support.size() > 1)
        ASSERT(next_anchor_support.size() == 1);
        crt_anchor = next_anchor;
        crt_allele = next_anchor_support.begin()->first;
    }
} // Read_Merger::extend_haploid_support

void Read_Merger::merge_reads(const Allele_Read_Support& allele_support)
{
}

pair< Read_Entry_CBPtr, Read_Entry_CBPtr >
Read_Merger::split_read(Read_Entry_CBPtr re_cbptr, const Allele_Anchor& l_anchor, const Allele_Anchor& r_anchor)
{
    LOG("Read_Merger", debug) << ptree("begin")
        .put("re_bptr", re_cbptr.to_int())
        .put("re_name", re_cbptr->name())
        .put("l_anchor", l_anchor)
        .put("r_anchor", r_anchor);
    ASSERT(re_cbptr);
    auto re_bptr = re_cbptr.unconst();
    ASSERT(l_anchor.ce_cbptr() == r_anchor.ce_cbptr());
    ASSERT(l_anchor.get_sibling(false) == r_anchor);
    auto ce_bptr = l_anchor.ce_cbptr().unconst();
    auto l_out_pos = l_anchor.is_endpoint()? 0 : l_anchor.mut_cbptr()->rf_start();
    (void)l_out_pos;
    auto l_in_pos = l_anchor.is_endpoint()? 0 : l_anchor.mut_cbptr()->rf_end();
    auto r_in_pos = r_anchor.is_endpoint()? ce_bptr->len() : r_anchor.mut_cbptr()->rf_start();
    auto r_out_pos = r_anchor.is_endpoint()? ce_bptr->len() : r_anchor.mut_cbptr()->rf_end();
    (void)r_out_pos;
    ASSERT(l_in_pos < r_in_pos);
    //auto match_len = r_in_pos - l_in_pos;
    auto rc_bptr = ce_bptr->chunk_cont().search_read(re_bptr).unconst();
    auto dir = rc_bptr->get_rc();
    ASSERT(rc_bptr);
    ASSERT(rc_bptr->get_c_start() <= l_out_pos);
    ASSERT(r_out_pos <= rc_bptr->get_c_end());
    // create new Read_Entry object that will hold the tail
    _g.re_cont().erase(re_bptr);
    auto re_p = Read_Entry::split(rc_bptr, l_in_pos, r_in_pos);
    _g.re_cont().insert(re_p.first);
    _g.re_cont().insert(re_p.second);
    _g.check(set< Read_Entry_CBPtr >{ re_p.first, re_p.second });
    LOG("Read_Merger", debug) << ptree("end");
    return (not dir ? re_p : make_pair(re_p.second, re_p.first));
    /*
    // cut chunk at the beginning of the right anchor
    auto rc_p = Read_Chunk::split(rc_bptr, r_in_pos, nullptr, true);
    ASSERT(rc_p.first->get_c_end() == r_in_pos);
    ASSERT(rc_p.second->get_c_start() == r_in_pos);
    ASSERT(not l_anchor.is_mutation()
           or (not rc_p.first->mut_ptr_cont().empty()
               and &*rc_p.first->mut_ptr_cont().rbegin() == l_anchor.mut_cbptr()));
    ASSERT(not r_anchor.is_mutation()
           or (not rc_p.second->mut_ptr_cont().empty()
               and &*rc_p.second->mut_ptr_cont().begin() == r_anchor.mut_cbptr()));
    // extend second chunk to cover the match length
    ASSERT(rc_p.second->get_c_start() == l_in_pos);
    */
} // Read_Merger::split_read

} // namespace MAC
