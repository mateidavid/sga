#include "Anchor_Support.hpp"
#include "Hap_Hop.hpp"

namespace MAC
{

Anchor_Chunk_Support get_anchor_chunk_support(const Allele_Anchor& anchor)
{
    Anchor_Chunk_Support res;
    if (anchor.is_mutation())
    {
        set< Read_Chunk_CBPtr > qr_set;
        set< Read_Chunk_CBPtr > full_rf_set;
        tie(qr_set, full_rf_set, ignore) = anchor.ce_cbptr()->mut_support(anchor.mut_cbptr());
        res[Allele_Specifier(false)] = move(full_rf_set);
        res[Allele_Specifier(true)] = move(qr_set);
    }
    else // is_endpoint
    {
        auto m = anchor.ce_cbptr()->out_chunks_dir(anchor.c_right(), 3, 1);
        if (not m.empty())
        {
            for (auto p : m)
            {
                res[Allele_Specifier(p.first)] = move(p.second);
            }
        }
        else
        {
            // out-degree 0
            auto pos = anchor.c_right()? anchor.ce_cbptr()->len() : 0;
            auto it_r = anchor.ce_cbptr()->chunk_cont().iintersect(pos, pos);
            auto it_rr = make_ref_range(it_r);
            res[Allele_Specifier(nullptr, true)].insert(it_rr.begin(), it_rr.end());
            ASSERT(not res.at(Allele_Specifier(nullptr, true)).empty());
        }
    }
    return res;
}

Anchor_Read_Support get_anchor_read_support(const Anchor_Chunk_Support& anchor_chunk_support, bool c_direction)
{
    Anchor_Read_Support res;
    for (const auto& p : anchor_chunk_support)
    {
        for (const auto rc_cbptr : p.second)
        {
            res[p.first][rc_cbptr->get_rc() != c_direction].insert(rc_cbptr->re_bptr());
        }
    }
    return res;
}

Anchor_Read_Support get_anchor_read_support(const Allele_Anchor& anchor, bool c_direction)
{
    return get_anchor_read_support(get_anchor_chunk_support(anchor), c_direction);
}

Anchor_Read_Support get_hop_read_support(const Hap_Hop_CBPtr& hh_cbptr, bool h_direction)
{
    return get_anchor_read_support(hh_cbptr->allele_anchor(), hh_cbptr->c_direction() != h_direction);
}

Anchor_Read_Support get_hop_read_support(const pair< Hap_Hop_CBPtr, bool >& p)
{
    return get_hop_read_support(p.first, p.second);
}

Allele_Read_Support get_allele_read_support(const Anchor_Read_Support& anchor_read_support,
                                            const Allele_Specifier& allele)
{
    return anchor_read_support.at(allele);
}

Allele_Read_Support get_allele_read_support(Anchor_Read_Support&& anchor_read_support,
                                            const Allele_Specifier& allele)
{
    return move(anchor_read_support[allele]);
}

Allele_Read_Support get_allele_read_support(const Allele_Anchor& anchor,
                                            const Allele_Specifier& allele,
                                            bool c_direction)
{
    return get_allele_read_support(get_anchor_read_support(anchor, c_direction), allele);
}

Allele_Read_Support collapsed_support(const Anchor_Read_Support& anchor_read_support)
{
    Allele_Read_Support res;
    for (const auto& v : anchor_read_support | ba::map_values)
    {
        for (int dir = 0; dir < 2; ++dir)
        {
            res[dir].insert(v[dir].begin(), v[dir].end());
        }
    }
    return res;
}

void subset_support(Anchor_Read_Support& anchor_read_support, const Allele_Read_Support& common_support, bool direction)
{
    for (auto& v : anchor_read_support | ba::map_values)
    {
        for (int dir = 0; dir < 2; ++dir)
        {
            for_each_advance(
                v[dir].begin(),
                v[dir].end(),
                [&] (Read_Entry_CBPtr re_cbptr) {
                    if (common_support[(dir + int(direction)) % 2].count(re_cbptr) == 0)
                    {
                        v[dir].erase(re_cbptr);
                    }
                });
        }
    }
}

Allele_Read_Support common_allele_support(const Allele_Read_Support& al1_read_support,
                                          const Allele_Read_Support& al2_read_support,
                                          bool direction)
{
    Allele_Read_Support res;
    for (int dir = 0; dir < 2; ++dir)
    {
        for (auto re_cbptr : al1_read_support[dir])
        {
            if (al2_read_support[(dir + int(direction)) % 2].count(re_cbptr))
            {
                res[dir].insert(re_cbptr);
            }
        }
    }
    return res;
}

} // namespace MAC
