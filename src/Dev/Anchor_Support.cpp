#include "Anchor_Support.hpp"
//#include "Hap_Hop.hpp"
#include "filter_cont.hpp"

namespace MAC
{

/*
Anchor_Chunk_Support get_anchor_chunk_support(const Allele_Anchor& anchor, unsigned min_edge_support)
{
    return anchor.chunk_support(min_edge_support);
}

Anchor_Read_Support get_anchor_read_support(const Anchor_Chunk_Support& anchor_chunk_support, bool c_direction)
{
    Anchor_Read_Support res;
    for (const auto& p : anchor_chunk_support)
    {
        for (int d = 0; d < 2; ++d)
        {
            for (auto rc_cbptr : p.second[d])
            {
                res[p.first][d != c_direction].insert(rc_cbptr->re_bptr());
            }
        }
    }
    return res;
}

Anchor_Read_Support get_anchor_read_support(const Allele_Anchor& anchor, unsigned min_edge_support, bool c_direction)
{
    return get_anchor_read_support(anchor.chunk_support(min_edge_support), c_direction);
}

Anchor_Read_Support get_hop_read_support(const Hap_Hop_CBPtr& hh_cbptr, bool h_direction)
{
    return get_anchor_read_support(hh_cbptr->allele_anchor(), hh_cbptr->c_direction() != h_direction);
}

Anchor_Read_Support get_hop_read_support(const pair< Hap_Hop_CBPtr, bool >& p)
{
    return get_hop_read_support(p.first, p.second);
}

Allele_Chunk_Support get_allele_chunk_support(const Anchor_Chunk_Support& anchor_chunk_support,
                                              const Allele_Specifier& allele)
{
    return anchor_chunk_support.at(allele);
}

Allele_Chunk_Support get_allele_chunk_support(Anchor_Chunk_Support&& anchor_chunk_support,
                                              const Allele_Specifier& allele)
{
    return move(anchor_chunk_support.at(allele));
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
                                            unsigned min_edge_support,
                                            const Allele_Specifier& allele,
                                            bool c_direction)
{
    auto anchor_rs = anchor.read_support(min_edge_support);
    if (not anchor_rs.count(allele)) return Allele_Read_Support();
    auto allele_rs = move(anchor_rs.at(allele));
    allele_rs.reverse(c_direction);
    return allele_rs;
    //return get_allele_read_support(get_anchor_read_support(anchor, min_edge_support, c_direction), allele);
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
*/

void subset_support(Anchor_Read_Support& anchor_read_support,
                    const Allele_Read_Support& common_support,
                    bool direction, bool remove_empty)
{
    for (auto& v : anchor_read_support | ba::map_values)
    {
        for (int dir = 0; dir < 2; ++dir)
        {
            filter_cont(
                v[dir],
                [&] (Read_Entry_CBPtr re_cbptr) {
                    return common_support[(dir + int(direction)) % 2].count(re_cbptr) > 0;
                });
        }
    }
    if (remove_empty)
    {
        filter_cont(
            anchor_read_support,
            [] (const Anchor_Read_Support::value_type& p) {
                return not (p.second[0].empty() and p.second[1].empty());
            });
    }
}

/*
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

RC_OSet get_oriented_chunks(Contig_Entry_CBPtr ce_cbptr, const RE_OSet& re_oset)
{
    RC_OSet res;
    for (int dir = 0; dir < 2; ++dir)
    {
        for (auto re_cbptr : re_oset[dir])
        {
            auto rc_cbptr = ce_cbptr->chunk_cont().search_read(re_cbptr);
            if (rc_cbptr)
            {
                res[dir].insert(rc_cbptr);
            }
        }
    }
    return res;
}
*/

} // namespace MAC
