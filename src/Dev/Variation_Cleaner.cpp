#include "Variation_Cleaner.hpp"


namespace MAC
{

void Variation_Cleaner::operator () () const
{
    unmap_clusters();
}

void Variation_Cleaner::unmap_clusters() const
{
    LOG("Variation_Cleaner", info) << ptree("begin");
    static const uint64_t visit_mask = 4;
    for (auto ce_bptr : _g.ce_cont() | referenced)
    {
        bitmask::reset(ce_bptr->tag(), visit_mask);
    }
    for (auto re_bptr : _g.re_cont() | referenced)
    {
        Size_Type r_pos = re_bptr->start();
        while (r_pos < re_bptr->end())
        {
            auto rc_cbptr = re_bptr->chunk_cont().get_chunk_with_pos(r_pos);
            ASSERT(rc_cbptr);
            Contig_Entry_BPtr ce_bptr = rc_cbptr->ce_bptr().unconst();
            if (not ce_bptr->is_normal() or bitmask::any(ce_bptr->tag(), visit_mask))
            {
                // skip this ce
                r_pos = rc_cbptr->get_r_end();
                continue;
            }

            auto p = find_mutation_cluster(ce_bptr, 0);
            if (p.first.empty() and not p.second[0])
            {
                r_pos = rc_cbptr->get_r_end();
                bitmask::set(ce_bptr->tag(), visit_mask);
                continue;
            }
            Range_Type rg(
                (p.second[0]
                 ? 0
                 : min_value_of(
                     p.first,
                     [] (Mutation_CBPtr mut_cbptr) { return mut_cbptr->rf_start(); })),
                (p.second[1]
                 ? ce_bptr->len()
                 : max_value_of(
                     p.first,
                     [] (Mutation_CBPtr mut_cbptr) { return mut_cbptr->rf_end(); })));
            LOG("Variation_Cleaner", debug) << ptree("unmap_ce_region")
                .put("ce_ptr", ce_bptr.to_int())
                .put("rg", rg);
            // save read entry positions to unmap in a vector
            //vector< tuple< Read_Entry_CBPtr, Range_Type > > unmap_re_pos_v;
            map< Read_Entry_BPtr, Range_Cont > unmap_re_set;
            for (auto unmap_rc_cbptr : ce_bptr->chunk_cont().iintersect(rg.begin(), rg.end()) | referenced)
            {
                // skip 0-length endpoint intersections
                // ONLY if cluster size is non-zero! otherwise, we might be trying to unmap an insertion
                if (rg.begin() < rg.end()
                    and (unmap_rc_cbptr->get_c_end() == rg.begin()
                         or unmap_rc_cbptr->get_c_start() == rg.end()))
                {
                    continue;
                }
                // transform cluster contig range into read range
                auto r_rg = unmap_rc_cbptr->mapped_range(rg, true, true, true);
                ASSERT(r_rg.begin() <= r_rg.end());
                if (r_rg.begin() == r_rg.end())
                {
                    // the unmapped range is spanned by a deletion in this read
                    continue;
                }
                LOG("Variation_Cleaner", debug) << ptree("unmap_re_region")
                    .put("re_ptr", unmap_rc_cbptr->re_bptr().to_int())
                    .put("rg", r_rg);
                //unmap_re_pos_v.emplace_back(unmap_rc_cbptr->re_bptr(), move(r_rg));
                unmap_re_set[unmap_rc_cbptr->re_bptr()].insert(r_rg);
            }
            _g.unmapper().unmap(move(unmap_re_set));
            r_pos = max(r_pos, re_bptr->start());
        }
    }
    _g.check_all();
    LOG("Variation_Cleaner", info) << ptree("end");
} // Variation_Cleaner::unmap_clusters

void Variation_Cleaner::sync_clusters() const
{
    LOG("Variation_Cleaner", info) << ptree("begin");
    for (auto ce_bptr : _g.ce_cont() | referenced) if (not ce_bptr->is_unmappable())
    {
        Size_Type start_pos = 0;
        while (true)
        {
            auto p = find_mutation_cluster(ce_bptr, start_pos);
            if (p.first.empty())
            {
                // cluster might have been left-right endpoint, but we don't care about those now
                break;
            }
            else if (p.first.size() == 1)
            {
                // cluster must involve at least one endpoint
                ASSERT(p.second[0] or p.second[1]);
                if (p.second[1])
                {
                    break;
                }
                else
                {
                    start_pos = p.first.begin()->rf_end();
                    continue;
                }
            }
            ASSERT(p.first.size() > 1);
            start_pos = max_value_of(
                p.first,
                [] (Mutation_CBPtr mut_cbptr) { return mut_cbptr->rf_end(); });
            sync_cluster(ce_cpbtr, p.first);
        } // while true
    } // for ce_bptr
    LOG("Variation_Cleaner", info) << ptree("end");
} // Variation_Cleaner::sync_clusters

pair< set< Mutation_CBPtr >, array< bool, 2 > >
Variation_Cleaner::find_mutation_cluster(Contig_Entry_CBPtr ce_cbptr, Size_Type start_pos) const
{
    pair< set< Mutation_CBPtr >, array< bool, 2 > > res({}, {{false, false}});
    auto& mut_set = res.first;
    auto& extend_to_end = res.second;
    Size_Type cluster_end = start_pos;

    auto well_separated = [&] (Size_Type pos1, Size_Type pos2) {
        return (pos1 + _min_num_3mers + 2 <= pos2
                and check_3mer_threshold(ce_cbptr->substr(pos1, pos2 - pos1)));
    };

    // find first mutation to consider
    Mutation_Cont::const_iterator mut_cit = ce_cbptr->mut_cont().begin();
    while (mut_cit != ce_cbptr->mut_cont().end()
           and mut_cit->rf_start() < start_pos)
    {
        ++mut_cit;
    }

    if (mut_cit == ce_cbptr->mut_cont().end())
    {
        // no mutations to consider
        extend_to_end[0] = not well_separated(start_pos, ce_cbptr->len());
        extend_to_end[1] = extend_to_end[0];
        return res;
    }

    extend_to_end[0] = not well_separated(start_pos, mut_cit->rf_start());
    for ( ; mut_cit != ce_cbptr->mut_cont().end(); ++mut_cit)
    {
        auto mut_cbptr = &*mut_cit;
        if (well_separated(cluster_end, mut_cbptr->rf_start()))
        {
            // mutation is well separated from any previous others
            if (mut_set.size() > 1 or (mut_set.size() == 1 and extend_to_end[0]))
            {
                // previous cluster was non-trivial; return it
                break;
            }
            else
            {
                // previous cluster was trivial; start a new one
                mut_set.clear();
            }
        }
        mut_set.insert(mut_cbptr);
        cluster_end = max(cluster_end, mut_cbptr->rf_end());
    }
    if (mut_cit == ce_cbptr->mut_cont().end())
    {
        // loop was stopped early; check separation with right end
        extend_to_end[1] = not well_separated(cluster_end, ce_cbptr->len());
    }
    if (mut_set.size() == 1 and not (extend_to_end[0] or extend_to_end[1]))
    {
        mut_set.clear();
    }
    return res;
} // Variation_Cleaner::find_mutation_cluster

bool Variation_Cleaner::check_3mer_threshold(const Seq_Proxy_Type& seq) const
{
    if (seq.size() < _min_num_3mers + 2)
    {
        return false;
    }
    set< Seq_Type > s;
    for (size_t i = 0; i + 2 < seq.size(); ++i)
    {
        Seq_Type kmer = seq.substr(i, 3);
        if (s.count(kmer) == 0)
        {
            if (s.size() + 1 >= _min_num_3mers)
            {
                return true;
            }
            s.insert(kmer);
        }
    }
    return false;
} // Variation_Cleaner::check_3mer_threshold

} // namespace MAC
