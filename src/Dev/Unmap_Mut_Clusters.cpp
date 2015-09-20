#include "Unmap_Mut_Clusters.hpp"


namespace MAC
{

void Unmap_Mut_Clusters::operator () (Graph& g) const
{
    LOG("Unmap_Mut_Clusters", info) << ptree("begin");
    static const uint64_t visit_mask = 4;
    for (auto ce_bptr : g.ce_cont() | referenced)
    {
        bitmask::reset(ce_bptr->tag(), visit_mask);
    }
    for (auto re_bptr : g.re_cont() | referenced)
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
            auto mut_cluster = find_mutation_cluster(ce_bptr);
            if (not mut_cluster)
            {
                r_pos = rc_cbptr->get_r_end();
                bitmask::set(ce_bptr->tag(), visit_mask);
                continue;
            }
            LOG("Unmap_Mut_Clusters", debug) << ptree("unmap_ce_region")
                .put("ce_ptr", ce_bptr.to_int())
                .put("rg", *mut_cluster);
            // save read entry positions to unmap in a vector
            //vector< tuple< Read_Entry_CBPtr, Range_Type > > unmap_re_pos_v;
            map< Read_Entry_BPtr, Range_Cont > unmap_re_set;
            for (auto unmap_rc_cbptr : ce_bptr->chunk_cont().iintersect(mut_cluster->begin(), mut_cluster->end()) | referenced)
            {
                // skip 0-length endpoint intersections
                // ONLY if cluster size is non-zero! otherwise, we might be trying to unmap an insertion
                if (mut_cluster->begin() < mut_cluster->end()
                    and (unmap_rc_cbptr->get_c_end() == mut_cluster->begin()
                         or unmap_rc_cbptr->get_c_start() == mut_cluster->end()))
                {
                    continue;
                }
                // transform cluster contig range into read range
                auto r_rg = unmap_rc_cbptr->mapped_range(mut_cluster.get(), true, true, true);
                ASSERT(r_rg.begin() <= r_rg.end());
                if (r_rg.begin() == r_rg.end())
                {
                    // the unmapped range is spanned by a deletion in this read
                    continue;
                }
                LOG("Unmap_Mut_Clusters", debug) << ptree("unmap_re_region")
                    .put("re_ptr", unmap_rc_cbptr->re_bptr().to_int())
                    .put("rg", r_rg);
                //unmap_re_pos_v.emplace_back(unmap_rc_cbptr->re_bptr(), move(r_rg));
                unmap_re_set[unmap_rc_cbptr->re_bptr()].insert(r_rg);
            }
            /*
            for (const auto& re_pos : unmap_re_pos_v)
            {
                g.unmap_re_region(get<0>(re_pos).unconst(), get<1>(re_pos));
            }
            */
            g.unmap_re_regions(move(unmap_re_set));
        }
    }
    LOG("Unmap_Mut_Clusters", info) << ptree("end");
}

optional< Range_Type >
Unmap_Mut_Clusters::find_mutation_cluster(Contig_Entry_CBPtr ce_cbptr) const
{
    optional< Range_Type > res;
    Mutation_CBPtr start_mut_cbptr = Mutation_Fact::new_elem(0, 0);
    Mutation_CBPtr end_mut_cbptr = Mutation_Fact::new_elem(ce_cbptr->len(), 0);
    Mutation_CBPtr last_mut_cbptr = start_mut_cbptr;
    vector< Mutation_CBPtr > tmp_v(1, end_mut_cbptr);
    for (auto mut_cbptr : join(
             ce_cbptr->mut_cont() | referenced,
             tmp_v
             )
        )
    {
        if (mut_cbptr->rf_start() < last_mut_cbptr->rf_end() + _min_num_3mers + 2
            or not check_3mer_threshold(ce_cbptr->substr(last_mut_cbptr->rf_end(), mut_cbptr->rf_start() - last_mut_cbptr->rf_end())))
        {
            // current mutation not well separated from last one
            if (res)
            {
                // extend region
                res->end() = max(res->end(), mut_cbptr->rf_end());
            }
            else
            {
                // initialize region and continue
                res = move(Range_Type(last_mut_cbptr->rf_start(), mut_cbptr->rf_end()));
            }
        }
        else if (res)
        {
            // region exists but mut_cbptr is outside of it
            break;
        }
        last_mut_cbptr = mut_cbptr;
    }
    Mutation_Fact::del_elem(start_mut_cbptr);
    Mutation_Fact::del_elem(end_mut_cbptr);
    return res;
}

bool Unmap_Mut_Clusters::check_3mer_threshold(const Seq_Proxy_Type& seq) const
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
}


} // namespace MAC
