#include "Hap_Map.hpp"


namespace MAC
{

Hap_Map::Hap_Map(const Graph& g)
{
    build(g);
}

Hap_Hop_BPtr Hap_Map::make_single_allele_hop(const Allele_Anchor& anchor, const Allele_Specifier& allele,
                                             bool direction)
{
    Hap_Entry_BPtr he_bptr = Hap_Entry_Fact::new_elem();
    _he_cont.insert(he_bptr);
    Hap_Hop_BPtr hh_bptr = Hap_Hop_Fact::new_elem(he_bptr, anchor, allele, direction);
    he_bptr->hh_cont().push_front(hh_bptr);
    _hh_set.insert(hh_bptr);
    check_he(he_bptr);
    return hh_bptr;
} // Hap_Map::make_single_allele_hop

map< Allele_Specifier, set< Hap_Hop_CBPtr > > Hap_Map::make_mutation_haps(Mutation_CBPtr mut_cbptr)
{
    map< Allele_Specifier, set< Hap_Hop_CBPtr > > res;
    Allele_Anchor anchor(mut_cbptr);
    for (int k = 0; k < 2; ++k)
    {
        Allele_Specifier allele(k == 1);
        Hap_Hop_BPtr hh_bptr = make_single_allele_hop(anchor, allele);
        res[allele].insert(hh_bptr);
    }
    return res;
} // Hap_Map::make_mutation_haps

map< Allele_Specifier, set< Hap_Hop_CBPtr > > Hap_Map::extend_endpoint_haps(const Allele_Anchor& anchor)
{
    ASSERT(anchor.is_endpoint());
    map< Allele_Specifier, set< Hap_Hop_CBPtr > > res;
    auto anchor_support = anchor.support();
    for (auto p : anchor_support)
    {
        Allele_Specifier allele = p.first;
        ASSERT(allele.is_endpoint());
        // compute corresponding anchor at the other endpoint of the edge
        Allele_Anchor mirror_anchor(allele.ce_next_cbptr(), anchor.c_right() != allele.same_orientation());
        // compute corresponding allele at the other endpoint
        Allele_Specifier mirror_allele(anchor.ce_cbptr(), allele.same_orientation());
        // check if mirror_anchor has haplotypes
        set< Hap_Hop_BPtr > haps_to_extend;
        auto it_p = _hh_set.equal_range(mirror_anchor, detail::Hap_Hop_Comparator());
        for (auto it = it_p.first; it != it_p.second; ++it)
        {
            auto hh_bptr = &*it;
            if (hh_bptr->allele_specifier() != mirror_allele) continue;
            haps_to_extend.insert(hh_bptr);
        }
        // check: mirror anchor was visited iff hh_set contains a hop at the mirror allele
        ASSERT((it_p.first == it_p.second) == haps_to_extend.empty());
        if (haps_to_extend.empty())
        {
            // mirror allele not visited yet
            // construct new haplotype containing single hop
            Hap_Hop_BPtr hh_bptr = make_single_allele_hop(anchor, allele);
            res[allele].insert(hh_bptr);
        }
        else
        {
            // haplotypes exist at mirror allele
            // we extend them to this allele
            for (auto hh_bptr : haps_to_extend)
            {
                ASSERT(hh_bptr);
                Hap_Entry_BPtr he_bptr = hh_bptr->he_cbptr().unconst();
                ASSERT(he_bptr);
                ASSERT(not he_bptr->hh_cont().empty());
                ASSERT(hh_bptr == he_bptr->hh_cont().front()
                       or hh_bptr == he_bptr->hh_cont().back());
                bool is_hap_end = hh_bptr == he_bptr->hh_cont().back();
                // create new hop for current allele
                Hap_Hop_BPtr hh_new_bptr = Hap_Hop_Fact::new_elem(he_bptr, anchor, allele,
                                                                  is_hap_end == anchor.c_right());
                if (is_hap_end)
                {
                    he_bptr->hh_cont().insert_after(hh_bptr, hh_new_bptr);
                }
                else
                {
                    he_bptr->hh_cont().insert_before(hh_bptr, hh_new_bptr);
                }
                _hh_set.insert(hh_new_bptr);
                res[allele].insert(hh_new_bptr);
                check_he(he_bptr);
            }
        }
        ASSERT(res.count(allele) == 1);
        ASSERT(not res[allele].empty());
    } // for (p : anchor_support)
    ASSERT(res.size() == anchor_support.size());
    return res;
} // Hap_Map::extend_haps_at_endpoint

void Hap_Map::build(const Graph& g)
{
    for (auto ce_cbptr : g.ce_cont() | referenced)
    {
        if (not ce_cbptr->is_normal()) continue;

        // construct left endpoint haps
        Allele_Anchor prev_anchor(ce_cbptr, false);
        auto prev_haps = extend_endpoint_haps(prev_anchor);
        // construct disconnected haps at every mutation
        for (auto mut_cbptr : ce_cbptr->mut_cont() | referenced)
        {
            Allele_Anchor anchor(mut_cbptr);
            auto haps = make_mutation_haps(mut_cbptr);
            // connect with haps at previous anchor
            connect_unique(prev_anchor, anchor, prev_haps, haps);
            prev_anchor = move(anchor);
            prev_haps = move(haps);
        }
        // construct right endpoint haps
        Allele_Anchor re_anchor(ce_cbptr, true);
        auto re_haps = extend_endpoint_haps(re_anchor);
        // connect with haps at previous anchor
        connect_unique(prev_anchor, re_anchor, prev_haps, re_haps);
    }
} // Hap_Map::build

void Hap_Map::connect_unique(const Allele_Anchor& a1, const Allele_Anchor& a2,
                             map< Allele_Specifier, set< Hap_Hop_CBPtr > >& a1_hops,
                             map< Allele_Specifier, set< Hap_Hop_CBPtr > >& a2_hops)
{
    // check there is 1 hop per allele, and all hops are terminal
    auto check_hop_terminal = [] (Hap_Hop_CBPtr hh_cbptr) {
        return (hh_cbptr == hh_cbptr->he_cbptr()->hh_cont().front()
                or hh_cbptr == hh_cbptr->he_cbptr()->hh_cont().back());
    };
    static_cast< void >(check_hop_terminal);
    ASSERT(all_of(a1_hops.begin(), a1_hops.end(),
                  [&check_hop_terminal] (const pair< Allele_Specifier, set< Hap_Hop_CBPtr > >& p) {
                      return p.second.size() == 1 and check_hop_terminal(*p.second.begin());
                  }));
    ASSERT(all_of(a2_hops.begin(), a2_hops.end(),
                  [&check_hop_terminal] (const pair< Allele_Specifier, set< Hap_Hop_CBPtr > >& p) {
                      return p.second.size() == 1 and check_hop_terminal(*p.second.begin());
                  }));
    auto a1_support = a1.support();
    ASSERT(a1_support.size() == a1_hops.size());
    auto a2_support = a2.support();
    ASSERT(a2_support.size() == a2_hops.size());
    auto connect_map = Allele_Anchor::connect(a1_support, a2_support);
    // create map of unique allele connections
    map< Allele_Specifier, Allele_Specifier > unique_conn_map[2];
    auto make_unique_conn_map = [&] (bool a1_to_a2) {
        for (auto a1_allele : (a1_to_a2? a1_support : a2_support) | ba::map_keys)
        {
            Allele_Specifier unique_a2_allele;
            int d = 0;
            for (auto a2_allele : (a1_to_a2? a2_support : a1_support) | ba::map_keys)
            {
                if (connect_map.count(a1_to_a2? make_pair(a1_allele, a2_allele) : make_pair(a2_allele, a1_allele)) > 0)
                {
                    ++d;
                    if (d > 1) break;
                    unique_a2_allele = a2_allele;
                }
            }
            if (d != 1) continue;
            unique_conn_map[static_cast< int >(not a1_to_a2)][a1_allele] = unique_a2_allele;
        }
    };
    make_unique_conn_map(true);
    make_unique_conn_map(false);
    for (auto p : unique_conn_map[0])
    {
        auto a1_allele = p.first;
        auto a2_allele = p.second;
        if (unique_conn_map[1].count(a2_allele) != 1) continue;
        ASSERT(unique_conn_map[1][a2_allele] == a1_allele);
        // unique connection: a1_allele to a2_allele
        auto a1_hop_cbptr = *a1_hops[a1_allele].begin();
        auto a2_hop_cbptr = *a2_hops[a2_allele].begin();
        auto he1_bptr = a1_hop_cbptr->he_cbptr().unconst();
        auto he2_bptr = a2_hop_cbptr->he_cbptr().unconst();
        // if they are on the same haplotype, this consitututes a haplotype cycle
        if (he1_bptr == he2_bptr)
        {
            logger("Hap_Map", debug) << ptree("hap_cycle")
                .put("a1_hop", a1_hop_cbptr->to_ptree()).put("a2_hop", a2_hop_cbptr->to_ptree());
            continue;
        }
        if (a1_hop_cbptr->c_direction() != a2_hop_cbptr->c_direction())
        {
            // reverse haplotype at a2_hop to enable merge
            he2_bptr->hh_cont().reverse();
        }
        ASSERT(a1_hop_cbptr->c_direction() == a2_hop_cbptr->c_direction());
        // merge haplotypes
        if (a1_hop_cbptr->c_direction() == false)
        {
            he1_bptr->hh_cont().splice_right(he2_bptr->hh_cont(), he1_bptr);
            _he_cont.erase(he2_bptr);
            Hap_Entry_Fact::del_elem(he2_bptr);
        }
        else
        {
            he2_bptr->hh_cont().splice_right(he1_bptr->hh_cont(), he2_bptr);
            _he_cont.erase(he1_bptr);
            Hap_Entry_Fact::del_elem(he1_bptr);
        }
    }
} // Hap_Map::connect_unique

void Hap_Map::check_he(Hap_Entry_CBPtr he_cbptr) const
{
    static_cast< void >(he_cbptr);
#ifndef BOOST_DISABLE_ASSERTS
    ASSERT(he_cbptr);
    Hap_Hop_CBPtr last_hh_cbptr;
    for (auto hh_cbptr : he_cbptr->hh_cont() | referenced)
    {
        ASSERT(hh_cbptr);
        ASSERT(hh_cbptr->he_cbptr() == he_cbptr);
        if (last_hh_cbptr)
        {
            if (last_hh_cbptr->ce_cbptr() == hh_cbptr->ce_cbptr())
            {
                ASSERT(last_hh_cbptr->c_direction() == hh_cbptr->c_direction());
                // anchors must be in an order corresponding to the direction
                ASSERT(not hh_cbptr->c_direction()?
                       last_hh_cbptr->allele_anchor() < hh_cbptr->allele_anchor()
                       : hh_cbptr->allele_anchor() < last_hh_cbptr->allele_anchor());
            }
            else
            {
                // these must be mirror hops
                ASSERT(last_hh_cbptr->allele_anchor().is_endpoint());
                ASSERT(hh_cbptr->allele_anchor().is_endpoint());
                bool same_orientation = last_hh_cbptr->c_direction() == hh_cbptr->c_direction();
                ASSERT(last_hh_cbptr->allele_specifier() == Allele_Specifier(hh_cbptr->ce_cbptr(), same_orientation));
                ASSERT(hh_cbptr->allele_specifier() == Allele_Specifier(last_hh_cbptr->ce_cbptr(), same_orientation));
                auto last_hh_allele_support = last_hh_cbptr->allele_anchor().support();
                auto hh_allele_support = hh_cbptr->allele_anchor().support();
                ASSERT(last_hh_allele_support.count(last_hh_cbptr->allele_specifier()) > 1);
                ASSERT(last_hh_allele_support.count(last_hh_cbptr->allele_specifier())
                       == hh_allele_support.count(hh_cbptr->allele_specifier()));
            }
        }
        last_hh_cbptr = hh_cbptr;
    }
#endif
} // Hap_Map::check_he

void Hap_Map::check(const set< Hap_Entry_CBPtr >& he_set) const
{
    static_cast< void >(he_set);
#ifndef BOOST_DISABLE_ASSERTS
    for (auto he_cbptr : he_set)
    {
        check_he(he_cbptr);
    }
#endif
} // Hap_Map::check

void Hap_Map::check_all() const
{
#ifndef BOOST_DISABLE_ASSERTS
    for (auto he_cbptr : _he_cont | referenced)
    {
        check_he(he_cbptr);
    }
#endif
} // Hap_Map::check_all


} // namespace MAC
