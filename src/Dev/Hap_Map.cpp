#include "Hap_Map.hpp"


namespace MAC
{

Hap_Hop_BPtr Hap_Map::make_single_hop_hap(const Allele_Anchor& anchor, const Allele_Specifier& allele,
                                          bool c_direction)
{
    Hap_Entry_BPtr he_bptr = Hap_Entry_Fact::new_elem();
    he_cont().insert(he_bptr);
    Hap_Hop_BPtr hh_bptr = Hap_Hop_Fact::new_elem(he_bptr, anchor, allele, c_direction);
    he_bptr->hh_cont().push_front(hh_bptr);
    hh_set().insert(hh_bptr);
    check_he(he_bptr);
    return hh_bptr;
} // Hap_Map::make_single_hop_hap

map< Allele_Specifier, set< Hap_Hop_CBPtr > >
Hap_Map::make_mutation_haps(Mutation_CBPtr mut_cbptr)
{
    LOG("hap_map", debug) << ptree("make_mutation_haps").put("mut_ptr", mut_cbptr.to_int());
    map< Allele_Specifier, set< Hap_Hop_CBPtr > > res;
    Allele_Anchor anchor(mut_cbptr);
    auto anchor_support = anchor.support();
    for (auto p : anchor_support)
    {
        Allele_Specifier allele = p.first;
        Hap_Hop_BPtr hh_bptr = make_single_hop_hap(anchor, allele, false);
        res[allele].insert(hh_bptr);
    }
    return res;
} // Hap_Map::make_mutation_haps

map< Allele_Specifier, set< Hap_Hop_CBPtr > >
Hap_Map::make_endpoint_haps(Contig_Entry_CBPtr ce_cbptr, bool c_right)
{
    LOG("hap_map", debug) << ptree("make_endpoint_haps")
        .put("ce_ptr", ce_cbptr.to_int())
        .put("c_right", c_right);
    map< Allele_Specifier, set< Hap_Hop_CBPtr > > res;
    Allele_Anchor anchor(ce_cbptr, c_right);
    auto anchor_support = anchor.support();
    for (auto p : anchor_support)
    {
        Allele_Specifier allele = p.first;
        Hap_Hop_BPtr hh_bptr = make_single_hop_hap(anchor, allele, false);
        res[allele].insert(hh_bptr);
    }
    return res;
} // Hap_Map::make_endpoint_haps

void Hap_Map::initialize_haps(const Graph& g)
{
    ASSERT(he_cont().empty());
    ASSERT(hh_set().empty());
    for (auto ce_cbptr : g.ce_cont() | referenced)
    {
        if (not ce_cbptr->is_normal()) continue;
        make_endpoint_haps(ce_cbptr, false);
        make_endpoint_haps(ce_cbptr, true);
        for (auto mut_cbptr : ce_cbptr->mut_cont() | referenced)
        {
            make_mutation_haps(mut_cbptr);
        }
    }
} // Hap_Map::initialize_haps

void Hap_Map::connect_endpoint_haps(const Graph& g)
{
    for (auto ce_cbptr : g.ce_cont() | referenced)
    {
        if (not ce_cbptr->is_normal()) continue;
        for (int c_right = 0; c_right < 2; ++c_right)
        {
            Allele_Anchor anchor(ce_cbptr, c_right);
            auto p = hh_set().equal_range(anchor);
            for (auto it = p.first; it != p.second; ++it)
            {
                auto hh_bptr = &*it;
                ASSERT(hh_bptr->allele_anchor() == anchor);
                if (not size_one(hh_bptr->he_cbptr()->hh_cont()))
                {
                    continue;
                }
                ASSERT(size_one(hh_bptr->he_cbptr()->hh_cont()));
                // fetch mirror hap
                auto ce_next_cbptr = hh_bptr->allele_specifier().ce_next_cbptr();
                bool same_orientation = hh_bptr->allele_specifier().same_orientation();
                Allele_Anchor mirror_anchor(ce_next_cbptr, c_right != same_orientation);
                Allele_Specifier mirror_allele(ce_cbptr, same_orientation);
                auto mirror_hh_m = hh_set().find_allele(mirror_anchor, mirror_allele, false);
                ASSERT(mirror_hh_m.size() == 1);
                auto mirror_hh_bptr = mirror_hh_m.begin()->first.unconst();
                ASSERT(size_one(mirror_hh_bptr->he_cbptr()->hh_cont()));
                //
                // connect hh_bptr and mirror_hh_bptr into single he
                //
                merge_haps(hh_bptr, mirror_hh_bptr,
                           not c_right, not c_right == same_orientation);
            }
        }
    }
} // Hap_Map::connect_endpoint_haps

void Hap_Map::merge_haps(Hap_Hop_BPtr hh1_bptr, Hap_Hop_BPtr hh2_bptr, bool c1_dir, bool c2_dir)
{
    auto he1_bptr = hh1_bptr->he_cbptr().unconst();
    auto he2_bptr = hh2_bptr->he_cbptr().unconst();
    ASSERT(he1_bptr != he2_bptr);
    if (hh1_bptr->c_direction() != c1_dir)
    {
        he1_bptr->hh_cont().reverse();
    }
    if (hh2_bptr->c_direction() != c2_dir)
    {
        he2_bptr->hh_cont().reverse();
    }
    he1_bptr->hh_cont().splice_right(he2_bptr->hh_cont(), he1_bptr);
    he_cont().erase(he2_bptr);
    Hap_Entry_Fact::del_elem(he2_bptr);
} // Hap_Map::merge_haps

void Hap_Map::connect_consecutive_anchors(const Allele_Anchor& a1)
{
    ASSERT(not a1.is_endpoint() or not a1.c_right());
    Allele_Anchor a2 = a1.get_sibling(false);
    auto a1_hops = hh_set().find_anchor(a1, true);
    auto a2_hops = hh_set().find_anchor(a2, false);
    // check that all haps end at a1&a2
    auto check_hops_terminal = [] (const Hap_Hop_Set::find_anchor_type& hops) {
        return all_of(
            hops.begin(),
            hops.end(),
            [] (const Hap_Hop_Set::find_anchor_type::value_type& v) {
                return all_of(
                    v.second.begin(),
                    v.second.end(),
                    [] (const pair< Hap_Hop_CBPtr, bool >& p) {
                        return Hap_Map::is_end_hop(p.first, p.second);
                    });
            });
    };
    ASSERT(check_hops_terminal(a1_hops));
    ASSERT(check_hops_terminal(a2_hops));
    //Hap_Tree a1_tree(a1, a1_hops);
    //Hap_Tree a2_tree(a2, a2_hops);
}

/*
map< Allele_Specifier, set< Hap_Hop_CBPtr > > Hap_Map::extend_endpoint_haps(const Allele_Anchor& anchor)
{
    ASSERT(anchor.is_endpoint());
    LOG("hap_map", debug) << ptree("extend_endpoint_haps")
        .put("ce_ptr", anchor.ce_cbptr().to_int())
        .put("c_right", anchor.c_right());
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
        auto it_p = hh_set().equal_range(mirror_anchor);
        for (auto it = it_p.first; it != it_p.second; ++it)
        {
            auto hh_bptr = &*it;
            ASSERT(hh_bptr->allele_anchor() == mirror_anchor);
            if (hh_bptr->allele_specifier() != mirror_allele) continue;
            haps_to_extend.insert(hh_bptr);
        }
        // check: mirror anchor was visited iff hh_set contains a hop at the mirror allele
        ASSERT((it_p.first == it_p.second) == haps_to_extend.empty());
        if (haps_to_extend.empty())
        {
            // mirror allele not visited yet
            // construct new haplotype containing single hop
            Hap_Hop_BPtr hh_bptr = make_single_hop_hap(anchor, allele, false);
            res[allele].insert(hh_bptr);
        }
        else
        {
            // haplotypes exist at mirror allele
            // we extend them to this allele
            // during this stage, there can be at most one haplotype to extend
            ASSERT(haps_to_extend.size() == 1);
            for (auto hh_bptr : haps_to_extend)
            {
                ASSERT(hh_bptr);
                Hap_Entry_BPtr he_bptr = hh_bptr->he_cbptr().unconst();
                ASSERT(he_bptr);
                ASSERT(not he_bptr->hh_cont().empty());
                // create new hop for current allele
                Hap_Hop_BPtr hh_new_bptr = Hap_Hop_Fact::new_elem(he_bptr, anchor, allele,
                                                                  hh_bptr->c_direction() == allele.same_orientation());
                // direction and endpoint are enough to identify if this should be first or last hop
                if (hh_new_bptr->c_direction() == anchor.c_right())
                {
                    ASSERT(Hap_Map::is_end_hop(hh_bptr, true));
                    he_bptr->hh_cont().insert_after(hh_bptr, hh_new_bptr);
                }
                else
                {
                    ASSERT(Hap_Map::is_end_hop(hh_bptr, false));
                    he_bptr->hh_cont().insert_before(hh_bptr, hh_new_bptr);
                }
                hh_set().insert(hh_new_bptr);
                res[allele].insert(hh_new_bptr);
                check_he(he_bptr);
            }
        }
        ASSERT(res.count(allele) == 1);
        ASSERT(not res[allele].empty());
    } // for (p : anchor_support)
    ASSERT(res.size() == anchor_support.size());
    return res;
} // Hap_Map::extend_endpoint_haps
*/

void Hap_Map::build(const Graph& g)
{
    LOG("hap_map", info) << ptree("build_start");
    initialize_haps(g);
    check_all();
    connect_endpoint_haps(g);
    check_all();
    for (auto ce_cbptr : g.ce_cont() | referenced)
    {
        if (not ce_cbptr->is_normal()) continue;
        for (Allele_Anchor anchor = Allele_Anchor(ce_cbptr, false);
             anchor != Allele_Anchor(ce_cbptr, right);
             anchor = anchor.get_sibling(false))
        {
            connect_consecutive_anchors(anchor);
        }
    }
/*
    clear_and_dispose();
    for (auto ce_cbptr : g.ce_cont() | referenced)
    {
        if (not ce_cbptr->is_normal()) continue;
        LOG("hap_map", debug) << ptree("build_loop").put("ce_ptr", ce_cbptr.to_int());
        // do not allow overlapping mutations during haplotype building
        ASSERT(ce_cbptr->separated_mutations(10, true));
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
        hh_set().check();
    }
*/
    LOG("hap_map", info) << ptree("build_end");
} // Hap_Map::build

/*
void Hap_Map::connect_unique(const Allele_Anchor& a1, const Allele_Anchor& a2,
                             map< Allele_Specifier, set< Hap_Hop_CBPtr > >& a1_hops,
                             map< Allele_Specifier, set< Hap_Hop_CBPtr > >& a2_hops)
{
    ASSERT(a1.get_sibling(true) == a2);
    // check there is 1 hop per allele, and all hops are terminal
    ASSERT(all_of(a1_hops.begin(), a1_hops.end(),
                  [] (const pair< Allele_Specifier, set< Hap_Hop_CBPtr > >& p) {
                      return p.second.size() == 1 and Hap_Map::is_end_hop(*p.second.begin(), not (*p.second.begin())->c_direction());
                  }));
    ASSERT(all_of(a2_hops.begin(), a2_hops.end(),
                  [] (const pair< Allele_Specifier, set< Hap_Hop_CBPtr > >& p) {
                      return p.second.size() == 1 and Hap_Map::is_end_hop(*p.second.begin(), (*p.second.begin())->c_direction());
                  }));
    auto a1_support = a1.support();
    auto a2_support = a2.support();
    ASSERT(a1_support.size() == a1_hops.size());
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
        // if they are on the same haplotype, this constitutes a haplotype cycle
        if (he1_bptr == he2_bptr)
        {
            LOG("hap_map", debug) << ptree("connect_unique__hap_cycle")
                .put("a1_hop", a1_hop_cbptr->to_ptree()).put("a2_hop", a2_hop_cbptr->to_ptree());
            continue;
        }
        LOG("hap_map", debug) << ptree("connect_unique__new_connection")
            .put("a1_anchor", a1.to_ptree())
            .put("a2_anchor", a2.to_ptree())
            .put("a1_allele", a1_allele.to_ptree())
            .put("a2_allele", a2_allele.to_ptree());
        if (a1_hop_cbptr->c_direction() != a2_hop_cbptr->c_direction())
        {
            // reverse haplotype at a2_hop to enable merge
            LOG("hap_map", debug) << ptree("connect_unique__reverse_hap")
                .put("he_ptr", he2_bptr.to_int());
            he2_bptr->hh_cont().reverse();
        }
        ASSERT(a1_hop_cbptr->c_direction() == a2_hop_cbptr->c_direction());
        // merge haplotypes
        if (a1_hop_cbptr->c_direction() == false)
        {
            he1_bptr->hh_cont().splice_right(he2_bptr->hh_cont(), he1_bptr);
            _he_cont.erase(he2_bptr);
            Hap_Entry_Fact::del_elem(he2_bptr);
            check_he(he1_bptr);
        }
        else
        {
            he2_bptr->hh_cont().splice_right(he1_bptr->hh_cont(), he2_bptr);
            _he_cont.erase(he1_bptr);
            Hap_Entry_Fact::del_elem(he1_bptr);
            check_he(he2_bptr);
        }
    }
} // Hap_Map::connect_unique
*/

void Hap_Map::dispose(Hap_Entry_BPtr he_bptr)
{
    he_bptr->hh_cont().clear_and_dispose([this] (Hap_Hop_BPtr hh_bptr) {
            this->hh_set().erase(hh_bptr);
            Hap_Hop_Fact::del_elem(hh_bptr);
        });
    Hap_Entry_Fact::del_elem(he_bptr);
}

void Hap_Map::clear_and_dispose()
{
    he_cont().clear_and_dispose([this] (Hap_Entry_BPtr he_bptr) {
            this->dispose(he_bptr);
        });
}

void Hap_Map::dump_consecutive_anchor_pair_stats(ostream& os, const Allele_Anchor& a1, const Allele_Anchor& a2) const
{
    ASSERT(a1.get_sibling(true) == a2);
    auto a1_support = a1.support();
    auto a2_support = a2.support();
    auto connect_map = Allele_Anchor::connect(a1_support, a2_support);
    auto a1_hops_rg = hh_set().equal_range(a1);
    auto a2_hops_rg = hh_set().equal_range(a2);
    auto make_terminal_hop_set = [] (const Allele_Anchor& anchor, const decltype(a1_hops_rg)& rg, bool right_anchor) {
        static_cast< void >(anchor);
        set< Hap_Hop_CBPtr > res;
        for (auto it = rg.first; it != rg.second; ++it)
        {
            auto hh_cbptr = &*it;
            ASSERT(it->allele_anchor() == anchor);
            if (Hap_Map::is_end_hop(hh_cbptr, right_anchor == hh_cbptr->c_direction()))
            {
                res.insert(hh_cbptr);
            }
        }
        return res;
    };
    auto a1_terminal_hops = make_terminal_hop_set(a1, a1_hops_rg, false);
    auto a2_terminal_hops = make_terminal_hop_set(a2, a2_hops_rg, true);
    ASSERT(a1_terminal_hops.size() <= a1_support.size());
    ASSERT(a2_terminal_hops.size() <= a2_support.size());

    Size_Type dist = Allele_Anchor::dist(a1, a2);
    os << "consec.anchors\t" << a1 << "\t" << a2 << "\t" << dist << "\t"
       << a1_support.size() << "\t" << a2_support.size() << "\t" << connect_map.size() << "\t"
       << a1_terminal_hops.size() << "\t" << a2_terminal_hops.size() << "\t";
    if (a1_support.size() == 2 and a2_support.size() == 2)
    {
        list< size_t > l;
        for (auto al1 : a1_support | ba::map_keys)
        {
            for (auto al2 : a2_support | ba::map_keys)
            {
                auto p = make_pair(al1, al2);
                if (connect_map.count(p))
                {
                    l.push_back(connect_map.at(p).size());
                }
                else
                {
                    l.push_back(0);
                }
            }
        }
        bool first = true;
        for (const auto& e : l)
        {
            if (not first) os << ",";
            first = false;
            os << e;
        }
    }
    else
    {
        os << ".";
    }
    os << endl;
}

void Hap_Map::dump_stats(ostream& os, const Graph& g) const
{
    os << "consec.anchors\ta1\ta2\tdist\tnum.alleles.a1\tnum.alleles.a2\tsize.conn\tnum.term.hops.a1\tnum.term.hops.a2\t2-2-degs" << endl;
    for (auto ce_cbptr : g.ce_cont() | referenced)
    {
        if (not ce_cbptr->is_normal()) continue;
        Allele_Anchor last_anchor(ce_cbptr, false);
        for (auto mut_cbptr : ce_cbptr->mut_cont() | referenced)
        {
            Allele_Anchor anchor(mut_cbptr);
            // last_anchor <-> anchor
            dump_consecutive_anchor_pair_stats(os, last_anchor, anchor);
            last_anchor = move(anchor);
        }
        Allele_Anchor re_anchor(ce_cbptr, true);
        // last_anchor <-> re_anchor
        dump_consecutive_anchor_pair_stats(os, last_anchor, re_anchor);
    }
}

void Hap_Map::check_he(Hap_Entry_CBPtr he_cbptr) const
{
    static_cast< void >(he_cbptr);
#ifndef DISABLE_ASSERTS
    ASSERT(he_cbptr);
    he_cbptr->hh_cont().check();
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
                // anchors must be consecutive in the given direction
                ASSERT(last_hh_cbptr->allele_anchor().get_sibling(not hh_cbptr->c_direction())
                       == hh_cbptr->allele_anchor());
            }
            else
            {
                // these must be mirror hops
                ASSERT(last_hh_cbptr->allele_anchor().is_endpoint());
                ASSERT(hh_cbptr->allele_anchor().is_endpoint());
                bool same_orientation = last_hh_cbptr->c_direction() == hh_cbptr->c_direction();
                ASSERT(last_hh_cbptr->allele_specifier() == Allele_Specifier(hh_cbptr->ce_cbptr(), same_orientation));
                ASSERT(hh_cbptr->allele_specifier() == Allele_Specifier(last_hh_cbptr->ce_cbptr(), same_orientation));
                ASSERT(last_hh_cbptr->allele_anchor().c_right() == not last_hh_cbptr->c_direction());
                auto last_hh_allele_support = last_hh_cbptr->allele_anchor().support();
                auto hh_allele_support = hh_cbptr->allele_anchor().support();
                ASSERT(last_hh_allele_support.count(last_hh_cbptr->allele_specifier()) == 1);
                ASSERT(hh_allele_support.count(hh_cbptr->allele_specifier()) == 1);
                ASSERT(last_hh_allele_support.at(last_hh_cbptr->allele_specifier()).size() > 1);
                ASSERT(last_hh_allele_support.at(last_hh_cbptr->allele_specifier()).size()
                       == hh_allele_support.at(hh_cbptr->allele_specifier()).size());
            }
        }
        last_hh_cbptr = hh_cbptr;
    }
#endif
} // Hap_Map::check_he

void Hap_Map::check(const set< Hap_Entry_CBPtr >& he_set) const
{
    static_cast< void >(he_set);
#ifndef DISABLE_ASSERTS
    for (auto he_cbptr : he_set)
    {
        check_he(he_cbptr);
    }
#endif
} // Hap_Map::check

void Hap_Map::check_all() const
{
#ifndef DISABLE_ASSERTS
    he_cont().check();
    hh_set().check();
    for (auto he_cbptr : he_cont() | referenced)
    {
        check_he(he_cbptr);
    }
#endif
} // Hap_Map::check_all


} // namespace MAC
