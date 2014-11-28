#include "Hap_Map.hpp"


namespace MAC
{

Hap_Map::Hap_Map(const Graph& g)
{
    build(g);
}

map< Allele_Specifier, set< Hap_Hop_CBPtr > > extend_haps_at_endpoint(const Allele_Anchor& anchor)
{
    ASSERT(anchor.is_endpoint());
    map< Allele_Specifier, set< Hap_Hop_CBPtr > > res;
    auto anchor_support = anchor.support();
    for (auto p : anchor_support)
    {
        Allele_Specifier allele = p.first;
        ASSERT(allele.is_endpoint());
        // compute corresponding anchor at the other endpoint of the edge
        Allele_Anchor mirror_anchor(allele.ce_cbptr(), anchor.c_right() != allele.same_strand());
        // compute corresponding allele at the other endpoint
        Allele_Specifier mirror_allele(anchor.ce_cbptr(), allele.same_strand());
        // check if mirror_anchor has haplotypes
        set< Hap_Hop_BPtr > haps_to_extend;
        auto it_p = _hh_set.equal_range(mirror_anchor, detail::Hap_Hop_Comparator);
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
            Hap_Entry_BPtr he_bptr = Hap_Entry_Fact::new_elem();
            _he_cont.insert(he_bptr);
            Hap_Hop_BPtr hh_bptr = Hap_Hop_Fact::new_elem(he_bptr, anchor, allele, false);
            he_bptr->hh_cont().push_front(hh_bptr);
            _hh_set.insert(hh_bptr);
            res[allele].insert(hh_bptr);
            check_he(he_bptr);
        }
        else
        {
            // haplotypes exist at mirror allele
            // we extend them to this allele
            for (auto hh_bptr : haps_to_extend)
            {
                ASSERT(hh_bptr);
                ASSERT(hh_bptr->he_cbptr());
                ASSERT(not hh_bptr->he_cbptr().hh_cont().empty());
                ASSERT(hh_bptr == &hh_bptr->he_cbptr().hh_cont().front()
                       or hh_bptr == &hh_bptr->he_cbptr().hh_cont().back());
                Hap_Entry_BPtr he_bptr = hh_bptr->he_cbptr().unconst();
                bool is_hap_end = hh_bptr == &hh_bptr->he_cbptr().hh_cont().back();
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
    for (auto ce_cbptr : g.ce_cont())
    {
        if (not ce_cptr->is_normal()) continue;

        // construct left and right endpoint anchors
        Allele_Anchor le_anchor(ce_cbptr, false);
        Allele_Anchor re_anchor(ce_cbptr, true);
        // extend existing haplotypes into this contig
        auto le_haps = extend_haps_at_endpoint(le_anchor);
        auto re_haps = extend_haps_at_endpoint(re_anchor);
        //TODO
    }
} // Hap_Map::build

void Hap_Map::check_he(Hap_Entry_CBPtr he_cbptr) const
{
    static_cast< void >(he_cbptr);
#ifndef BOOST_DISABLE_ASSERTS
    ASSERT(he_cbptr);
    Hap_Hop_CBPtr last_hh_cbptr;
    for (auto hh_cbptr : he_cbptr.chunk_cont() | referenced)
    {
        ASSERT(hh_cbptr);
        ASSERT(hh_cbptr->he_cbptr() == he_cbptr());
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
}

void Hap_Map::check(const set< Hap_Entry_CBPtr >& he_set) const
{
    static_cast< void >(he_set);
#ifndef BOOST_DISABLE_ASSERTS
    for (auto he_cbptr : he_set)
    {
        check_he(he_cbptr);
    }
#endif
}

void Hap_Map::check_all() const
{
#ifndef BOOST_DISABLE_ASSERTS
    for (auto he_cbptr : _he_cont | referenced)
    {
        check_he(he_cbptr);
    }
#endif
}


} // namespace MAC
