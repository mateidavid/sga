#ifndef __HAP_MAP_HPP
#define __HAP_MAP_HPP

#include "Graph.hpp"
#include "Hap_Entry_Cont.hpp"
#include "Hap_Hop_Set.hpp"


namespace MAC
{

class Hap_Map
{
public:
    Hap_Map(const Graph& g) { build(g); }

    DEFAULT_DEF_CTOR(Hap_Map)
    DELETE_COPY_CTOR(Hap_Map)
    DELETE_MOVE_CTOR(Hap_Map)
    DELETE_COPY_ASOP(Hap_Map)
    DELETE_MOVE_ASOP(Hap_Map)

    ~Hap_Map()
    {
        ASSERT(he_cont().empty());
        ASSERT(hh_set().empty());
    }

    GETTER(Hap_Entry_Cont, he_cont, _he_cont)
    GETTER(Hap_Hop_Set, hh_set, _hh_set)

    void build(const Graph& g);

    void dispose(Hap_Entry_BPtr he_bptr);
    void clear_and_dispose();

    static bool is_end_hop(Hap_Hop_CBPtr hh_cbptr, bool right_end)
    {
        ASSERT(hh_cbptr);
        ASSERT(hh_cbptr->he_cbptr());
        ASSERT(not hh_cbptr->he_cbptr()->hh_cont().empty());
        return hh_cbptr == (right_end? hh_cbptr->he_cbptr()->hh_cont().back()
                            : hh_cbptr->he_cbptr()->hh_cont().front());
    }
    static bool is_end_hop(Hap_Hop_CBPtr hh_cbptr)
    {
        return is_end_hop(hh_cbptr, false) or is_end_hop(hh_cbptr, true);
    }

    void dump_stats(ostream& os, const Graph& g) const;

private:
    Hap_Hop_BPtr make_single_allele_hop(const Allele_Anchor& anchor, const Allele_Specifier& allele,
                                        bool direction = false);
    map< Allele_Specifier, set< Hap_Hop_CBPtr > > make_mutation_haps(Mutation_CBPtr mut_cbptr);
    map< Allele_Specifier, set< Hap_Hop_CBPtr > > extend_endpoint_haps(const Allele_Anchor& anchor);
    void connect_unique(const Allele_Anchor& a1, const Allele_Anchor& a2,
                        map< Allele_Specifier, set< Hap_Hop_CBPtr > >& a1_hops,
                        map< Allele_Specifier, set< Hap_Hop_CBPtr > >& a2_hops);

    void dump_consecutive_anchor_pair_stats(ostream& os, const Allele_Anchor& a1, const Allele_Anchor& a2) const;

    void check_he(Hap_Entry_CBPtr) const;
    void check(const set< Hap_Entry_CBPtr >&) const;
    void check_all() const;

    Hap_Entry_Fact _he_fact;
    Hap_Hop_Fact _hh_fact;

    Hap_Entry_Cont _he_cont;
    Hap_Hop_Set _hh_set;

}; // class Hap_Map

} // namespace MAC


#endif
