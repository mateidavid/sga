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
    Hap_Map(const Graph& g);

private:
    Hap_Hop_BPtr make_single_allele_hop(const Allele_Anchor& anchor, const Allele_Specifier& allele,
                                        bool direction = false);
    map< Allele_Specifier, set< Hap_Hop_CBPtr > > make_mutation_haps(Mutation_CBPtr mut_cbptr);
    map< Allele_Specifier, set< Hap_Hop_CBPtr > > extend_endpoint_haps(const Allele_Anchor& anchor);
    void connect_unique(const Allele_Anchor& a1, const Allele_Anchor& a2,
                        map< Allele_Specifier, set< Hap_Hop_CBPtr > >& a1_hops,
                        map< Allele_Specifier, set< Hap_Hop_CBPtr > >& a2_hops);
    void build(const Graph& g);

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
