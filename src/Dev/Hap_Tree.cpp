#include "Hap_Tree.hpp"

namespace MAC
{

Hap_Tree_Node *
make_hap_tree(const Hap_Hop_Set::find_anchor_type& haps, const Allele_Anchor::read_support_type& anchor_support)
{
    ASSERT(not haps.empty());
    ASSERT(not ba::values(haps).begin()->empty());
    Hap_Tree_Node * node = new Hap_Tree_Node;
    node->parent = nullptr;
    node->anchor = ba::values(haps).begin()->begin()->first->allele_anchor();
    for (const auto& v : haps)
    {
        Allele_Specifier allele = v.first;
        ASSERT(not v.second.empty());
        const array< set< Read_Entry_CBPtr >, 2 >& allele_support = anchor_support.at(v.first);
        ASSERT(not allele_support[0].empty() or not allele_support[1].empty());
        Hap_Tree_Node * child = nullptr;
        if (v.second.size() == 1)
        {
            // no more branching past this allele
            child = new Hap_Tree_Node;
            child->he_bptr = v.second.begin()->first->he_cbptr().unconst();
            child->h_direction = v.second.begin()->second;
        }
        else
        {
            // compute next hops
            Hap_Hop_Set::find_anchor_type next_haps;
            for (const auto& p : v.second)
            {
                Hap_Hop_CBPtr next_hh_cbptr = p.first->he_cbptr()->hh_cont().get_sibling(p.first, p.second);
                ASSERT(next_hh_cbptr);
                next_haps[next_hh_cbptr->allele_specifier()].insert(make_pair(next_hh_cbptr, p.second));
            }
            ASSERT(all_of(
                       next_haps.begin(),
                       next_haps.end(),
                       [] (const Hap_Hop_Set::find_anchor_type::value_type& v) {
                           return all_of(
                               v.second.begin(),
                               v.second.end(),
                               [] (const pair< Hap_Hop_CBPtr, bool >& p) {
                                   return p.first->allele_anchor() == next_anchor;
                               });
                       }));

            child = make_hap_tree(next_haps);
        }
        child->parent = node;
        node->edge_m[allele] = child;
    }
    return node;
}

} // namespace MAC
