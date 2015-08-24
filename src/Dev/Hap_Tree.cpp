#include "Hap_Tree.hpp"

namespace MAC
{

Hap_Tree_Node *
make_hap_tree(const Hap_Hop_Set::find_anchor_type& haps)
{
    ASSERT(not haps.empty());
    Hap_Tree_Node * node = new Hap_Tree_Node;
    node->parent = nullptr;
    if (haps.size() == 1 and haps.begin()->second.size() == 1)
    {
        // single hap left
        node->he_bptr = haps.begin()->second.begin()->first->he_cbptr().unconst();
        node->h_direction = haps.begin()->second.begin()->second;
    }
    else
    {
        node->anchor = haps.begin()->second.begin()->first->allele_anchor();
        for (const auto& v : haps)
        {
            Allele_Specifier allele = v.first;
            ASSERT(not v.second.empty());
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
                child = make_hap_tree(next_haps);
            }
            child->parent = node;
            node->edge_m[allele] = child;
        }
    }
    return node;
}

} // namespace MAC
