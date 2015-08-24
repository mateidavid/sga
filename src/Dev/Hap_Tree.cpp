#include "Hap_Tree.hpp"

Hap_Tree_Node *
make_hap_tree(const Allele_Anchor& anchor, bool c_direction,
              const Hap_Hop_Set::find_anchor_type& haps)
{
    ASSERT(not haps.empty());
    Hap_Tree_Node * node = new Hap_Tree_Node;
    node->parent = nullptr;
    if (haps.size() == 1 and haps.begin()->second.size() == 1)
    {
        // single hap left
        node->he_bptr = haps.begin()->second.begin()->he_cbptr().unconst();
    }
    else
    {
        for (const auto& p : haps)
        {
            Allele_Specifier allele = p.first;
            ASSERT(not p.second.empty());
            if (p.second.size() == 1)
            {
                // no more branching past this allele
                Hap_Tree_Node * child = new Hap_Tree_Node;
                child->parent = node;
                child->he_bptr = p.second.begin()->he_cbptr().unconst();
                node->edge_m[allele] = child;
            }
            else
            {
                // compute next hops
                Hap_Hop_Set::find_anchor_type next_haps;
                for (auto hh_cbptr : p.second)
                {
                    Hap_Hop_CBPtr next_hh_cbptr = hh_cbptr->he_cbptr()->hh_cont().get_sibling(
                        hh_cbptr, hh_cbptr->c_direction() != c_direction);
                    ASSERT(next_hh_cbptr);
                    
                }

            }
        }
    }
    return res;
}
