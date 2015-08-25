#ifndef __HAP_TREE_HPP
#define __HAP_TREE_HPP

#include <vector>
#include <map>
#include "Hap_Map.hpp"

namespace MAC
{

struct Hap_Tree_Node
{
    Hap_Tree_Node * parent;
    Allele_Anchor anchor;
    map< Allele_Specifier, Hap_Tree_Node * > edge_m;
    Hap_Entry_BPtr he_bptr;
    bool h_direction;
}; // class Hap_Tree_Node

Hap_Tree_Node *
make_hap_tree(const Hap_Hop_Set::find_anchor_type& haps, const Allele_Anchor::read_support_type& anchor_support);

} // namespace MAC

#endif
