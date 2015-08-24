#ifndef __HAP_TREE_HPP
#define __HAP_TREE_HPP

struct Hap_Tree_Node
{
    Hap_Tree_Node * parent;
    Allele_Anchor anchor;
    bool c_direction;
    map< Allele_Specifier, Hap_Tree_Node * > edge_m;
    Hap_Entry_BPtr he_bptr;
}; // class Hap_Tree_Node


Hap_Tree_Node *
make_hap_tree(const Allele_Anchor& anchor,
              const map< Allele_Specifier, vector< Hap_Hop_CBPtr > >& haps,
              bool c_direction);

#endif
