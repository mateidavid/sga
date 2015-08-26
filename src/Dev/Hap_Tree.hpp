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
    Allele_Read_Support allele_support;
    set< pair< Hap_Hop_CBPtr, bool > > hops;
}; // class Hap_Tree_Node

typedef map< Hap_Tree_Node *, set< Hap_Tree_Node * > > Hap_Tree_Connections;

Hap_Tree_Node *
make_hap_tree(const Hap_Hop_Set::find_anchor_type& haps, const Anchor_Read_Support& anchor_support);

bool is_branch_end(Hap_Tree_Node * node);
bool is_branch_start(Hap_Tree_Node * node);
Hap_Tree_Node * get_branch_end(Hap_Tree_Node * node);
Hap_Tree_Node * get_branch_start(Hap_Tree_Node * node);
vector< Hap_Tree_Node * > get_branch_start_descendants(Hap_Tree_Node * node);
vector< Hap_Tree_Node * > get_branch_start_ancestors(Hap_Tree_Node * node);
Hap_Tree_Node * get_least_common_branch_start_ancestor(Hap_Tree_Node * n1, Hap_Tree_Node * n2);

void
add_root_connections(Hap_Tree_Node * a1_root,
                     Hap_Tree_Node * a2_root,
                     const vector< Hap_Tree_Node * >& a1_branch_starts,
                     const vector< Hap_Tree_Node * >& a2_branch_starts,
                     Hap_Tree_Connections& a1_to_a2_conn,
                     Hap_Tree_Connections& a2_to_a1_conn);

void
add_branch_connections(const vector< Hap_Tree_Node * >& a1_branch_starts,
                       const vector< Hap_Tree_Node * >& a2_branch_starts,
                       Hap_Tree_Connections& a1_to_a2_conn,
                       Hap_Tree_Connections& a2_to_a1_conn,
                       unsigned min_re_support);

void
erase_contained_connections(Hap_Tree_Connections& a1_to_a2_conn,
                            Hap_Tree_Connections& a2_to_a1_conn);

void
expand_overlapping_connections(Hap_Tree_Connections& a1_to_a2_conn,
                               Hap_Tree_Connections& a2_to_a1_conn);

void
expand_connections_branch_ends(Hap_Tree_Connections& a1_to_a2_conn,
                               Hap_Tree_Connections& a2_to_a1_conn);

void
implement_connections(Hap_Tree_Node * a1_root,
                      Hap_Tree_Node * a2_root,
                      const Hap_Tree_Connections& a1_to_a2_conn,
                      const Hap_Tree_Connections& a2_to_a1_conn,
                      Hap_Map& hm);

} // namespace MAC

#endif
