#include "Hap_Tree.hpp"
#include "Anchor_Support.hpp"

namespace MAC
{

Hap_Tree_Node *
make_hap_tree(const Hap_Hop_Set::find_anchor_type& haps, const Anchor_Read_Support& anchor_support)
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
        auto allele_support = get_allele_read_support(move(anchor_support), v.first);
        Hap_Tree_Node * child = nullptr;
        if (v.second.size() == 1)
        {
            // no more branching past this allele
            child = new Hap_Tree_Node;
            child->allele_support = allele_support;
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
            Allele_Anchor next_anchor = ba::values(next_haps).begin()->begin()->first->allele_anchor();
            ASSERT(all_of(
                       ba::values(next_haps),
                       [&] (const Hap_Hop_Set::find_anchor_type::mapped_type& s) {
                           return all_of(
                               s,
                               [&] (const pair< Hap_Hop_CBPtr, bool >& p) {
                                   return p.first->allele_anchor() == next_anchor;
                               });
                       }));
            // compute support of next anchor, remove reads not reaching current allele
            auto next_anchor_support = get_hop_read_support(*ba::values(next_haps).begin()->begin());
            subset_support(next_anchor_support, allele_support, false);
            child = make_hap_tree(next_haps, next_anchor_support);
        }
        child->parent = node;
        child->hops = v.second;
        node->edge_m[allele] = child;
    }
    return node;
}

void check(Hap_Tree_Node * node)
{
#ifndef DISABLE_ASSERTS
    ASSERT(node);
    ASSERT((node->edge_m.size() == 0) == (node->hops.size() == 1));
#endif
}

bool is_branch_end(Hap_Tree_Node * node)
{
    check(node);
    // root is a branch end
    if (not node->parent) return true;
    // otherwise, node is a branch end iff degree != 1
    return node->edge_m.size() != 1;
}

bool is_branch_start(Hap_Tree_Node * node)
{
    check(node);
    // root is a branch start
    if (not node->parent) return true;
    // otherwise, this is a branch start iff parent is a branch end
    return is_branch_end(node->parent);
}

Hap_Tree_Node *
get_branch_end(Hap_Tree_Node * node)
{
    check(node);
    // root is a branch end
    while (not is_branch_end(node))
    {
        ASSERT(node->edge_m.size() == 1);
        node = node->edge_m.begin()->second;
        check(node);
    }
    return node;
}

Hap_Tree_Node *
get_branch_start(Hap_Tree_Node * node)
{
    check(node);
    while (not is_branch_start(node))
    {
        node = node->parent;
        check(node);
    }
    return node;
}


vector< Hap_Tree_Node * >
get_branch_start_descendants(Hap_Tree_Node * node)
{
    check(node);
    node = get_branch_end(node);
    vector< Hap_Tree_Node * > res;
    for (const auto child : node->edge_m | ba::map_values)
    {
        ASSERT(is_branch_start(child));
        res.push_back(child);
        auto tmp = get_branch_start_descendants(child);
        res.insert(res.end(), tmp.begin(), tmp.end());
    }
    return res;
}

vector< Hap_Tree_Node * >
get_branch_start_ancestors(Hap_Tree_Node * node)
{
    check(node);
    vector< Hap_Tree_Node * > res;
    if (not node->parent)
    {
        return res;
    }
    node = get_branch_start(node->parent);
    while (node->parent)
    {
        ASSERT(is_branch_start(node));
        res.push_back(node);
        node = get_branch_start(node->parent);
    }
    return res;
}

Hap_Tree_Node *
get_least_common_branch_start_ancestor(Hap_Tree_Node * n1, Hap_Tree_Node * n2)
{
    auto n1_ancestors = get_branch_start_ancestors(n1);
    auto n2_ancestors = get_branch_start_ancestors(n2);
    auto it1 = n1_ancestors.rbegin();
    auto it2 = n2_ancestors.rbegin();
    while (it1 != n1_ancestors.rend()
           and it1 != n2_ancestors.rend()
           and *it1 == *it2)
    {
        ++it1;
        ++it2;
    }
    if (it1 != n1_ancestors.rbegin())
    {
        return *prev(it1);
    }
    else
    {
        return (*n1_ancestors.rbegin())->parent;
    }
}

void
add_root_connections(Hap_Tree_Node * a1_root,
                     Hap_Tree_Node * a2_root,
                     const vector< Hap_Tree_Node * >& a1_branch_starts,
                     const vector< Hap_Tree_Node * >& a2_branch_starts,
                     Hap_Tree_Connections& a1_to_a2_conn,
                     Hap_Tree_Connections& a2_to_a1_conn)
{
    auto make_root_connections = [] (Hap_Tree_Node * b1_root,
                                     const vector< Hap_Tree_Node * >& b2_branch_starts,
                                     Hap_Tree_Connections& b1_to_b2_conn,
                                     Hap_Tree_Connections& b2_to_b1_conn) {
        if (b1_root->edge_m.size() == 1)
        {
            auto node1 = b1_root->edge_m.begin()->second;
            for (auto node2 : b2_branch_starts)
            {
                // add connection b1:node1 <-> b2:node2
                b1_to_b2_conn[node1].insert(node2);
                b2_to_b1_conn[node2].insert(node1);
            }
        }
    };
    make_root_connections(a1_root, a2_branch_starts, a1_to_a2_conn, a2_to_a1_conn);
    make_root_connections(a2_root, a1_branch_starts, a2_to_a1_conn, a1_to_a2_conn);
}

void
add_branch_connections(const vector< Hap_Tree_Node * >& a1_branch_starts,
                       const vector< Hap_Tree_Node * >& a2_branch_starts,
                       Hap_Tree_Connections& a1_to_a2_conn,
                       Hap_Tree_Connections& a2_to_a1_conn,
                       unsigned min_re_support)
{
    for (auto node1 : a1_branch_starts)
    {
        for (auto node2 : a2_branch_starts)
        {
            auto common_support = common_allele_support(node1->allele_support, node2->allele_support, true);
            if (common_support[0].size() + common_support[1].size() >= min_re_support)
            {
                // add connection node1<->node2
                a1_to_a2_conn[node1].insert(node2);
                a2_to_a1_conn[node2].insert(node1);
            }
        }
    }
}

bool connection_contained(Hap_Tree_Node * node1, Hap_Tree_Node * node2, const Hap_Tree_Connections& a1_to_a2_conn)
{
    auto node1_branch_end = get_branch_end(node1);
    auto node2_branch_end = get_branch_end(node2);
    auto node1_descendants = get_branch_start_descendants(node1_branch_end);
    node1_descendants.push_back(node1);
    auto node2_descendants = get_branch_start_descendants(node2_branch_end);
    node1_descendants.push_back(node2);
    for (auto u1 : node1_descendants)
    {
        for (auto u2 : node2_descendants)
        {
            if (u1 == node1 and u2 == node2) continue;
            if (a1_to_a2_conn.count(u1) and a1_to_a2_conn.at(u1).count(u2))
            {
                return true;
            }
        }
    }
    return false;
}

void
erase_contained_connections(Hap_Tree_Connections& a1_to_a2_conn,
                            Hap_Tree_Connections& a2_to_a1_conn)
{
    for_each_it_advance(
        a1_to_a2_conn.begin(),
        a1_to_a2_conn.end(),
        [&] (Hap_Tree_Connections::iterator it)
        {
            auto node1 = it->first;
            for_each_advance(
                it->second.begin(),
                it->second.end(),
                [&] (Hap_Tree_Node * node2)
                {
                    // if connection node1<->node2 is completely contained, remove it
                    if (connection_contained(node1, node2, a1_to_a2_conn))
                    {
                        it->second.erase(node2);
                        a2_to_a1_conn.at(node2).erase(node1);
                        if (a2_to_a1_conn.at(node2).empty())
                        {
                            a2_to_a1_conn.erase(node2);
                        }
                    }
                });
            if (it->second.empty())
            {
                a1_to_a2_conn.erase(node1);
            }
        });
}

bool expand_connection(Hap_Tree_Node * n1,
                       Hap_Tree_Node * n2,
                       Hap_Tree_Connections& a1_to_a2_conn,
                       Hap_Tree_Connections& a2_to_a1_conn)
{
    // look for connections between an ancestors of n1 and descendants of n2
    auto n1_ancestors = get_branch_start_ancestors(n1);
    auto n2_descendants_v = get_branch_start_descendants(n2);
    set< Hap_Tree_Node * > n2_descendants(n2_descendants_v.begin(), n2_descendants_v.end());
    Hap_Tree_Node * m = nullptr;
    for (auto u1 : n1_ancestors)
    {
        if (a1_to_a2_conn.count(u1) == 0) continue;
        for (auto u2 : a1_to_a2_conn.at(u1))
        {
            if (n2_descendants.count(u2) == 0) continue;
            // found connection u1<->u2
            if (m)
            {
                m = get_least_common_branch_start_ancestor(m, u2);
                if (m == n2) return false;
            }
            else
            {
                m = u2;
            }
        }
    }
    // expand connection n1->n2 to n1->m
    ASSERT(n2_descendants.count(m));
    a1_to_a2_conn.at(n1).erase(n2);
    a1_to_a2_conn.at(n1).insert(m);
    a2_to_a1_conn.at(n2).erase(n1);
    if (a2_to_a1_conn.at(n2).empty()) a2_to_a1_conn.erase(n2);
    a2_to_a1_conn[m].insert(n1);
    return true;
}

void
expand_overlapping_connections(Hap_Tree_Connections& a1_to_a2_conn,
                               Hap_Tree_Connections& a2_to_a1_conn)
{
    while (bool done = false)
    {
        done = true;
        for (auto node1 : a1_to_a2_conn | ba::map_keys)
        {
            for (auto node2 : a1_to_a2_conn.at(node1))
            {
                if (expand_connection(node1, node2, a1_to_a2_conn, a2_to_a1_conn)
                    or expand_connection(node1, node2, a1_to_a2_conn, a2_to_a1_conn))
                {
                    done = false;
                    erase_contained_connections(a1_to_a2_conn, a2_to_a1_conn);
                    break;
                }
            }
            if (not done) break;
        }
    }
}

map< Hap_Tree_Node *, set< pair< Hap_Hop_CBPtr, bool > > >
get_initial_node_haps(Hap_Tree_Node * a1_root,
                      const Hap_Tree_Connections& a1_to_a2_conn,
                      Hap_Map& hm)
{
    map< Hap_Tree_Node *, set< pair< Hap_Hop_CBPtr, bool > > > res;
    for (auto n1 : a1_to_a2_conn | ba::map_keys)
    {
        auto m1 = get_branch_end(n1);
        ASSERT(not m1->hops.empty());
        Hap_Hop_CBPtr crt_hh_cbptr;
        bool h_direction;
        tie(crt_hh_cbptr, h_direction) = *m1->hops.begin();
        if (m1->hops.size() > 1)
        {
            // for internal branch_end nodes, we create new haplotypes
            crt_hh_cbptr = hm.duplicate_haplotype(crt_hh_cbptr, not h_direction);
        }
        Hap_Hop_CBPtr head_hh_cbptr = (not h_direction
                                       ? crt_hh_cbptr->he_cbptr()->hh_cont().front()
                                       : crt_hh_cbptr->he_cbptr()->hh_cont().back());
        ASSERT(head_hh_cbptr->allele_anchor() == a1_root->anchor);
        res[n1].insert(make_pair(head_hh_cbptr, h_direction));
    }
    return res;
}

void
multiply_node_haps(map< Hap_Tree_Node *, set< pair< Hap_Hop_CBPtr, bool > > >& haps_m,
                   const Hap_Tree_Connections& a1_to_a2_conn,
                   Hap_Map& hm)
{
    ASSERT(a1_to_a2_conn.size() == haps_m.size());
    for (auto n1 : a1_to_a2_conn | ba::map_keys)
    {
        unsigned cnt = a1_to_a2_conn.at(n1).size();
        ASSERT(not haps_m.at(n1).empty());
        Hap_Hop_CBPtr hh_cbptr;
        bool h_direction;
        tie(hh_cbptr, h_direction) = *haps_m.at(n1).begin();
        cnt -= haps_m.at(n1).size();
        while (cnt)
        {
            Hap_Hop_CBPtr new_hh_cbptr = hm.duplicate_haplotype(hh_cbptr, h_direction);
            haps_m.at(n1).insert(make_pair(new_hh_cbptr, h_direction));
            --cnt;
        }
        ASSERT(a1_to_a2_conn.at(n1).size() == haps_m.at(n1).size());
    }
}

void
implement_connections(Hap_Tree_Node * a1_root,
                      Hap_Tree_Node * a2_root,
                      const Hap_Tree_Connections& a1_to_a2_conn,
                      const Hap_Tree_Connections& a2_to_a1_conn,
                      Hap_Map& hm)
{
    // first: create haplotypes at connected nodes
    auto a1_node_haps = get_initial_node_haps(a1_root, a1_to_a2_conn, hm);
    auto a2_node_haps = get_initial_node_haps(a2_root, a2_to_a1_conn, hm);
    // second: multiply haplotypes to match connection number
    multiply_node_haps(a1_node_haps, a1_to_a2_conn, hm);
    multiply_node_haps(a2_node_haps, a2_to_a1_conn, hm);
    // finally: glue haplotypes
    for (auto n1 : a1_to_a2_conn | ba::map_keys)
    {
        for (auto n2 : a1_to_a2_conn.at(n1))
        {
            // implement n1<->n2
            auto p1 = *a1_node_haps.at(n1).begin();
            a1_node_haps.at(n1).erase(p1);
            if (a1_node_haps.at(n1).empty()) a1_node_haps.erase(n1);
            auto p2 = *a2_node_haps.at(n2).begin();
            a2_node_haps.at(n2).erase(p2);
            if (a2_node_haps.at(n2).empty()) a2_node_haps.erase(n2);
            hm.merge_haps(p1.first.unconst(), p2.first.unconst(), false, false);
        }
    }
    ASSERT(a1_node_haps.empty());
    ASSERT(a2_node_haps.empty());
}

} // namespace MAC
