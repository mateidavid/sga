#include "Allele_Anchor.hpp"
#include "Contig_Entry.hpp"

namespace MAC
{

auto
Allele_Anchor::connect(
    const Anchor_Read_Support& a1_support,
    const Anchor_Read_Support& a2_support,
    bool same_st)
    -> anchor_connect_type
{
    anchor_connect_type res;
    for (const auto& p1 : a1_support)
    {
        for (const auto& p2 : a2_support)
        {
            Allele_Read_Support s;
            for (int d = 0; d < 2; ++d)
            {
                set_intersection(p1.second[d].begin(), p1.second[d].end(),
                                 p2.second[d == same_st].begin(), p2.second[d == same_st].end(),
                                 inserter(s[d], s[d].end()));
            }
            res.insert(make_pair(make_pair(p1.first, p2.first), move(s)));
        }
    }
    return res;
} // Allele_Anchor::connect

Allele_Anchor Allele_Anchor::get_sibling(bool c_direction) const
{
    if (is_endpoint())
    {
        // cannot go past left or right endpoints
        ASSERT(c_right() == c_direction);
        if (not c_right())
        {
            // anchor is left endpoint; next is first mutation (if one exists)
            if (not ce_cbptr()->mut_cont().empty())
            {
                auto mut_cbptr = &*ce_cbptr()->mut_cont().begin();
                return Allele_Anchor(mut_cbptr, ce_cbptr());
            }
            else
            {
                return Allele_Anchor(ce_cbptr(), true);
            }
        }
        else
        {
            // anchor is right endpoint; next is last mutation (if one exists)
            if (not ce_cbptr()->mut_cont().empty())
            {
                auto mut_cbptr = &*ce_cbptr()->mut_cont().rbegin();
                return Allele_Anchor(mut_cbptr, ce_cbptr());
            }
            else
            {
                return Allele_Anchor(ce_cbptr(), false);
            }
        }
    }
    else // current anchor is a mutation
    {
        auto it = ce_cbptr()->mut_cont().iterator_to(*mut_cbptr());
        if (c_direction)
        {
            // find previous mutation if one exists
            if (it != ce_cbptr()->mut_cont().begin())
            {
                --it;
                return Allele_Anchor(&*it, ce_cbptr());
            }
            else
            {
                return Allele_Anchor(ce_cbptr(), false);
            }
        }
        else
        {
            // find next mutation if one exists
            ++it;
            if (it != ce_cbptr()->mut_cont().end())
            {
                return Allele_Anchor(&*it, ce_cbptr());
            }
            else
            {
                return Allele_Anchor(ce_cbptr(), true);
            }
        }
    }
} // Allele_Anchor::get_sibling

Anchor_Chunk_Support
Allele_Anchor::chunk_support(unsigned min_edge_support) const
{
    min_edge_support = max(min_edge_support, 1u);
    auto orient_by_strand = [] (Read_Chunk_CBPtr rc_cbptr) { return rc_cbptr->get_rc(); };
    Anchor_Chunk_Support res;
    if (is_mutation())
    {
        set< Read_Chunk_CBPtr > qr_set;
        set< Read_Chunk_CBPtr > full_rf_set;
        tie(qr_set, full_rf_set, ignore) = ce_cbptr()->mut_support(mut_cbptr());
        res.insert(make_pair(Allele_Specifier(false),
                             Allele_Chunk_Support(full_rf_set, orient_by_strand)));
        res.insert(make_pair(Allele_Specifier(true),
                             Allele_Chunk_Support(qr_set, orient_by_strand)));
    }
    else // is_endpoint
    {
        auto m = ce_cbptr()->out_chunks_dir(c_right(), 3, min_edge_support - 1);
        if (not m.empty())
        {
            for (auto& p : m)
            {
                res.insert(make_pair(Allele_Specifier(p.first),
                                     Allele_Chunk_Support(p.second, orient_by_strand)));
            }
        }
        else
        {
            // out-degree 0
            auto pos = c_right()? ce_cbptr()->len() : 0;
            auto it_r = ce_cbptr()->chunk_cont().iintersect(pos, pos) | referenced;
            ASSERT(it_r.begin() != it_r.end());
            res.insert(make_pair(Allele_Specifier(nullptr, true),
                                 Allele_Chunk_Support(it_r.begin(), it_r.end(), orient_by_strand)));
        }
    }
    return res;
} // Allele_Anchor::chunk_support

Anchor_Read_Support
Allele_Anchor::read_support(unsigned min_edge_support) const
{
    Anchor_Read_Support res;
    auto cs = chunk_support(min_edge_support);
    for (auto& p : cs)
    {
        res.insert(make_pair(p.first, move(p.second.reads())));
    }
    return res;
} // Allele_Anchor::read_support

/// Distance between consecutive allele anchors.
Size_Type Allele_Anchor::dist(const Allele_Anchor& lhs, const Allele_Anchor& rhs)
{
    ASSERT(lhs.get_sibling(false) == rhs);
    return (rhs.is_endpoint()? rhs.ce_cbptr()->len() : rhs.mut_cbptr()->rf_start())
        - (lhs.is_endpoint()? 0 : lhs.mut_cbptr()->rf_end());
} // Allele_Anchor::dist

/// Comparator for storage a tree.
bool operator < (const Allele_Anchor& lhs, const Allele_Anchor& rhs)
{
    // order by Contig_Entry bptr first
    if (lhs.ce_cbptr() != rhs.ce_cbptr())
    {
        return lhs.ce_cbptr() < rhs.ce_cbptr();
    }
    // next, order left endpoint before mutations before right endpoint
    else if (lhs.is_endpoint())
    {
        if (rhs.is_endpoint())
        {
            return lhs.c_right() < rhs.c_right();
        }
        else
        {
            return lhs.c_right() == false;
        }
    }
    else // lhs.is_mutation()
    {
        if (rhs.is_endpoint())
        {
            return rhs.c_right() == true;
        }
        else
        {
            // order by Mutation rf_start, then by bptr
            if (lhs.mut_cbptr()->rf_start() != rhs.mut_cbptr()->rf_start())
            {
                return lhs.mut_cbptr()->rf_start() < rhs.mut_cbptr()->rf_start();
            }
            else
            {
                return lhs.mut_cbptr() < rhs.mut_cbptr();
            }
        }
    }
} // operator <

bool operator == (const Allele_Anchor& lhs, const Allele_Anchor& rhs)
{
    if (lhs.is_endpoint() != rhs.is_endpoint())
    {
        return false;
    }
    else if (lhs.is_endpoint()) // both endpoints
    {
        return lhs.ce_cbptr() == rhs.ce_cbptr() and lhs.c_right() == rhs.c_right();
    }
    else // both mutations
    {
        return lhs.mut_cbptr() == rhs.mut_cbptr();
    }
} // operator ==

ptree Allele_Anchor::to_ptree() const
{
    if (is_endpoint())
    {
        return ptree().put("is_endpoint", true)
            .put("ce_cbptr", ce_cbptr().to_int())
            .put("c_right", c_right());
    }
    else
    {
        return ptree().put("is_endpoint", false)
            .put("ce_cbptr", ce_cbptr().to_int())
            .put("mut_cbptr", mut_cbptr().to_int());
    }
} // Allele_Anchor::to_ptree

ostream& operator << (ostream& os, const Allele_Anchor& a)
{
    if (a.is_endpoint())
    {
        os << "(" << a.ce_cbptr().to_int() << "," << a.c_right() << ")";
    }
    else
    {
        os << a.mut_cbptr().to_int();
    }
    return os;
} // operator <<

} // namespace MAC
