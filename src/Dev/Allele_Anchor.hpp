#ifndef __ALLELE_ANCHOR_HPP
#define __ALLELE_ANCHOR_HPP

#include "Mutation.hpp"
#include "Contig_Entry.hpp"
#include "Allele_Specifier.hpp"

namespace MAC
{

/** An allele anchor is either a Mutation or a Contig_Entry edge */
class Allele_Anchor
{
public:
    DEFAULT_DEF_CTOR(Allele_Anchor)
    DEFAULT_COPY_CTOR(Allele_Anchor)
    DEFAULT_MOVE_CTOR(Allele_Anchor)
    DEFAULT_COPY_ASOP(Allele_Anchor)
    DEFAULT_MOVE_ASOP(Allele_Anchor)

    explicit Allele_Anchor(Mutation_CBPtr mut_cbptr)
        : _mut_cbptr(mut_cbptr)
    {
        ASSERT(mut_cbptr);
        ASSERT(not mut_cbptr->chunk_ptr_cont().empty());
        _ce_cbptr = mut_cbptr->chunk_ptr_cont().begin()->chunk_cbptr()->ce_bptr();
    }

    Allele_Anchor(Contig_Entry_CBPtr ce_cbptr, bool c_right)
        : _ce_cbptr(ce_cbptr), _c_right(c_right)
    {
        ASSERT(ce_cbptr);
    }

    bool is_mutation() const { return _mut_cbptr; }
    bool is_endpoint() const { return not is_mutation(); }

    GETTER(Contig_Entry_CBPtr, ce_cbptr, _ce_cbptr)
    GETTER(Mutation_CBPtr, mut_cbptr, _mut_cbptr)
    GETTER(bool, c_right, _c_right)

    /** Group read chunks supporting various alleles for this anchor.
     */
    typedef map< Allele_Specifier, set< Read_Chunk_CBPtr > > allele_support_type;
    allele_support_type support() const
    {
        allele_support_type res;
        if (is_mutation())
        {
            auto p = ce_cbptr()->mut_support(mut_cbptr());
            res[Allele_Specifier(false)] = move(p.first);
            res[Allele_Specifier(true)] = move(p.second);
        }
        else // is_endpoint
        {
            auto m = ce_cbptr()->out_chunks_dir(c_right(), 3, 1);
            for (auto p : m)
            {
                res[Allele_Specifier(p.first)] = move(p.second);
            }
        }
        return res;
    }

    /** Find Read_Entry objects supporting each pairs of alleles at the given anchors.
     * For every pair of alleles at the given anchors, find all Read_Entry objects
     * observing those 2 alleles on the same or different strand.
     * @param a1_support_m Support map for anchor a1.
     * @param a2_support_m Support map for anchor a2.
     * @param same_st Bool; if true, find reads supporting pairs of alleles on the same strand;
     * if false, on different strands.
     * NOTE: Using same_st==false only makes sense when a1 & a2 are on different contigs.
     * @return A map with keys: are pairs of alleles, values: sets of Read_Entry objects.
     */
    typedef map< pair< Allele_Specifier, Allele_Specifier >, set< Read_Entry_CBPtr > > anchor_connect_type;
    static anchor_connect_type
    connect(const allele_support_type& a1_support_m, const allele_support_type& a2_support_m, bool same_st = true)
    {
        anchor_connect_type res;
        // for each anchor, construct a map of the form
        //   [ Allele_Specifier ] -> set< [ Read_Entry, strand ] >
        auto make_allele_spec_re_set_map = [] (const allele_support_type& m) {
            map< Allele_Specifier, set< pair< Read_Entry_CBPtr, bool > > > r;
            for (const auto& p : m)
            {
                Allele_Specifier allele_spec = p.first;
                auto rg = ba::transform(p.second, [] (Read_Chunk_CBPtr rc_cbptr) {
                        return make_pair(rc_cbptr->re_bptr(), rc_cbptr->get_rc());
                    });
                r.insert(make_pair(allele_spec, set< pair< Read_Entry_CBPtr, bool > >(rg.begin(), rg.end())));
            }
            return r;
        };
        auto a1_allele_spec_re_set_map = make_allele_spec_re_set_map(a1_support_m);
        auto a2_allele_spec_re_set_map = make_allele_spec_re_set_map(a2_support_m);
        // consider every pair of alleles
        for (const auto& p1 : a1_allele_spec_re_set_map)
        {
            for (const auto& p2 : a2_allele_spec_re_set_map)
            {
                // iterate through Read_Entry objects supporting the allele at a1 (p1.first),
                // find those supporting the allele at a2 (p2.first)
                for (const auto& q : p1.second)
                {
                    Read_Entry_CBPtr re_cbptr = q.first;
                    bool a1_neg_st = q.second;
                    bool a2_neg_st = (a1_neg_st == same_st);
                    if (p2.second.count(make_pair(re_cbptr, a2_neg_st)) > 0)
                    {
                        res[make_pair(p1.first, p2.first)].insert(re_cbptr);
                    }
                }
            }
        }
        return res;
    }

    /// Comparator for storage a tree.
    friend bool operator < (const Allele_Anchor& lhs, const Allele_Anchor& rhs)
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
    }
    friend bool operator == (const Allele_Anchor& lhs, const Allele_Anchor& rhs)
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
    }
    friend bool operator != (const Allele_Anchor& lhs, const Allele_Anchor& rhs) { return !(lhs == rhs); }
    friend bool operator <= (const Allele_Anchor& lhs, const Allele_Anchor& rhs) { return lhs == rhs or lhs < rhs; }
    friend bool operator >  (const Allele_Anchor& lhs, const Allele_Anchor& rhs) { return !(lhs <= rhs); }
    friend bool operator >= (const Allele_Anchor& lhs, const Allele_Anchor& rhs) { return !(lhs < rhs); }

    boost::property_tree::ptree to_ptree() const
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
    }

private:
    Contig_Entry_CBPtr _ce_cbptr;
    Mutation_CBPtr _mut_cbptr;
    bool _c_right;

}; // class Allele_Anchor

} // namespace MAC


#endif
