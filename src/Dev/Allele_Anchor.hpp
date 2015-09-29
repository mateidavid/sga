#ifndef __ALLELE_ANCHOR_HPP
#define __ALLELE_ANCHOR_HPP

#include "Mutation.hpp"
#include "Contig_Entry.hpp"
#include "Allele_Specifier.hpp"
#include "OSet.hpp"

namespace MAC
{

typedef OSet< Read_Chunk_CBPtr > Allele_Chunk_Support;
typedef map< Allele_Specifier, Allele_Chunk_Support > Anchor_Chunk_Support;
typedef OSet< Read_Entry_CBPtr > Allele_Read_Support;
typedef map< Allele_Specifier, Allele_Read_Support > Anchor_Read_Support;

/** An allele anchor is either a Mutation or a Contig_Entry edge */
class Allele_Anchor
{
public:
    DEFAULT_DEF_CTOR(Allele_Anchor);

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
    typedef map< pair< Allele_Specifier, Allele_Specifier >, Allele_Read_Support > anchor_connect_type;
    static anchor_connect_type
    connect(const Anchor_Read_Support& a1_support_m,
            const Anchor_Read_Support& a2_support_m,
            bool same_st = true);

    /** Get sibling anchor.
     * PRE: Sibling anchor must exist in the given direction.
     * @param c_direction Bool; if false: get anchor to the right, else to the left
     * @return Sibling anchor.
     */
    Allele_Anchor get_sibling(bool c_direction) const;

    Anchor_Chunk_Support chunk_support(unsigned min_edge_support = 0) const;

    /// Distance between consecutive allele anchors.
    static Size_Type dist(const Allele_Anchor& lhs, const Allele_Anchor& rhs)
    {
        ASSERT(lhs.get_sibling(false) == rhs);
        return (rhs.is_endpoint()? rhs.ce_cbptr()->len() : rhs.mut_cbptr()->rf_start())
            - (lhs.is_endpoint()? 0 : lhs.mut_cbptr()->rf_end());
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

    ptree to_ptree() const;
    friend ostream& operator << (ostream& os, const Allele_Anchor& a);

private:
    Contig_Entry_CBPtr _ce_cbptr;
    Mutation_CBPtr _mut_cbptr;
    bool _c_right;
}; // class Allele_Anchor

} // namespace MAC


#endif
