#ifndef __ALLELE_ANCHOR_HPP
#define __ALLELE_ANCHOR_HPP

#include "Allele_Specifier.hpp"
#include "MAC_sets.hpp"

namespace MAC
{

/** An allele anchor is either a Mutation or a Contig_Entry edge */
class Allele_Anchor
{
public:
    DEFAULT_DEF_CTOR(Allele_Anchor);

    explicit Allele_Anchor(Mutation_CBPtr mut_cbptr, Contig_Entry_CBPtr ce_cbptr)
        : _ce_cbptr(ce_cbptr), _mut_cbptr(mut_cbptr), _c_right(true) {}

    Allele_Anchor(Contig_Entry_CBPtr ce_cbptr, bool c_right)
        : _ce_cbptr(ce_cbptr), _mut_cbptr(nullptr), _c_right(c_right) {}

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
    Anchor_Read_Support read_support(unsigned min_edge_support = 0) const;

    /// Distance between consecutive allele anchors.
    static Size_Type dist(const Allele_Anchor& lhs, const Allele_Anchor& rhs);

    friend bool operator <  (const Allele_Anchor& lhs, const Allele_Anchor& rhs);
    friend bool operator == (const Allele_Anchor& lhs, const Allele_Anchor& rhs);
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
