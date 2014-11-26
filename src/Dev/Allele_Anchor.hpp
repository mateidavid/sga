#ifndef __ALLELE_ANCHOR_HPP
#define __ALLELE_ANCHOR_HPP

#include "Mutation.hpp"
#include "Contig_Entry.hpp"

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

    /** Allele specifier type.
     * When the anchor is a mutation, an allele is specified by a bool:
     * rf=false; qr=true.
     * When the anchor is an endpoint, an allele is specified by the
     * destination of the edge: (ce_next, same_orientation)
     */
    typedef boost::variant< bool, pair< Contig_Entry_CBPtr, bool > > allele_specifier_type;

    bool is_mutation() const { return _mut_cbptr; }
    bool is_endpoint() const { return not is_mutation(); }

    GETTER(Contig_Entry_CBPtr, ce_cbptr, _ce_cbptr)
    GETTER(Mutation_CBPtr, mut_cbptr, _mut_cbptr)
    GETTER(bool, c_right, _c_right)

    /** Group read chunks supporting various alleles for this anchor.
     */
    typedef map< allele_specifier_type, set< Read_Chunk_CBPtr > > allele_support_type;
    allele_support_type support() const
    {
        allele_support_type res;
        if (is_mutation())
        {
            // first, get support for qr allele
            for (auto mca_cbptr : mut_cbptr()->chunk_ptr_cont() | referenced)
            {
                res[allele_specifier_type(true)].insert(mca_cbptr->chunk_cbptr());
            }
            // next, get support for rf allele
            // add all chunks fully spanning mutation
            for (auto rc_cbptr : ce_cbptr()->chunk_cont().iintersect(mut_cbptr()->rf_start(), mut_cbptr()->rf_end()) | referenced)
            {
                // if chunk observes qr allele, skip it
                if (res[allele_specifier_type(true)].count(rc_cbptr) > 0) { continue; }
                // if the rf allele is empty, we require >=1bp mapped on either side
                if ((mut_cbptr()->rf_len() > 0
                     and rc_cbptr->get_c_start() <= _mut_cbptr->rf_start()
                     and rc_cbptr->get_c_end() >= _mut_cbptr->rf_end())
                    or (mut_cbptr()->rf_len() == 0
                        and rc_cbptr->get_c_start() < _mut_cbptr->rf_start()
                        and rc_cbptr->get_c_end() > _mut_cbptr->rf_end()))
                {
                    res[allele_specifier_type(false)].insert(rc_cbptr);
                }
            }
        }
        else // is_endpoint
        {
            auto m = ce_cbptr()->out_chunks_dir(c_right(), 3);
            for (auto p : m)
            {
                res[allele_specifier_type(p.first)] = move(p.second);
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
    typedef map< pair< allele_specifier_type, allele_specifier_type >, set< Read_Entry_CBPtr > > anchor_connect_type;
    static anchor_connect_type
    connect(const allele_support_type& a1_support_m, const allele_support_type& a2_support_m, bool same_st = true)
    {
        anchor_connect_type res;
        // for each anchor, construct a map of the form
        //   [ allele_specifier_type ] -> set< [ Read_Entry, strand ] >
        auto make_allele_spec_re_set_map = [] (const allele_support_type& m) {
            map< allele_specifier_type, set< pair< Read_Entry_CBPtr, bool > > > r;
            for (const auto& p : m)
            {
                allele_specifier_type allele_spec = p.first;
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

private:
    Contig_Entry_CBPtr _ce_cbptr;
    Mutation_CBPtr _mut_cbptr;
    bool _c_right;

}; // class Allele_Anchor

} // namespace MAC


#endif
