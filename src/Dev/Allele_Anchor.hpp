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

    Allele_Anchor(Mutation_CBPtr mut_cbptr, bool rf_allele)
        : _mut_cbptr(mut_cbptr), _rf_allele(rf_allele)
    {
        ASSERT(mut_cbptr);
        ASSERT(not mut_cbptr->chunk_ptr_cont().empty());
        _ce_cbptr = mut_cbptr->chunk_ptr_cont().begin()->chunk_cbptr()->ce_bptr();
    }

    Allele_Anchor(Contig_Entry_CBPtr ce_cbptr, Contig_Entry_CBPtr ce_next_cbptr,
                  bool c_right, bool same_orientation)
        : _ce_cbptr(ce_cbptr), _ce_next_cbptr(ce_next_cbptr),
          _c_right(c_right), _same_orientation(same_orientation)
    {
        ASSERT(ce_cbptr);
        ASSERT(ce_next_cbptr);
    }

    bool is_mutation_allele() const { return _mut_cbptr; }
    bool is_contig_edge() const { return _ce_next_cbptr; }

    GETTER(Mutation_CBPtr, mut_cbptr, _mut_cbptr)
    GETTER(bool, rf_allele, _rf_allele)
    GETTER(Contig_Entry_CBPtr, ce_cbptr, _ce_cbptr)
    GETTER(Contig_Entry_CBPtr, ce_next_cbptr, _ce_next_cbptr)
    GETTER(bool, c_right, _c_right)
    GETTER(bool, same_orientation, _same_orientation)

    /** Return the set of read chunks (in ce_cbptr()) supporting this anchor.
     */
    set< Read_Chunk_CBPtr > support() const
    {
        set< Read_Chunk_CBPtr > res;
        if (is_mutation_allele())
        {
            if (rf_allele())
            {
                // add all chunks fully spanning mutation
                // if the rf allele is empty, we require >=1bp mapped on either side
                for (rc_cbptr : ce_cbptr()->iintersect(mut_cbptr()->rf_start(), mut_cbptr()->rf_end()) | referenced)
                {
                    if ((mut_cbptr()->rf_len() > 0
                         and rc_cbptr->c_start() <= _mut_cbptr->rf_start()
                         and rc_cbptr->c_end() >= _mut_cbptr->rf_end())
                        or (mut_cbptr()->rf_len() == 0
                            and rc_cbptr->c_start() < _mut_cbptr->rf_start()
                            and rc_cbptr->c_end() > _mut_cbptr->rf_end()))
                    {
                        res.insert(rc_cbptr);
                    }
                }
                // then, remove the chunks observing the alternate allele
                auto qr_res = Allele_Anchor(mut_cbptr(), false).support();
                for (rc_cbptr : qr_res)
                {
                    // if mut is not an insertion, every chunk observing qr allele should already be in res
                    ASSERT(not (mut_cbptr()->rf_len() > 0) or res.count(rc_cbptr) == 1);
                    // if mut is an insertion, chunks observing qr allele should be in res iff they span >=1bp in both directions
                    ASSERT(not (mut_cbptr()->rf_len() == 0)
                           or ((res.count(rc_cbptr) == 1)
                               == (rc_cbptr->c_start() < _mut_cbptr->rf_start()
                                   and rc_cbptr->c_end() > _mut_cbptr->rf_end())));
                    res.erase(rc_cbptr);
                }
            }
            else // is_mutation_allele && qr_allele
            {
                for (auto mca_cbptr : mut_cbptr()->chunk_ptr_cont() | referenced)
                {
                    res.insert(mca_cbptr->chunk_cbptr());
                }
            }
        }
        else // is_contig_edge
        {
            auto m = ce_cbptr()->out_chunks_dir(c_right(), 3);
            auto k = make_tuple(ce_next_cbptr(), same_orientation());
            ASSERT(m.count(k) == 1);
            res = move(m[k]);
        }
        return res;
    }

    /** Return (indirectly) the set of Read_Entry supporting 2 anchors.
     * For every Read_Entry supporting both anchors, return a pait of chunks,
     * the first supporting the first anchor, the second supporting the second anchor.
     * @param a1 First anchor.
     * @param a1 First anchor.
     * @param same_st Bool; if true, consider reads supporting the anchors on the same strand;
     * if false, on different strands.
     */
    static set< pair< Read_Chunk_CBPtr, Read_Chunk_CBPtr > >
    connect(const Allele_Anchor& a1, const Allele_Anchor& a2, bool same_st = true)
    {
        set< pair< Read_Chunk_CBPtr, Read_Chunk_CBPtr > > res;
        // for each anchor, compute support
        auto res1 = a1.support();
        auto res2 = a2.support();
        // for each anchor, construct maps of the form [Read_Entry, strand] -> Read_Chunk
        auto re_orientation_rg_1 = ba::transform(res1, [] (Read_Chunk_CBPtr rc_cbptr) {
                return make_pair(make_pair(rc_cbptr->re_bptr(), rc_cbptr->get_rc()), rc_cbptr);
            });
        auto re_orientation_rg_2 = ba::transform(res2, [] (Read_Chunk_CBPtr rc_cbptr) {
                return make_pair(make_pair(rc_cbptr->re_bptr(), rc_cbptr->get_rc()), rc_cbptr);
            });
        map< pair< Read_Entry_CBPtr, bool >, Read_Chunk_CBPtr >
            re_orientation_map_1(re_orientation_rg_1.begin(), re_orientation_rg_1.end());
        map< pair< Read_Entry_CBPtr, bool >, Read_Chunk_CBPtr >
            re_orientation_map_2(re_orientation_rg_2.begin(), re_orientation_rg_2.end());
        // iterate through map of a1 support
        for (auto p1 : re_orientation_map_1 | ba::map_keys)
        {
            Read_Entry_CBPtr re_cbptr = p1.first;
            bool a1_neg_strand = p1.second;
            bool a2_neg_strand = (a1_neg_strand != same_st);
            auto p2 = make_pair(re_cbptr, a2_neg_strand);
            if (re_orientation_map_2.count(p2) == 0)
            {
                continue;
            }
            // found pair of chunks of the same read supporting each anchor
            res.insert(make_pair(re_orientation_map_1[p1], re_orientation_map_2[p2]));
        }
        return res;
    }

private:
    Mutation_CBPtr _mut_cbptr;
    bool _rf_allele;

    Contig_Entry_CBPtr _ce_cbptr;
    Contig_Entry_CBPtr _ce_next_cbptr;
    bool _c_right;
    bool _same_orientation;

}; // class Allele_Anchor

} // namespace MAC


#endif
