#ifndef __UNMAPPER_HPP
#define __UNMAPPER_HPP

#include "Contig_Entry.hpp"
#include "Read_Entry.hpp"

namespace MAC
{

class Unmapper
{
public:
    typedef set< Contig_Entry_BPtr > ce_set_type;
    typedef map< Read_Entry_BPtr, Range_Cont > re_set_type;

    explicit Unmapper(Graph* g_p) : _g_p(g_p) {}

    void unmap(Read_Chunk_BPtr rc_bptr) const
    {
        ce_set_type ce_set;
        re_set_type re_set;
        if (rc_bptr->get_c_start() == 0 and rc_bptr->get_c_end() == rc_bptr->ce_bptr()->len())
        {
            ce_set.insert(rc_bptr->ce_bptr());
        }
        else
        {
            Range_Type rg(rc_bptr->get_r_start(), rc_bptr->get_r_end());
            re_set[rc_bptr->re_bptr()].insert(rg);
        }
        _unmap_loop(move(ce_set), move(re_set));
    }
    void unmap(Read_Entry_BPtr re_bptr, Range_Cont&& rg_cont) const
    {
        re_set_type re_set;
        re_set[re_bptr] = move(rg_cont);
        _unmap_loop(ce_set_type(), move(re_set));
    }
    void unmap(ce_set_type&& ce_set) const { _unmap_loop(move(ce_set), re_set_type()); }
    void unmap(re_set_type&& re_set) const { _unmap_loop(ce_set_type(), move(re_set)); }

    /**
     * Unmap read chunks not supported by any other reads.
     */
    void unmap_single_chunks() const;

    /**
     * Unmap terminal read regions not supported by any other reads.
     */
    void unmap_single_terminal_regions() const;

    /**
     * Repeatedly unmap short contigs with length < min_len and out-degree > max_deg.
     */
    void unmap_short_contigs(unsigned min_len, unsigned max_deg) const;

    /** Unmap indels adjacent to homopolymers. */
    void unmap_homopolymer_indels(unsigned min_len = 6) const;

private:
    Graph* _g_p;

    /**
     * Main unmapper loop.
     * Repeatedly unmap contigs, read regions, or extend unmapped regions.
     */
    void _unmap_loop(ce_set_type&& unmap_ce_set, re_set_type&& unmap_re_set) const;

    /**
     * Unmap contig.
     * All unmapped chunks are added as unmapped regions to extend.
     */
    void _unmap_ce(Contig_Entry_BPtr ce_bptr, re_set_type& extend_re_set) const;

    /**
     * Unmap read region.
     * If there is any chunk not yet unmapped in this region,
     * the chunk is cut to the extent of the region,
     * its contig is added to the unmap_ce_set,
     * and the remaining read region is added back to the unmap_re_set.
     */
    void _unmap_re_region(Read_Entry_BPtr re_bptr, Range_Type rg,
                          ce_set_type& unmap_ce_set, re_set_type& unmap_re_set) const;

    /**
     * Extend unmapped read region.
     * If the chunk next to the unmappable region is itself unmappable, it is merged in.
     * If the unmappable region is ended by a mappable chunk smaller than unmap_trigger_len,
     * that chunk is added to the unmap_ce_set.
     */
    void _extend_unmappable_re_region(Read_Entry_BPtr re_bptr, Range_Type rg,
                                      re_set_type& unmap_re_set) const;
}; // class Unmapper

} // namespace MAC

#endif
