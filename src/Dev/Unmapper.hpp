#ifndef __UNMAPPER_HPP
#define __UNMAPPER_HPP

#include "Graph.hpp"

namespace MAC
{

class Unmapper
{
public:
    Unmapper(Graph& g) : _g(g);

    void unmap_chunk(Read_Chunk_BPtr rc_bptr)
    {
        rc_bptr = g.trim_contig_to_chunk(rc_bptr);
        set< Contig_Entry_BPtr > unmap_ce_set({rc_bptr->ce_bptr});
        set< pair< Read_Entry_BPtr, Range_Type > > unmap_re_set;
        _unmap_loop(unmap_ce_set, unmap_re_set);
    }
    void unmap_regions(set< pair< Read_Entry_BPtr, Range_Type > >& unmap_re_set)
    {
        set< Contig_Entry_BPtr > unmap_ce_set;
        _unmap_loop(unmap_ce_set, unmap_re_set);
    }

private:
    Graph& _g;

    /**
     * Main unmapper loop.
     * Repeatedly unmap contigs, read regions, or extend unmapped regions.
     */
    void _unmap_loop(set< Contig_Entry_BPtr >& unmap_ce_set,
                    set< pair< Read_Entry_BPtr, Range_Type > >& unmap_re_set);

    /**
     * Unmap contig.
     * All unmapped chunks are added as unmapped regions to extend.
     */
    void _unmap_ce(Contig_Entry_BPtr ce_bptr,
                   set< pair< Read_Entry_BPtr, Range_Type > >& extend_re_set);

    /**
     * Unmap read region.
     * If there is any chunk not yet unmapped in this region,
     * the chunk is cut to the extent of the region,
     * its contig is added to the unmap_ce_set,
     * and the remaining read region is added back to the unmap_re_set.
     */
    void _unmap_re_region(const pair< Read_Entry_BPtr, Range_Type >& p,
                          set< Contig_Entry_BPtr >& unmap_ce_set,
                          set< pair< Read_Entry_BPtr, Range_Type > > unmap_re_set);

    /**
     * Extend unmapped read region.
     * If the chunk next to the unmappable region is itself unmappable, it is merged in.
     * If the unmappable region is ended by a mappable chunk smaller than unmap_trigger_len,
     * that chunk is added to the unmap_ce_set.
     */
    void _extend_unmappable_re_region(const pair< Read_Entry_BPtr, Range_Type >& p,
                                      set< Contig_Entry_BPtr >& unmap_ce_set);

}; // class Unmapper

} // namespace MAC

#endif
