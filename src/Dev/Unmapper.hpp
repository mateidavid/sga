#ifndef __UNMAPPER_HPP
#define __UNMAPPER_HPP

#include "Graph.hpp"

namespace MAC
{

class Unmapper
{
public:
    typedef set< Contig_Entry_BPtr > ce_set_type;
    typedef map< Read_Entry_BPtr, Range_Cont > re_set_type;

    Unmapper(Graph& g) : _g(g) {}

    void unmap_chunk(Read_Chunk_BPtr rc_bptr)
    {
        ce_set_type unmap_ce_set;
        re_set_type unmap_re_set;
        if (rc_bptr->get_c_start() == 0 and rc_bptr->get_c_end() == rc_bptr->ce_bptr()->len())
        {
            unmap_ce_set.insert(rc_bptr->ce_bptr());
        }
        else
        {
            Range_Type rg(rc_bptr->get_r_start(), rc_bptr->get_r_end());
            unmap_re_set[rc_bptr->re_bptr()].insert(rg);
        }
        _unmap_loop(move(unmap_ce_set), move(unmap_re_set));
    }
    void unmap_re_regions(re_set_type&& unmap_re_set)
    {
        ce_set_type unmap_ce_set;
        _unmap_loop(move(unmap_ce_set), move(unmap_re_set));
    }
    void unmap_re_regions(Read_Entry_BPtr re_bptr, Range_Cont&& rg_cont)
    {
        ce_set_type unmap_ce_set;
        re_set_type unmap_re_set;
        unmap_re_set[re_bptr] = move(rg_cont);
        _unmap_loop(move(unmap_ce_set), move(unmap_re_set));
    }

private:
    Graph& _g;

    /**
     * Main unmapper loop.
     * Repeatedly unmap contigs, read regions, or extend unmapped regions.
     */
    void _unmap_loop(ce_set_type&& unmap_ce_set, re_set_type&& unmap_re_set);

    /**
     * Unmap contig.
     * All unmapped chunks are added as unmapped regions to extend.
     */
    void _unmap_ce(Contig_Entry_BPtr ce_bptr, re_set_type& extend_re_set);

    /**
     * Unmap read region.
     * If there is any chunk not yet unmapped in this region,
     * the chunk is cut to the extent of the region,
     * its contig is added to the unmap_ce_set,
     * and the remaining read region is added back to the unmap_re_set.
     */
    void _unmap_re_region(Read_Entry_BPtr re_bptr, const Range_Type& rg,
                          ce_set_type& unmap_ce_set, re_set_type& unmap_re_set);

    /**
     * Extend unmapped read region.
     * If the chunk next to the unmappable region is itself unmappable, it is merged in.
     * If the unmappable region is ended by a mappable chunk smaller than unmap_trigger_len,
     * that chunk is added to the unmap_ce_set.
     */
    void _extend_unmappable_re_region(Read_Entry_BPtr re_bptr, const Range_Type& rg,
                                      re_set_type& unmap_re_set);

}; // class Unmapper

} // namespace MAC

#endif
