#ifndef __READ_MERGER_HPP
#define __READ_MERGER_HPP

#include "Graph.hpp"
//#include "Anchor_Support.hpp"

namespace MAC
{

class Read_Merger
{
public:
    Read_Merger(Graph& g, unsigned max_discordant_support, unsigned merged_weight)
        : _g(g), _max_discordant_support(max_discordant_support), _merged_weight(merged_weight) {}

    void operator () () const;

private:
    Graph& _g;
    unsigned _max_discordant_support;
    unsigned _merged_weight;

    void extend_haploid_support(
        const Allele_Anchor& anchor, const Allele_Specifier& allele, bool c_direction,
        RE_DSet& allele_support) const;
    void merge_reads(Contig_Entry_BPtr ce_bptr, const RE_DSet& common_support) const;

    pair< Read_Entry_CBPtr, Read_Entry_CBPtr >
    split_read(Read_Entry_CBPtr re_cbptr, const Allele_Anchor& l_anchor, const Allele_Anchor& r_anchor) const;

    Read_Chunk_BPtr merge_contig_chunks(const RC_DSet& crt_chunks, Read_Entry_BPtr m_re_bptr) const;
    Read_Chunk_BPtr merge_unmappable_chunks(const RC_DSet& crt_chunks, Read_Entry_BPtr m_re_bptr) const;
    pair< RC_DSet, RC_DSet > advance_chunks(const RC_DSet& crt_rc_dset, bool direction) const;
}; // class Read_Merger

} // namespace MAC

#endif
