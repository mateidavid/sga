#ifndef __READ_MERGER_HPP
#define __READ_MERGER_HPP

#include "Graph.hpp"
#include "Anchor_Support.hpp"

namespace MAC
{

class Read_Merger
{
public:
    Read_Merger(Graph& g, unsigned max_discordant_support) : _g(g), _max_discordant_support(max_discordant_support) {}

    void operator () ();

private:
    Graph& _g;
    unsigned _max_discordant_support;

    void extend_haploid_support(
        const Allele_Anchor& anchor, const Allele_Specifier& allele, bool c_direction,
        RE_OSet& allele_support);
    void merge_reads(Contig_Entry_BPtr ce_bptr, const RE_OSet& common_support);

    pair< Read_Entry_CBPtr, Read_Entry_CBPtr >
    split_read(Read_Entry_CBPtr re_cbptr, const Allele_Anchor& l_anchor, const Allele_Anchor& r_anchor);

    Read_Chunk_BPtr merge_contig_chunks(const RC_OSet& crt_chunks, Read_Entry_BPtr m_re_bptr);
    Read_Chunk_BPtr merge_unmappable_chunks(const RC_OSet& crt_chunks, Read_Entry_BPtr m_re_bptr);
    pair< RC_OSet, RC_OSet > advance_chunks(const RC_OSet& crt_rc_oset, bool direction);
}; // class Read_Merger

} // namespace MAC

#endif
