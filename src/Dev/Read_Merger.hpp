#ifndef __READ_MERGER_HPP
#define __READ_MERGER_HPP

#include "Graph.hpp"

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

    struct Traversal_Struct
    {
        Allele_Anchor anchor;
        Allele_Specifier allele;
        bool c_direction;
        Anchor_Chunk_Support chunk_support;
        RC_DSet unmappable_chunks;
        ptree to_ptree() const;
    };
    typedef list< Traversal_Struct > Traversal_List;

    bool extend_haploid_layout(Traversal_List& l, Traversal_List::iterator it) const;
    bool extend_haploid_layout_dir(Traversal_List& l, Traversal_List::iterator it, bool e_direction) const;
    void split_diverging_reads(Traversal_List& l, Traversal_List::iterator it) const;
    void merge_reads(Traversal_List& l, Traversal_List::iterator it) const;

    Read_Chunk_BPtr merge_contig_chunks(const RC_DSet& crt_chunks, Read_Entry_BPtr m_re_bptr, bool c_direction) const;
    Read_Chunk_BPtr merge_unmappable_chunks(const RC_DSet& crt_chunks, Read_Entry_BPtr m_re_bptr) const;
    pair< RC_DSet, RC_DSet > advance_chunks(const RC_DSet& crt_rc_dset, bool direction) const;
    void filter_duplicate_reads(RC_DSet& rc_dset, set< Read_Entry_CBPtr >& duplicate_re_set) const;
}; // class Read_Merger

} // namespace MAC

#endif
