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
        const Allele_Anchor& _anchor, const Allele_Specifier& _allele, bool _c_direction,
        Allele_Read_Support& allele_support);
    void merge_reads(const Allele_Read_Support& allele_support);

    pair< Read_Entry_CBPtr, Read_Entry_CBPtr >
    split_read(Read_Entry_CBPtr re_cbptr, const Allele_Anchor& l_anchor, const Allele_Anchor& r_anchor);

}; // class Read_Merger

} // namespace MAC

#endif
