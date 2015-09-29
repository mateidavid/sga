#ifndef __UNMAP_MUT_CLUSTERS_HPP
#define __UNMAP_MUT_CLUSTERS_HPP

#include "Graph.hpp"


namespace MAC
{

class Unmap_Mut_Clusters
{
public:
    Unmap_Mut_Clusters(Graph& g, size_t min_num_3mers = 8)
        : _g(g), _min_num_3mers(min_num_3mers) {}

    void operator() () const;

private:
    /** Search for 2 Mutations separated by less than the minimum number of 3mers.
     * @param ce_cbptr Contig_Entry to search in.
     * @return A single tuple (c_start,c_end) of such a range, if one exists.
     */
    optional< Range_Type >
    find_mutation_cluster(Contig_Entry_CBPtr ce_cbptr) const;

    /** Check that the sequence given contains the minimum number of 3mers. */
    bool check_3mer_threshold(const Seq_Proxy_Type& seq) const;

    Graph& _g;
    size_t _min_num_3mers;
};

} // namespace MAC


#endif
