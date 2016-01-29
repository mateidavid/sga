#ifndef __VARIATION_CLEANER_HPP
#define __VARIATION_CLEANER_HPP

#include "Graph.hpp"


namespace MAC
{

class Variation_Cleaner
{
public:
    Variation_Cleaner(Graph& g, size_t min_num_3mers = 8)
        : _g(g), _min_num_3mers(min_num_3mers) {}

    void operator() () const;

    void unmap_clusters() const;
    void sync_clusters() const;

private:
    /**
     * Find clusters of non-separated mutations.
     * 2 mutations are not well separated if the non-mutated region between them
     * does not contain the given minimum number of 3mers.
     * @param ce_cbptr Contig_Entry to search in.
     * @param start_pos Minimum start position to search for clusters.
     * @return A pair of: (1) a set of non-separated mutations; and (2) 2 bools,
     * signaling if there is no spearation between first/last mutation and the
     * left/right endpoint (or between the endpoints themselves, if no mutations are present).
     */
    pair< set< Mutation_CBPtr >, array< bool, 2 > >
    find_mutation_cluster(Contig_Entry_CBPtr ce_cbptr, Size_Type start_pos) const;

    /** Check that the sequence given contains the minimum number of 3mers. */
    bool check_3mer_threshold(const Seq_Proxy_Type& seq) const;

    Graph& _g;
    size_t _min_num_3mers;
}; // class Variation_Cleaner

} // namespace MAC

#endif
