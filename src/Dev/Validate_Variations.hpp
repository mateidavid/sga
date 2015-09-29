#ifndef __VALIDATE_VARIATIONS_HPP
#define __VALIDATE_VARIATIONS_HPP

#include "Graph.hpp"
#include "BWT.h"
#include "BWTAlgorithms.h"


namespace MAC
{

class Validate_Variations
{
public:
    Validate_Variations(Graph& g) : Validate_Variations(g, 10, 2, 2, false) {}
    Validate_Variations(Graph& g,
                        size_t flank_size,
                        size_t min_graph_support_to_skip,
                        size_t min_read_support_to_validate,
                        bool check_both_strands)
        : _g(g),
          _flank_size(flank_size),
          _min_graph_support_to_skip(min_graph_support_to_skip),
          _min_read_support_to_validate(min_read_support_to_validate),
          _check_both_strands(check_both_strands)
    {}

    void operator () () const;

private:
    bool check_sequence(const Seq_Type& s) const;
    bool validate_mutation_allele(Mutation_CBPtr mut_cbptr, const set< Read_Chunk_CBPtr >& rc_set) const;
    bool validate_edge(Read_Chunk_CBPtr rc_cbptr) const;

    Graph& _g;
    size_t _flank_size;
    size_t _min_graph_support_to_skip;
    size_t _min_read_support_to_validate;
    bool _check_both_strands;
}; // class Validate_Variations

} // namespace MAC


#endif
