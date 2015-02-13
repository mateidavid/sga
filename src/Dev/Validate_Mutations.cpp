#include "Validate_Mutations.hpp"

namespace MAC
{

void Validate_Mutations::operator () (Graph& g, const BWTIndexSet& index_set) const
{
    LOG("Validate_Mutations", info) << ptree("Validate_Mutations__start");
    // search for well-separated mutations with low support
    for (auto ce_bptr : g.ce_cont() | referenced)
    {
        size_t min_mut_start = 0;
        for (auto mut_bptr : ce_bptr->mut_cont() | referenced)
        {
            // check no previous mutation overlaps current one
            if (mut_bptr->rf_start() < min_mut_start)
            {
                min_mut_start = max(min_mut_start, mut_bptr->rf_end());
                continue;
            }
            min_mut_start = mut_bptr->rf_end();
            // check next mutation doesn't overlap it either
            {
                auto next_mut_it = next(ce_bptr->mut_cont().iterator_to(*mut_bptr));
                if (next_mut_it != ce_bptr->mut_cont().end() and next_mut_it->rf_start() < min_mut_start)
                {
                    min_mut_start = max(min_mut_start, next_mut_it->rf_end());
                    continue;
                }
            }
            // found a separated mutation
            set< Read_Chunk_CBPtr > rf_set;
            set< Read_Chunk_CBPtr > qr_set;
            for (auto mca_cbptr : mut_bptr->chunk_ptr_cont() | referenced)
            {
                qr_set.insert(mca_cbptr->chunk_cbptr());
            }
            auto iint_res = ce_bptr->chunk_cont().iintersect(mut_bptr->rf_start(), mut_bptr->rf_end());
            for (auto rc_cbptr : iint_res | referenced)
            {
                // if chunk observes qr allele, skip it
                if (qr_set.count(rc_cbptr)) continue;
                rf_set.insert(rc_cbptr);
            }
            //TODO
            //IPR
        }
    }
    LOG("Validate_Mutations", info) << ptree("Validate_Mutations__end");
}

} // namespace MAC
