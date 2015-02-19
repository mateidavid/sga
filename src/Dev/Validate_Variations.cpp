#include "Validate_Variations.hpp"
#include "BWTAlgorithms.h"

namespace MAC
{

void Validate_Variations::operator () (Graph& g, const BWTIndexSet& index_set) const
{
    LOG("Validate_Variations", info) << ptree("Validate_Variations__start");
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
            // place chunks into bins:
            //   qr_set: observe mutation allele
            //   part_rf_set: observe reference allele but do not span fully
            //   full_rf_set: observe reference allele, spanning fully
            set< Read_Chunk_CBPtr > qr_set;
            set< Read_Chunk_CBPtr > part_rf_set;
            set< Read_Chunk_CBPtr > full_rf_set;
            for (auto mca_cbptr : mut_bptr->chunk_ptr_cont() | referenced)
            {
                qr_set.insert(mca_cbptr->chunk_cbptr());
            }
            auto iint_res = ce_bptr->chunk_cont().iintersect(mut_bptr->rf_start(), mut_bptr->rf_end());
            for (auto rc_cbptr : iint_res | referenced)
            {
                if (qr_set.count(rc_cbptr))
                {
                    continue;
                }
                if (rc_cbptr->get_c_start() <= mut_bptr->rf_start()
                    and mut_bptr->rf_end() <= rc_cbptr->get_c_end())
                {
                    full_rf_set.insert(rc_cbptr);
                }
                else
                {
                    part_rf_set.insert(rc_cbptr);
                }
            }
            Range_Type c_rg(mut_bptr->rf_start(), mut_bptr->rf_end());
            bool validated_qr = (qr_set.size() >= _min_graph_support_to_skip)
                or validate_allele(c_rg, qr_set, index_set);
            bool validated_rf = (full_rf_set.size() >= _min_graph_support_to_skip)
                or validate_allele(c_rg, full_rf_set, index_set);
            LOG("Validate_Variations", debug) << ptree("validation_result")
                .put("mut_ptr", mut_bptr.to_int())
                .put("validated_qr", validated_qr)
                .put("validated_rf", validated_rf);

            //TODO
        }
    }
    LOG("Validate_Variations", info) << ptree("Validate_Variations__end");
} // Validate_Variations::operator ()

bool Validate_Variations::validate_allele(const Range_Type& c_rg,
                                         const set< Read_Chunk_CBPtr >& rc_set,
                                         const BWTIndexSet& index_set) const
{
    if (rc_set.size() == 0)
    {
        LOG("Validate_Variations", debug) << ptree("validate_allele__empty_support");
        return false;
    }
    // find a read supporting the allele with a large enough flank
    Range_Type r_rg;
    auto rc_cbptr_it = find_if(rc_set.begin(), rc_set.end(), [&] (Read_Chunk_CBPtr rc_cbptr) {
            r_rg = rc_cbptr->mapped_range(c_rg, true, true, true);
            return (r_rg.start() >= _flank_size and
                    r_rg.end() + _flank_size <= rc_cbptr->re_bptr()->len());
        });
    if (rc_cbptr_it == rc_set.end())
    {
        LOG("Validate_Variations", debug) << ptree("validate_allele__no_flank_support");
        return false;
    }
    Read_Chunk_CBPtr rc_cbptr = *rc_cbptr_it;
    Seq_Type s = rc_cbptr->re_bptr()->get_seq().substr(r_rg.start() - _flank_size,
                                                       2 * _flank_size + r_rg.len());
    if (_check_both_strands)
    {
        auto cnt_0 = BWTAlgorithms::countSequenceOccurrencesSingleStrand(s, index_set);
        auto cnt_1 = BWTAlgorithms::countSequenceOccurrencesSingleStrand(s.revcomp(), index_set);
        return cnt_0 >= _min_read_support_to_validate and cnt_1 >= _min_read_support_to_validate;
    }
    else
    {
        auto cnt = BWTAlgorithms::countSequenceOccurrences(s, index_set);
        return cnt >= _min_read_support_to_validate;
    }
} // Validate_Variations::validate_allele

} // namespace MAC
