#include "Validate_Variations.hpp"
#include "BWTAlgorithms.h"

namespace MAC
{

void Validate_Variations::operator () () const
{
    if (not _g.aux_index_set().pBWT)
    {
        LOG("Validate_Variations", error) << "auxiliary BWT index is required" << endl;
        exit(EXIT_FAILURE);
    }
    LOG("Validate_Variations", info) << ptree("begin");
    // search for well-separated mutations with low support
    for (auto ce_bptr : _g.ce_cont() | referenced)
    {
        size_t min_mut_start = 0;
        LOG("Validate_Variations", debug1) << ptree("loop").put("ce_bptr", ce_bptr.to_int());
        auto mut_cit = ce_bptr->mut_cont().cbegin();
        while (mut_cit != ce_bptr->mut_cont().cend())
        {
            auto mut_cbptr = &*mut_cit;
            ++mut_cit;
            // check no previous mutation overlaps current one
            if (mut_cbptr->rf_start() < min_mut_start)
            {
                LOG("Validate_Variations", debug) << ptree("validation_skip")
                    .put("mut_ptr", mut_cbptr.to_int());
                min_mut_start = max(min_mut_start, mut_cbptr->rf_end());
                continue;
            }
            min_mut_start = mut_cbptr->rf_end() + 1;
            // check next mutation doesn't overlap it either
            if (mut_cit != ce_bptr->mut_cont().end() and mut_cit->rf_start() < min_mut_start)
            {
                LOG("Validate_Variations", debug) << ptree("validation_skip")
                    .put("mut_ptr", mut_cbptr.to_int());
                min_mut_start = max(min_mut_start, mut_cit->rf_end());
                continue;
            }
            // found a separated mutation
            // place chunks into bins:
            //   qr_set: observe mutation allele
            //   part_rf_set: observe reference allele but do not span fully
            //   full_rf_set: observe reference allele, spanning fully
            set< Read_Chunk_CBPtr > qr_set;
            set< Read_Chunk_CBPtr > part_rf_set;
            set< Read_Chunk_CBPtr > full_rf_set;
            bool validated_qr;
            bool validated_rf;
            auto validate_allele_pair = [&] () {
                tie(qr_set, full_rf_set, part_rf_set) = ce_bptr->mut_support(mut_cbptr);
                validated_qr = ((qr_set.size() >= _min_graph_support_to_skip)
                                or validate_mutation_allele(mut_cbptr, qr_set));
                validated_rf = ((full_rf_set.size() >= _min_graph_support_to_skip)
                                or validate_mutation_allele(mut_cbptr, full_rf_set));
            };
            validate_allele_pair();

            LOG("Validate_Variations", info) << ptree("validation_result.mutation")
                .put("ce_ptr", ce_bptr.to_int())
                .put("mut_ptr", mut_cbptr.to_int())
                .put("validated_qr", validated_qr)
                .put("validated_rf", validated_rf);

            if (validated_qr and validated_rf)
            {
                // both validated; nothing to do
                continue;
            }
            if (not validated_qr and not validated_rf)
            {
                // neither validated; we don't to anything for now; perhaps unmap this region?
                LOG("Validate_Variations", info) << ptree("neither_allele_validated")
                    .put("mut_bptr", mut_cbptr.to_int());
                continue;
            }
            if (not validated_rf)
            {
                ASSERT(validated_qr);
                // we swap the mutation allele with the reference allele
                mut_cbptr = Contig_Entry::swap_mutation_alleles(ce_bptr, mut_cbptr).unconst();
                if (not mut_cbptr)
                {
                    LOG("Validate_Variations", info) << ptree("mutation_disappeared");
                    continue;
                }
                // re-run validation
                validate_allele_pair();
            }
            ASSERT(validated_rf);
            ASSERT(not validated_qr);
            // we edit the reads that contain the mutation allele
            set< Read_Entry_CBPtr > re_set;
            Range_Type c_rg(mut_cbptr->rf_start(), mut_cbptr->rf_end());
            for (auto mca_bptr : mut_cbptr->chunk_ptr_cont() | referenced)
            {
                auto rc_bptr = mca_bptr->chunk_cbptr().unconst();
                auto re_bptr = rc_bptr->re_bptr();
                auto r_rg = rc_bptr->mapped_range(c_rg, true, true, true);
                LOG("Validate_Variations", info) << ptree("erasing_nonvalidated_allele")
                    .put("re_bptr", re_bptr.to_int())
                    .put("r_rg_start", r_rg.begin())
                    .put("r_rg_end", r_rg.end());
                ASSERT(r_rg.len() == mut_cbptr->seq_len());
                re_bptr->edit(rc_bptr, r_rg.begin(),
                              mut_cbptr->seq().revcomp(rc_bptr->get_rc()),
                              ce_bptr->substr(c_rg.begin(), c_rg.len()));
                re_set.insert(re_bptr);
            }
            // remove the Mutation
            ce_bptr->mut_cont().erase(mut_cbptr);
            mut_cbptr.unconst()->chunk_ptr_cont().clear_and_dispose();
            Mutation_Fact::del_elem(mut_cbptr);
            _g.check(re_set);
        }
        _g.check(set< Contig_Entry_CBPtr >({ ce_bptr }));
    }
    _g.check_all();
    // search for inter-contig edges with low support
    for (auto ce_bptr : _g.ce_cont() | referenced)
    {
        if (ce_bptr->is_unmappable()) continue;
        for (int c_direction = 0; c_direction < 2; ++c_direction)
        {
            //auto oc = ce_bptr->out_chunks_dir(not c_direction, 3);
            auto oc = ce_bptr->out_chunks(c_direction, true);
            for (const auto& p : oc)
            {
                ASSERT(not p.second.empty());
                if (p.second.size() > 1) continue;
                Read_Chunk_CBPtr rc_cbptr = (not p.second[0].empty()? *p.second[0].begin() : *p.second[1].begin());
                int r_direction = (c_direction + rc_cbptr->get_rc()) % 2;
                if (r_direction) continue;
                bool res = validate_edge(rc_cbptr);
                LOG("Validate_Variations", info) << ptree("validation_result.edge")
                    .put("ce_bptr", ce_bptr.to_int())
                    .put("c_direction", c_direction)
                    .put("next_ce_bptr", p.first.ce_next_cbptr().to_int())
                    .put("same_orientation", p.first.same_orientation())
                    .put("validated", res);
                if (not res)
                {
                    _g.split_read(rc_cbptr);
                }
            }
        }
    }
    _g.check_all();
    LOG("Validate_Variations", info) << ptree("end");
} // Validate_Variations::operator ()

bool Validate_Variations::check_sequence(const Seq_Type& s) const
{
    if (_check_both_strands)
    {
        auto cnt_0 = BWTAlgorithms::countSequenceOccurrencesSingleStrand(s, _g.aux_index_set());
        auto cnt_1 = BWTAlgorithms::countSequenceOccurrencesSingleStrand(s.revcomp(), _g.aux_index_set());
        return cnt_0 >= _min_read_support_to_validate and cnt_1 >= _min_read_support_to_validate;
    }
    else
    {
        auto cnt = BWTAlgorithms::countSequenceOccurrences(s, _g.aux_index_set());
        return cnt >= _min_read_support_to_validate;
    }
} // Validate_Variations::check_sequence

bool
Validate_Variations::validate_mutation_allele(
    Mutation_CBPtr mut_cbptr, const set< Read_Chunk_CBPtr >& rc_set) const
{
    if (rc_set.size() == 0)
    {
        LOG("Validate_Variations", debug) << ptree("empty_support");
        return false;
    }
    // find a read supporting the allele with a large enough flank
    Range_Type c_rg(mut_cbptr->rf_start(), mut_cbptr->rf_end());
    Range_Type r_rg;
    auto rc_cbptr_it = find_if(rc_set.begin(), rc_set.end(), [&] (Read_Chunk_CBPtr rc_cbptr) {
            r_rg = rc_cbptr->mapped_range(c_rg, true, true, true); // TODO: CHECK
            return (r_rg.begin() >= rc_cbptr->re_bptr()->start() + _flank_size and
                    r_rg.end() + _flank_size <= rc_cbptr->re_bptr()->end());
        });
    if (rc_cbptr_it == rc_set.end())
    {
        LOG("Validate_Variations", debug) << ptree("no_flank_support");
        return false;
    }
    Read_Chunk_CBPtr rc_cbptr = *rc_cbptr_it;
    Seq_Type s = rc_cbptr->re_bptr()->get_seq().substr(
        r_rg.begin() - rc_cbptr->re_bptr()->start() - _flank_size,
        2 * _flank_size + r_rg.len());
    return check_sequence(s);
} // Validate_Variations::validate_mutation_allele

bool
Validate_Variations::validate_edge(Read_Chunk_CBPtr rc_cbptr) const
{
    LOG("Validate_Variations", debug) << ptree("begin")
        .put("re_bptr", rc_cbptr->re_bptr().to_int())
        .put("rc_bptr", rc_cbptr.to_int());
    Read_Entry_CBPtr re_cbptr = rc_cbptr->re_bptr();
    Read_Chunk_CBPtr next_rc_cbptr = re_cbptr->chunk_cont().get_sibling(rc_cbptr, true, true);
    bool have_unmappable_region = false;
    ASSERT(next_rc_cbptr);
    if (next_rc_cbptr->ce_bptr()->is_unmappable())
    {
        have_unmappable_region = true;
        next_rc_cbptr = re_cbptr->chunk_cont().get_sibling(next_rc_cbptr, true, true);
        ASSERT(next_rc_cbptr);
        ASSERT(not next_rc_cbptr->ce_bptr()->is_unmappable());
    }
    Seq_Type seq = re_cbptr->get_seq();
    Size_Type s1_start = (rc_cbptr->get_r_end() >= re_cbptr->start() + _flank_size
                          ? rc_cbptr->get_r_end() - _flank_size
                          : re_cbptr->start());
    Size_Type s1_end = (rc_cbptr->get_r_end() + _flank_size <= re_cbptr->end()
                        ? rc_cbptr->get_r_end() + _flank_size
                        : re_cbptr->end());
    Size_Type s2_start = (next_rc_cbptr->get_r_start() >= re_cbptr->start() + _flank_size
                          ? next_rc_cbptr->get_r_start() - _flank_size
                          : 0);
    Size_Type s2_end = (next_rc_cbptr->get_r_start() + _flank_size <= re_cbptr->end()
                        ? next_rc_cbptr->get_r_start() + _flank_size
                        : re_cbptr->end());
    bool res = (check_sequence(seq.substr(s1_start - re_cbptr->start(), s1_end - s1_start))
                and (not have_unmappable_region
                     or check_sequence(seq.substr(s2_start - re_cbptr->start(), s2_end - s2_start))));
    LOG("Validate_Variations", debug) << ptree("end").put("res", res);
    return res;
} // Validate_Variations::validate_edge

} // namespace MAC
