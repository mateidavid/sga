#include "Read_Merger.hpp"

#include <seqan/store.h>
#include <seqan/consensus.h>

#include "filter_cont.hpp"

namespace MAC
{

void replace_chunk(Read_Merger::Traversal_List& l, Read_Merger::Traversal_List::iterator it, bool e_direction,
                   Read_Chunk_CBPtr rc_cbptr, Read_Chunk_CBPtr new_rc_cbptr, bool d)
{
    LOG("Read_Merger", debug1) << ptree("begin")
        .put("l", l)
        .put("rc_bptr", rc_cbptr.to_int())
        .put("new_rc_bptr", new_rc_cbptr.to_int());
    while (it->chunk_support.at(it->allele)[d].count(rc_cbptr))
    {
        it->chunk_support.at(it->allele)[d].erase(rc_cbptr);
        it->chunk_support.at(it->allele)[d].insert(new_rc_cbptr);
        if (not e_direction? it == prev(l.end()) : it == l.begin())
        {
            break;
        }
        it = not e_direction? next(it) : prev(it);
    }
    LOG("Read_Merger", debug1) << ptree("end")
        .put("l", l);
} // replace_chunk

void filter_duplicate_reads(RC_DSet& rc_dset, set< Read_Entry_CBPtr >& duplicate_re_set)
{
    map< Read_Entry_CBPtr, unsigned > re_count;
    for (int d = 0; d < 2; ++d)
    {
        for (auto rc_cbptr : rc_dset[d])
        {
            auto re_cbptr = rc_cbptr->re_bptr();
            if (not re_count.count(re_cbptr))
            {
                re_count[re_cbptr] = 0;
            }
            ++re_count.at(re_cbptr);
        }
    }
    for (auto& p : re_count)
    {
        if (p.second == 1) continue;
        duplicate_re_set.insert(p.first);
    }
    rc_dset.filter([&] (Read_Chunk_CBPtr rc_cbptr, bool) { return re_count.at(rc_cbptr->re_bptr()) == 1; });
} // filter_duplicate_reads

pair< RC_DSet, RC_DSet > advance_chunks(const RC_DSet& crt_rc_dset, bool c_direction)
{
    RC_DSet next_rc_dset;
    RC_DSet unmappable_rc_dset;
    for (int d = 0; d < 2; ++d)
    {
        for (auto rc_cbptr : crt_rc_dset[d])
        {
            ASSERT(rc_cbptr);
            bool r_direction = (c_direction + rc_cbptr->get_rc()) % 2;
            auto next_rc_cbptr = rc_cbptr->re_bptr()->chunk_cont().get_sibling(
                rc_cbptr, true, not r_direction);
            if (not next_rc_cbptr) continue;
            if (next_rc_cbptr->ce_bptr()->is_unmappable())
            {
                // skip one unmappable chunk if necessary
                unmappable_rc_dset[d].insert(next_rc_cbptr);
                auto next_next_rc_cbptr = next_rc_cbptr->re_bptr()->chunk_cont().get_sibling(
                    next_rc_cbptr, true, not r_direction);
                LOG("Read_Merger", debug1) << ptree()
                    .put("rc_bptr", rc_cbptr)
                    .put("next_rc_bptr", next_rc_cbptr)
                    .put("next_next_rc_bptr", next_next_rc_cbptr);
                if (not next_next_rc_cbptr)
                {
                    LOG("Read_Merger", warning) << ptree("terminal_unmappable_chunk")
                        .put("re_name", rc_cbptr->re_bptr()->name());
                    continue;
                }
                next_rc_cbptr = next_next_rc_cbptr;
            }
            else
            {
                LOG("Read_Merger", debug1) << ptree()
                    .put("rc_bptr", rc_cbptr)
                    .put("next_rc_bptr", next_rc_cbptr);
            }
            next_rc_dset[d].insert(next_rc_cbptr);
        }
    }
    return make_pair(next_rc_dset, unmappable_rc_dset);
} // advance_chunks

bool read_is_contained(Read_Entry_CBPtr re_cbptr)
{
    LOG("Read_Merger", debug) << ptree("begin")
        .put("re_bptr", re_cbptr.to_int())
        .put("re_name", re_cbptr->name());
    // get first chunk from this RE
    auto rc_cbptr = &*re_cbptr->chunk_cont().begin();
    auto ce_cbptr = rc_cbptr->ce_bptr();
    bool c_direction = rc_cbptr->get_rc();
    // get all chunks intersecting it
    auto rg = ce_cbptr->chunk_cont().iintersect(
        rc_cbptr->get_c_start(), rc_cbptr->get_c_end()) | referenced;
    RC_DSet rc_dset(rg, [] (Read_Chunk_CBPtr other_rc_cbptr) { return other_rc_cbptr->get_rc(); });
    ASSERT(rc_dset[rc_cbptr->get_rc()].count(rc_cbptr));
    rc_dset[rc_cbptr->get_rc()].erase(rc_cbptr);
    while (true)
    {
        LOG("Read_Merger", debug) << ptree("loop")
            .put("rc_bptr", rc_cbptr.to_int())
            .put("ce_bptr", ce_cbptr.to_int())
            .put("c_direction", c_direction);
        // filter rc_set, keeping only those chunks that extend rc_cbptr
        rc_dset.filter([&] (Read_Chunk_CBPtr other_rc_cbptr, bool) {
                return Read_Chunk::is_contained_in(rc_cbptr, other_rc_cbptr);
            });
        if (rc_dset.empty())
        {
            // no chunks fully contain rc_cbptr
            LOG("Read_Merger", debug) << ptree("end.no.no_chunks_contain_current_chunk");
            return false;
        }
        // find (contig,orientation) that follows rc_cbptr
        auto oc = ce_cbptr->out_chunks(c_direction, true);
        Allele_Specifier allele(nullptr, true);
        for (auto& p : oc)
        {
            if (p.second[rc_cbptr->get_rc()].count(rc_cbptr))
            {
                allele = p.first;
                break;
            }
        }
        if (not allele.ce_next_cbptr())
        {
            // rc_cbptr is the last chunk in the read
            ASSERT(rc_cbptr == &*re_cbptr->chunk_cont().rbegin());
            LOG("Read_Merger", debug) << ptree("end.yes.read_end")
                .put("other_re_bptr", (*rc_dset[rc_dset[0].empty()].begin())->re_bptr().to_int())
                .put("other_re_name", (*rc_dset[rc_dset[0].empty()].begin())->re_bptr()->name());
            return true;
        }
        // keep only chunks in rc_set that follow the same edge
        rc_dset.filter([&] (Read_Chunk_CBPtr other_rc_cbptr, bool) {
                return oc.at(allele)[other_rc_cbptr->get_rc()].count(other_rc_cbptr) > 0;
            });
        if (rc_dset.empty())
        {
            // no chunks continue to the same contig
            LOG("Read_Merger", debug) << ptree("end.no.no_chunks_to_next_contig");
            return false;
        }
        // advance rc_cbptr
        rc_cbptr = re_cbptr->chunk_cont().get_sibling(rc_cbptr, true, true);
        ASSERT(rc_cbptr);
        if (rc_cbptr->ce_bptr()->is_unmappable())
        {
            rc_cbptr = re_cbptr->chunk_cont().get_sibling(rc_cbptr, true, true);
            ASSERT(rc_cbptr);
        }
        // advance rc_dset
        rc_dset = advance_chunks(rc_dset, c_direction).first;
        rc_dset.reverse(not allele.same_orientation());
        // update vars
        ce_cbptr = rc_cbptr->ce_bptr();
        ASSERT(ce_cbptr == allele.ce_next_cbptr());
        c_direction = (c_direction + 1 + allele.same_orientation()) % 2;
    }
} // read_is_contained

void Read_Merger::operator () () const
{
    merge_haploid_alleles();
    remove_contained();
}

void Read_Merger::merge_haploid_alleles () const
{
    LOG("Read_Merger", info) << ptree("begin")
        .put("max_discordant_support", _max_discordant_support)
        .put("merged_weight", _merged_weight);
    //
    // Look for mutation alleles with copy number 1. When one is found:
    // - extract the reads that support it
    // - extend the reads past successive anchors as long as they are consistent
    //   (ie. support the same alleles)
    // - perform read splitting if necessary
    // - merge the corresponding reads
    // Restart the loop after each merge.
    //
    for (bool done = false; not done; )
    {
        done = true;
        for (auto ce_bptr : _g.ce_cont() | referenced)
        {
            if (not ce_bptr->is_normal()) continue;
            for (auto mut_bptr : ce_bptr->mut_cont() | referenced)
            {
                if (mut_bptr->copy_num(0) != 1 and mut_bptr->copy_num(1) != 1) continue;
                Allele_Anchor anchor(mut_bptr, ce_bptr);
                auto anchor_support = anchor.chunk_support();
                for (int al = 0; al < 2; ++al)
                {
                    if (mut_bptr->copy_num(al) != 1) continue;
                    Allele_Specifier allele(al);
                    if (anchor_support.at(allele).size() <= 1) continue;
                    //
                    // found a mutation allele with copy number 1
                    //
                    LOG("Read_Merger", debug) << ptree("haploid_mutation")
                        .put("anchor", anchor)
                        .put("allele", allele)
                        .put("allele_support", anchor_support.at(allele));
                    // initialize layout structs
                    Traversal_List l;
                    l.emplace_back(Traversal_Struct{anchor, allele, false,
                                Anchor_Chunk_Support(), RC_DSet()});
                    auto merge_start_it = l.begin();
                    merge_start_it->chunk_support.insert(make_pair(allele, move(anchor_support.at(allele))));
                    // remove any chunk belonging to REs that visit allele more than once
                    LOG("Read_Merger", debug) << ptree("before_filter_duplicate_reads").put("l", l);
                    auto size_before = merge_start_it->chunk_support.at(allele).size();
                    set< Read_Entry_CBPtr > duplicate_re_set;
                    filter_duplicate_reads(merge_start_it->chunk_support.at(allele), duplicate_re_set);
                    auto size_after = merge_start_it->chunk_support.at(allele).size();
                    ASSERT(size_after <= size_before);
                    if (size_after < size_before)
                    {
                        LOG("Read_Merger", warning) << ptree("duplicate_reads_at_haploid_allele")
                            .put("anchor", anchor)
                            .put("allele", allele)
                            .put("duplicate_re_set", duplicate_re_set);
                    }
                    if (size_after <= 1)
                    {
                        LOG("Read_Merger", debug) << ptree("nothing_to_extend_after_filter").put("l", l);
                        continue;
                    }
                    // extend layout
                    bool res = extend_haploid_layout(l, merge_start_it);
                    LOG("Read_Merger", debug) << ptree("haploid_extension").put("res", res);
                    if (not res) continue;
                    // extension was successful: split diverging reads
                    LOG("Read_Merger", debug) << ptree("before_split").put("l", l);
                    split_diverging_reads(l, merge_start_it);
                    // finally, merge reads
                    LOG("Read_Merger", debug) << ptree("before_merge").put("l", l);
                    merge_reads(l, merge_start_it);
                }
            }
        }
    }
    LOG("Read_Merger", info) << ptree("end");
} // Read_Merger::merge_haploid_alleles

void Read_Merger::remove_contained() const
{
    LOG("Read_Merger", info) << ptree("begin");
    for_each_it_advance(
        _g.re_cont(),
        [&] (Read_Entry_Cont::iterator re_it)
        {
            auto re_cbptr = &*re_it;
            if (read_is_contained(re_cbptr))
            {
                LOG("Read_Merger", debug) << ptree("loop.contained_read")
                    .put("re_bptr", re_cbptr.to_int())
                    .put("re_name", re_cbptr->name());
                _g.remove_read(re_cbptr.unconst());
            }
        });
    LOG("Read_Merger", info) << ptree("end");
} // Read_Merger::remove_contained

bool Read_Merger::extend_haploid_layout(Traversal_List& l, Traversal_List::iterator it) const
{
    return (extend_haploid_layout_dir(l, it, false)
            and extend_haploid_layout_dir(l, it, true));
}

bool Read_Merger::extend_haploid_layout_dir(Traversal_List& l, Traversal_List::iterator it, bool e_direction) const
{
    LOG("Read_Merger", debug) << ptree("begin").put("e_direction", e_direction);
    while (true)
    {
        // direction crt_anchor -> next_anchor
        bool rel_c_direction = (it->c_direction + e_direction) % 2;
        LOG("Read_Merger", debug) << ptree("loop.start")
            .put("it", *it)
            .put("rel_c_direction", rel_c_direction);
        // compute next anchor, allele, c_direction, and chunk_support
        Allele_Anchor next_anchor;
        Allele_Specifier next_allele;
        bool next_c_direction;
        bool next_rel_c_direction;
        Anchor_Chunk_Support next_chunk_support;
        RC_DSet next_unmappable_chunks;
        if (it->anchor.is_endpoint() and it->anchor.c_right() != rel_c_direction)
        {
            LOG("Read_Merger", debug) << ptree("loop.contig_endpoint");
            //
            // this is the last anchor in the current contig direction
            // we jump to the associated anchor at the end of the contig edge
            //
            if (not it->allele.ce_next_cbptr())
            {
                // out-degree 0
                ASSERT(it->allele.same_orientation());
                LOG("Read_Merger", debug) << ptree("end.0_degree");
                return true;
            }
            next_rel_c_direction = (rel_c_direction == it->allele.same_orientation());
            next_c_direction = (next_rel_c_direction + e_direction) % 2;
            next_anchor = Allele_Anchor(it->allele.ce_next_cbptr(), next_rel_c_direction);
            next_allele = Allele_Specifier(it->anchor.ce_cbptr(), it->allele.same_orientation());
            auto cks_p = advance_chunks(it->chunk_support.at(it->allele), rel_c_direction);
            ASSERT(cks_p.first.size() == it->chunk_support.at(it->allele).size());
            next_chunk_support.insert(make_pair(next_allele, move(cks_p.first)));
            next_unmappable_chunks = move(cks_p.second);
        }
        else
        {
            LOG("Read_Merger", debug) << ptree("loop.not_contig_endpoint");
            //
            // next anchor is one the same contig
            //
            next_anchor = it->anchor.get_sibling(rel_c_direction);
            next_rel_c_direction = rel_c_direction;
            next_c_direction = it->c_direction;
            LOG("Read_Merger", debug) << ptree("loop").put("next_anchor", next_anchor);
            // compute support at next anchor
            next_chunk_support = next_anchor.chunk_support();
            for (auto& p : next_chunk_support)
            {
                // orient chunks per c_direction
                p.second.reverse(next_c_direction);
                // keep only chunks supporting current allele
                p.second.intersect(it->chunk_support.at(it->allele), false);
            }
            filter_cont(
                next_chunk_support,
                [] (Anchor_Chunk_Support::const_iterator cit) {
                    return not cit->second.empty();
                });
            LOG("Read_Merger", debug) << ptree("loop")
                .put("next_chunk_support", next_chunk_support);
            if (next_chunk_support.empty())
            {
                LOG("Read_Merger", debug) << ptree("end.support_end");
                return true;
            }
            else if (next_chunk_support.size() > 1)
            {
                LOG("Read_Merger", debug) << ptree("loop.multiple_next_allele_support")
                    .put("next_alleles", ba::keys(next_chunk_support));
                // check if there is a unique allele with support larger than the discordant threshold
                set< Allele_Specifier > next_alleles_strong_support;
                for (const auto& p : next_chunk_support)
                {
                    // each previously merged read counts as 5 when computing discordance
                    auto support_accum = [&] (unsigned s, Read_Chunk_CBPtr rc_cbptr) {
                        return s + (rc_cbptr->re_bptr()->name().substr(0, 5) == "merge"? _merged_weight : 1);
                    };
                    unsigned allele_support_size = 0;
                    allele_support_size = accumulate(p.second[0], allele_support_size, support_accum);
                    allele_support_size = accumulate(p.second[1], allele_support_size, support_accum);
                    ASSERT(allele_support_size > 0);
                    if (allele_support_size > _max_discordant_support)
                    {
                        next_alleles_strong_support.insert(p.first);
                    }
                }
                if (next_alleles_strong_support.empty())
                {
                    LOG("Read_Merger", debug) << ptree("end.no_allele_with_strong_support");
                    return false;
                }
                else if (next_alleles_strong_support.size() > 1)
                {
                    LOG("Read_Merger", debug) << ptree("end.multiple_alleles_with_strong_support")
                        .put("next_anchor", next_anchor)
                        .put("next_alleles_strong_support", next_alleles_strong_support);
                    return false;
                }
                else // next_alleles_strong_support.size() == 1
                {
                    next_allele = *next_alleles_strong_support.begin();
                    LOG("Read_Merger", debug) << ptree("loop.single_allele_with_strong_support")
                        .put("next_allele", next_allele);
                }
            }
            else // next_chunk_support.size() == 1
            {
                next_allele = next_chunk_support.begin()->first;
            }
        } // if next anchor on this or next CE

        Traversal_Struct ts{next_anchor, next_allele, next_c_direction,
                move(next_chunk_support), move(next_unmappable_chunks)};

        ASSERT(all_of(
                   ts.chunk_support,
                   [&] (const Anchor_Chunk_Support::value_type& p) {
                       return (all_of(
                                   p.second.at(0),
                                   [&] (Read_Chunk_CBPtr rc_cbptr) {
                                       return rc_cbptr->get_rc() == (0 + ts.c_direction) % 2;
                                   }) and
                               all_of(
                                   p.second.at(1),
                                   [&] (Read_Chunk_CBPtr rc_cbptr) {
                                       return rc_cbptr->get_rc() == (1 + ts.c_direction) % 2;
                                   }));
                   }));

        if (not e_direction)
        {
            ASSERT(it == prev(l.end()));
            l.emplace_back(move(ts));
            it = prev(l.end());
        }
        else
        {
            ASSERT(it == l.begin());
            l.emplace_front(move(ts));
            it = l.begin();
        }
    } // while
} // Read_Merger::extend_haploid_layout_dir

void Read_Merger::split_diverging_reads(Traversal_List& l, Traversal_List::iterator it) const
{
    ASSERT(it != l.end());
    for (int e_direction = 0; e_direction < 2; ++e_direction)
    {
        LOG("Read_Merger", debug) << ptree("begin").put("e_direction", e_direction);
        while (not e_direction? it != prev(l.end()) : it != l.begin())
        {
            auto prev_it = it;
            it = not e_direction? next(it) : prev(it);
            // direction prev_it -> it
            bool rel_c_direction = (it->c_direction + e_direction) % 2;
            ASSERT(not it->chunk_support.empty());
            if (it->chunk_support.size() > 1)
            {
                // layout may only diverge between consecutive anchors in the same contig
                ASSERT(it->anchor.ce_cbptr() == prev_it->anchor.ce_cbptr());
                LOG("Read_Merger", debug) << ptree("loop.divergence")
                    .put("prev_it", *prev_it)
                    .put("it", *it);
                // cut reads which observe an allele other than it->allele
                for (auto& p : it->chunk_support)
                {
                    if (p.first == it->allele) continue;
                    for (int d = 0; d < 2; ++d)
                    {
                        for (auto rc_cbptr : p.second[d])
                        {
                            ASSERT(rc_cbptr->get_rc() == (d + it->c_direction) % 2);
                            bool r_direction = rc_cbptr->get_rc();
                            auto re_cbptr = rc_cbptr->re_bptr();
                            // we make an overlapping cut of re_bptr between prev_anchor and crt_anchor
                            LOG("Read_Merger", debug) << ptree("loop.split_read")
                                .put("re_name", re_cbptr->name())
                                .put("re_bptr", re_cbptr.to_int())
                                .put("rc_bptr", rc_cbptr.to_int())
                                .put("l_anchor", not rel_c_direction? prev_it->anchor : it->anchor)
                                .put("r_anchor", not rel_c_direction? it->anchor : prev_it->anchor);
                            auto rp = _g.split_read(rc_cbptr, not rel_c_direction? prev_it->anchor : it->anchor);
                            auto new_rc_cbptr = (r_direction == rel_c_direction
                                                 ? &*rp.first->chunk_cont().rbegin()
                                                 : &*rp.second->chunk_cont().begin());
                            if (new_rc_cbptr != rc_cbptr)
                            {
                                replace_chunk(l, prev_it, not e_direction, rc_cbptr, new_rc_cbptr, d);
                            }
                        } // for rc_bptr
                    } // for d
                } // for p
                filter_cont(
                    it->chunk_support,
                    [&] (Anchor_Chunk_Support::const_iterator cit) {
                        return cit->first == it->allele;
                    });
            } // if chunk_support.size > 1
        } // while
        LOG("Read_Merger", debug) << ptree("end").put("e_direction", e_direction);
    } // for e_direction
    ASSERT(all_of(
               l,
               [] (const Traversal_Struct& s) {
                   return s.chunk_support.size() == 1 and s.chunk_support.count(s.allele);
               }));
} // Read_Merger::split_diverging_reads

void Read_Merger::merge_reads(Traversal_List& l, Traversal_List::iterator merge_start_it) const
{
    LOG("Read_Merger", debug) << ptree("begin")
        .put("traversal_list", l);

    // create new name, new re
    static unsigned merge_id = 0;
    ostringstream os;
    os << "merge-" << setfill('0') << setw(9) << merge_id++;
    auto m_re_bptr = Read_Entry_Fact::new_elem(string(os.str()), 0);
    deque< Read_Chunk_BPtr > m_chunk_cont;
    Read_Chunk_BPtr m_rc_bptr;
    set< Contig_Entry_BPtr > m_ce_set;
    //
    // create chunk corresponding to merge_start_it
    //
    m_rc_bptr = merge_contig_chunks(merge_start_it->chunk_support.at(merge_start_it->allele),
                                    m_re_bptr, merge_start_it->c_direction);
    m_chunk_cont.push_back(m_rc_bptr);
    //
    // traverse left&right
    //
    for (int e_direction = 0; e_direction < 2; ++e_direction)
    {
        auto m_chunk_inserter = [&] (Read_Chunk_BPtr rc_bptr) {
            not e_direction? m_chunk_cont.push_back(rc_bptr) : m_chunk_cont.push_front(rc_bptr);
        };
        auto it_last = [&] () { return not e_direction? prev(l.end()) : l.begin(); };
        auto it_advance = [&] (Traversal_List::iterator it) { return not e_direction? next(it) : prev(it); };

        auto it = merge_start_it;
        while (true)
        {
            // advance it to end of current contig
            bool rel_c_direction = (e_direction + it->c_direction) % 2;
            auto prev_it = it;
            while (it != it_last()
                   and (not it->anchor.is_endpoint() or it->anchor.c_right() == rel_c_direction))
            {
                ASSERT(it->anchor.ce_cbptr() == prev_it->anchor.ce_cbptr());
                it = it_advance(it);
            }
            if (it == it_last())
            {
                break;
            }
            ASSERT(it->anchor.ce_cbptr() == prev_it->anchor.ce_cbptr());
            ASSERT(prev_it == merge_start_it
                   or (prev_it->anchor.is_endpoint() and prev_it->anchor.c_right() == rel_c_direction));
            ASSERT(it->anchor.is_endpoint() and not it->anchor.c_right() == rel_c_direction);
            ASSERT(all_of(not e_direction? prev_it : it,
                          not e_direction? next(it) : next(prev_it),
                          [&] (const Traversal_Struct& ts) {
                              return (ts.anchor.ce_cbptr() == prev_it->anchor.ce_cbptr()
                                      and ts.c_direction == prev_it->c_direction
                                      and ts.chunk_support.size() == 1
                                      and ts.chunk_support.count(ts.allele)
                                      and (RC_DSet::intersect(
                                               ts.chunk_support.at(ts.allele),
                                               prev_it->chunk_support.at(prev_it->allele), false).size()
                                           == ts.chunk_support.at(ts.allele).size()));
                          }));
            prev_it = it;
            it = it_advance(it);
            rel_c_direction = (e_direction + it->c_direction) % 2;
            ASSERT(it->anchor.is_endpoint() and it->anchor.c_right() == rel_c_direction);
            ASSERT(it->chunk_support.at(it->allele).size() == prev_it->chunk_support.at(prev_it->allele).size());
            LOG("Read_Merger", debug) << ptree("loop.begin")
                .put("e_direction", e_direction)
                .put("it", *it);
            if (not it->unmappable_chunks.empty())
            {
                if (it->unmappable_chunks.size() < it->chunk_support.at(it->allele).size())
                {
                    // inconsistency:
                    // from the set of reads being merged,
                    // some contain an unmappable region but others do not;
                    // print an informative message
                    auto prev_ce_cbptr = prev_it->anchor.ce_cbptr();
                    // compute next ce_bptr
                    auto next_ce_cbptr = it->anchor.ce_cbptr();
                    // compute reads which have unmappable region
                    RE_DSet unmappable_re_dset;
                    for (int d = 0; d < 2; ++d)
                    {
                        for (auto rc_bptr : it->unmappable_chunks[d])
                        {
                            unmappable_re_dset[d].insert(rc_bptr->re_bptr());
                        }
                    }
                    // compute reads without unmappable region
                    RE_DSet no_unmappable_re_dset;
                    for (int d = 0; d < 2; ++d)
                    {
                        for (auto rc_bptr : it->chunk_support.at(it->allele)[d])
                        {
                            no_unmappable_re_dset[d].insert(rc_bptr->re_bptr());
                        }
                    }
                    no_unmappable_re_dset.subtract(unmappable_re_dset, false);
                    LOG("Read_Merger", warning) << ptree("unmappable_inconsistency")
                        .put("prev_ce_cbptr", prev_ce_cbptr.to_int())
                        .put("next_ce_cbptr", next_ce_cbptr.to_int())
                        .put("unmappable_re_dset", unmappable_re_dset)
                        .put("no_unmappable_re_dset", no_unmappable_re_dset);
                } // if unmappable_chunks.size < no_unmappable_chunks.size
                //
                // merge unmappable chunks between contigs
                //
                m_rc_bptr = merge_unmappable_chunks(it->unmappable_chunks, m_re_bptr);
                m_chunk_inserter(m_rc_bptr);
            }
            //
            // merge chunks from next contig
            //
            m_rc_bptr = merge_contig_chunks(it->chunk_support.at(it->allele), m_re_bptr, it->c_direction);
            m_chunk_inserter(m_rc_bptr);
        } // while true
    } // for e_direction

    // now all the new chunks are in m_chunk_cont;
    // compute the m_read offsets and add the chunks to m_re chunk_cont
    Size_Type pos = 0;
    for (auto rc_bptr : m_chunk_cont)
    {
        rc_bptr->r_start() = pos;
        pos += rc_bptr->r_len();
        m_re_bptr->chunk_cont().insert(rc_bptr);
    }
    m_re_bptr->len() = pos;
    _g.re_cont().insert(m_re_bptr);
    _g.check(set< Read_Entry_CBPtr >{ m_re_bptr });

    // finally, destroy the old reads
    for (int d = 0; d < 2; ++d)
    {
        for (auto rc_cbptr : merge_start_it->chunk_support.at(merge_start_it->allele)[d])
        {
            _g.remove_read(rc_cbptr->re_bptr().unconst());
        }
    }
    _g.check(set< Read_Entry_CBPtr >{ m_re_bptr });

    LOG("Read_Merger", info) << ptree("end");
} // Read_Merger::merge_reads

Read_Chunk_BPtr
Read_Merger::merge_contig_chunks(const RC_DSet& rc_dset, Read_Entry_BPtr m_re_bptr, bool c_direction) const
{
    LOG("Read_Merger", debug) << ptree("begin")
        .put("rc_dset", rc_dset);
    ASSERT(not (rc_dset[0].empty() and rc_dset[1].empty()));
    auto ce_bptr = not rc_dset[0].empty()
        ? (*rc_dset[0].begin())->ce_bptr()
        : (*rc_dset[1].begin())->ce_bptr();
    ASSERT(all_of(
               rc_dset[0],
               [&] (Read_Chunk_CBPtr rc_cbptr) {
                   return rc_cbptr->ce_bptr() == ce_bptr and rc_cbptr->get_rc() == (0 + c_direction) % 2;
               })
           and all_of(
               rc_dset[1],
               [&] (Read_Chunk_CBPtr rc_cbptr) {
                   return rc_cbptr->ce_bptr() == ce_bptr and rc_cbptr->get_rc() == (1 + c_direction) % 2;
               }));
    // initialize position map
    map< Read_Chunk_CBPtr, Read_Chunk_Pos > pos_map;
    for (int d = 0; d < 2; ++d)
    {
        for (auto rc_cbptr : rc_dset[d])
        {
            pos_map.insert(make_pair(rc_cbptr, rc_cbptr->get_start_pos()));
        }
    }
    // set read start to be 0; will be fixed once we have all chunks
    Size_Type r_pos = 0;
    // set contig start as minimum contig position
    Size_Type c_pos = min_value_of(
        pos_map,
        [] (const decltype(pos_map)::value_type& p) { return p.second.c_pos; });
    // initialize merged chunk
    auto m_rc_bptr = Read_Chunk_Fact::new_elem(r_pos, 0, c_pos, 0, c_direction);
    m_rc_bptr->re_bptr() = m_re_bptr;
    m_rc_bptr->ce_bptr() = ce_bptr;
    while (not pos_map.empty())
    {
        // all positions are before or at c_pos
        ASSERT(c_pos <= min_value_of(
                   pos_map,
                   [] (const decltype(pos_map)::value_type& p) { return p.second.c_pos; }));
        // check if there is any chunk with position at c_pos, that has a mutation next
        bool is_ins_next = any_of(
            ba::values(pos_map),
            [&] (const Read_Chunk_Pos& pos) {
                return (pos.c_pos == c_pos and pos.get_match_len() == 0 and not pos.past_last_mut() and pos.mut().is_ins());
            });
        // find active chunks: the ones whose position is at c_pos
        set< Read_Chunk_CBPtr > active_rc_set;
        for (const auto& p : pos_map)
        {
            if (p.second.c_pos != c_pos) continue;
            // exception: allow terminal chunks at their start position to not have an insertion
            if (is_ins_next
                and p.second == p.first->get_start_pos()
                and p.second.get_match_len() > 0
                and not p.first->re_bptr()->chunk_cont().get_sibling(p.first, false, false))
            {
                continue;
            }
            active_rc_set.insert(p.first);
        }
        ASSERT(not active_rc_set.empty());
        // compute next breakpoint
        auto match_len = min_value_of(
            ba::values(pos_map),
            [&] (const Read_Chunk_Pos& pos) {
                return (pos.c_pos == c_pos
                        ? pos.get_match_len()
                        : pos.c_pos - c_pos);
            });
        LOG("Read_Merger", debug) << ptree("loop.before_increment")
            .put("c_pos", c_pos)
            .put("is_ins_next", is_ins_next)
            .put("pos_map", pos_map)
            .put("active_rc_set", active_rc_set)
            .put("match_len", match_len);
        // advance active set
        if (match_len > 0)
        {
            // match stretch follows
            br::for_each(
                active_rc_set,
                [&] (Read_Chunk_CBPtr rc_cbptr) { pos_map.at(rc_cbptr).increment(c_pos + match_len, true); });
            Size_Type next_c_pos = min_value_of(
                ba::values(pos_map),
                [] (const Read_Chunk_Pos& pos) { return pos.c_pos; });
            ASSERT(c_pos < next_c_pos);
            r_pos += next_c_pos - c_pos;
            c_pos = next_c_pos;
        }
        else
        {
            // mutation follows
            auto mut_bptr = pos_map.at(*active_rc_set.begin()).mca_cit->mut_cbptr().unconst();
            // it must be the same in all active chunks
            ASSERT(all_of(
                       active_rc_set,
                       [&] (Read_Chunk_CBPtr rc_cbptr) {
                           return pos_map.at(rc_cbptr).mca_cit->mut_cbptr() == mut_bptr;
                       }));
            // none of the inactive chunks may be at a position inside the mutation
            ASSERT(all_of(
                       pos_map,
                       [&] (const decltype(pos_map)::value_type& p) {
                           return active_rc_set.count(p.first) or p.second.c_pos >= mut_bptr->rf_end();
                       }));
            // increment position for active chunks
            br::for_each(
                active_rc_set,
                [&] (Read_Chunk_CBPtr rc_cbptr) { pos_map.at(rc_cbptr).increment(); });
            // they now must all be at position right after mutation
            ASSERT(all_of(
                       active_rc_set,
                       [&] (Read_Chunk_CBPtr rc_cbptr) {
                           return pos_map.at(rc_cbptr).c_pos == mut_bptr->rf_end();
                       }));
            // add mutation to merged chunk
            auto mca_bptr = Mutation_Chunk_Adapter_Fact::new_elem(mut_bptr, m_rc_bptr);
            mut_bptr->chunk_ptr_cont().insert(mca_bptr);
            m_rc_bptr->mut_ptr_cont().push_back(mca_bptr);
            r_pos += mut_bptr->seq_len();
            c_pos += mut_bptr->rf_len();
        }
        LOG("Read_Merger", debug) << ptree("loop.after_increment")
            .put("c_pos", c_pos)
            .put("pos_map", range_to_ptree(pos_map, [] (const decltype(pos_map)::value_type& p) {
                        return ptree().put("rc_cbptr", p.first).put("pos", p.second);
                    }));
        // remove chunks that reached their end position
        filter_cont(
            pos_map,
            [&] (const decltype(pos_map)::value_type& p) { return p.second != p.first->get_end_pos(); });
        LOG("Read_Merger", debug) << ptree("loop.after_filter")
            .put("pos_map", range_to_ptree(pos_map, [] (const decltype(pos_map)::value_type& p) {
                        return ptree().put("rc_cbptr", p.first).put("pos", p.second);
                    }));
    } // while pos_map non-empty
    // fix r_len and c_len
    m_rc_bptr->r_len() = r_pos - m_rc_bptr->r_start();
    m_rc_bptr->c_len() = c_pos - m_rc_bptr->c_start();
    // add new chunk to the contig chunk cont
    ce_bptr->chunk_cont().insert(m_rc_bptr);
    LOG("Read_Merger", debug) << ptree("end");
    return m_rc_bptr;
} // Read_Merger::merge_contig_chunks

Read_Chunk_BPtr Read_Merger::merge_unmappable_chunks(const RC_DSet& unmappable_chunks, Read_Entry_BPtr m_re_bptr) const
{
    seqan::FragmentStore<> store;
    seqan::ConsensusAlignmentOptions options;
    options.useContigID = true;
    vector< Seq_Type > seq_v;
    for (int dir = 0; dir < 2; ++dir)
    {
        for (auto rc_bptr : unmappable_chunks[dir])
        {
            ASSERT(rc_bptr->ce_bptr()->is_unmappable());
            Seq_Type seq = rc_bptr->ce_bptr()->seq().revcomp(dir);
            auto id = seqan::appendRead(store, string("NNNNNNNNNNNNNNNNNNNN") + seq);
            seqan::appendAlignedRead(store, id, 0, 0, int(seq.size()));
            seq_v.emplace_back(move(seq));
        }
    }
    seqan::consensusAlignment(store, options);
    seqan::reAlignment(store, 0, 1, 10, false);
    ostringstream os;
    os << store.contigStore[0].seq;
    string res = os.str();
    auto prefix_end = res.find_first_not_of('N');
    ASSERT(prefix_end != string::npos);
    auto suffix_start = res.find_first_of('N', prefix_end);
    if (suffix_start != string::npos)
    {
        LOG("Read_Merger", warning) << ptree("consensus_alignment_failed")
            .put("unmappable_chunks", seq_v)
            .put("raw_result", vector< string >{res})
            .put("result", vector< string >{res.substr(prefix_end, suffix_start - prefix_end)});
    }
    auto m_ce_bptr = Contig_Entry_Fact::new_elem(Seq_Type(res.substr(prefix_end, suffix_start - prefix_end)));
    m_ce_bptr->set_unmappable();
    _g.ce_cont().insert(m_ce_bptr);
    auto m_rc_bptr = Read_Chunk_Fact::new_elem(0, m_ce_bptr->len(), 0, m_ce_bptr->len(), false);
    m_rc_bptr->re_bptr() = m_re_bptr;
    m_rc_bptr->ce_bptr() = m_ce_bptr;
    m_rc_bptr->set_unbreakable(true);
    m_ce_bptr->chunk_cont().insert(m_rc_bptr);
    LOG("Read_Merger", debug) << ptree("end")
        .put("unmappable_chunks", seq_v)
        .put("result", vector< string >{m_ce_bptr->seq()});
    return m_rc_bptr;
}

ptree Read_Merger::Traversal_Struct::to_ptree() const
{
    return ptree()
        .put("anchor", anchor)
        .put("allele", allele)
        .put("c_direction", c_direction)
        .put("chunk_support", chunk_support)
        .put("unmappable_chunks", unmappable_chunks);
} // Read_Merger::Traversal_Struct::to_ptree

} // namespace MAC
