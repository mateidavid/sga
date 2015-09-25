#include "Read_Merger.hpp"

#include <seqan/store.h>
#include <seqan/consensus.h>

#include "filter_cont.hpp"

namespace MAC
{

void Read_Merger::operator () ()
{
    LOG("Read_Merger", info) << ptree("begin");
    for (auto ce_bptr : _g.ce_cont() | referenced)
    {
        if (not ce_bptr->is_normal()) continue;
        for (auto mut_bptr : ce_bptr->mut_cont() | referenced)
        {
            Allele_Anchor anchor(mut_bptr);
            auto anchor_support = get_anchor_read_support(anchor, false);
            for (int al = 0; al < 2; ++al)
            {
                if (mut_bptr->copy_num(al) != 1) continue;
                // found a mutation allele with copy number 1
                // greb all reads supporting it
                Allele_Specifier allele(al);
                LOG("Read_Merger", debug) << ptree("haploid_mutation")
                    .put("anchor", anchor)
                    .put("allele", allele);
                auto allele_support = get_allele_read_support(move(anchor_support), allele);
                extend_haploid_support(anchor, allele, false, allele_support);
                extend_haploid_support(anchor, allele, true, allele_support);
                merge_reads(ce_bptr, allele_support);
            }
        }
    }
    LOG("Read_Merger", info) << ptree("end");
}

void Read_Merger::extend_haploid_support(
    const Allele_Anchor& init_anchor, const Allele_Specifier& init_allele, bool init_c_direction,
    RE_OSet& allele_support)
{
    LOG("Read_Merger", debug) << ptree("begin")
        .put("anchor", init_anchor)
        .put("allele", init_allele)
        .put("c_direction", init_c_direction)
        .put("allele_support", allele_support);
    auto crt_anchor = init_anchor;
    auto crt_allele = init_allele;
    auto crt_c_direction = init_c_direction;
    while (true)
    {
        LOG("Read_Merger", debug) << ptree("loop_start")
            .put("crt_anchor", crt_anchor)
            .put("crt_allele", crt_allele)
            .put("crt_c_direction", crt_c_direction);
        // compute next anchor
        Allele_Anchor next_anchor;
        Allele_Specifier next_allele;
        if (crt_anchor.is_endpoint() and crt_anchor.c_right() != crt_c_direction)
        {
            // this is the last anchor in the current contig direction
            // we jump to the associated anchor at the end of the contig edge
            if (not crt_allele.ce_next_cbptr())
            {
                // out-degree 0
                ASSERT(crt_allele.same_orientation());
                LOG("Read_Merger", debug) << ptree("loop_end__0_degree");
                return;
            }
            bool next_c_direction = crt_c_direction == crt_allele.same_orientation();
            next_anchor = Allele_Anchor(crt_allele.ce_next_cbptr(), next_c_direction);
            next_allele = Allele_Specifier(crt_anchor.ce_cbptr(), crt_allele.same_orientation());
            swap(crt_c_direction, next_c_direction);
            swap(crt_anchor, next_anchor);
            swap(crt_allele, next_allele);
            LOG("Read_Merger", debug) << ptree("loop_next__contig_endpoint");
            continue;
        }
        next_anchor = crt_anchor.get_sibling(crt_c_direction);
        LOG("Read_Merger", debug) << ptree("loop").put("next_anchor", next_anchor);
        // compute support at next anchor
        auto next_anchor_support = get_anchor_read_support(next_anchor, crt_c_direction);
        // remove all reads not currently being tracked
        for (auto& p : next_anchor_support)
        {
            for (int dir = 0; dir < 2; ++dir)
            {
                filter_cont(
                    p.second[dir],
                    [&] (Read_Entry_CBPtr re_cbptr) {
                        return allele_support[(dir + init_c_direction) % 2].count(re_cbptr) > 0;
                    });
            }
        }
        filter_cont(
            next_anchor_support,
            [] (Anchor_Read_Support::const_iterator cit) {
                return not (cit->second[0].empty() and cit->second[1].empty());
            });
        LOG("Read_Merger", debug) << ptree("loop")
            .put("next_anchor_support", next_anchor_support);
        if (next_anchor_support.empty())
        {
            LOG("Read_Merger", debug) << ptree("loop_end__support_end");
            break;
        }
        else if (next_anchor_support.size() > 1)
        {
            LOG("Read_Merger", debug) << ptree("loop__multiple_next_allele_support")
                .put("next_alleles", ba::keys(next_anchor_support));
            // check if there is a unique allele with support larger than the discordant threshold
            bool found_next_allele = false;
            for (const auto& p : next_anchor_support)
            {
                unsigned allele_support_size = p.second[0].size() + p.second[1].size();
                ASSERT(allele_support_size > 0);
                if (allele_support_size > _max_discordant_support)
                {
                    if (found_next_allele)
                    {
                        // this is a second allele with large support; we give up
                        LOG("Read_Merger", info) << ptree("loop_end__multiple_alleles_with_strong_support")
                            .put("next_allele_1", next_allele)
                            .put("next_allele_2", p.first);
                        allele_support[0].clear();
                        allele_support[1].clear();
                        return;
                    }
                    else
                    {
                        next_allele = p.first;
                        found_next_allele = true;
                    }
                }
            }
            if (not found_next_allele)
            {
                LOG("Read_Merger", info) << ptree("loop_end__no_allele_with_strong_support");
                allele_support[0].clear();
                allele_support[1].clear();
                return;
            }
            // support at next_anchor is ambiguous,
            // but a single allele has support greater than the discordant threshold
            LOG("Read_Merger", info) << ptree("loop")
                .put("next_allele", next_allele);
            for (auto& p : next_anchor_support)
            {
                unsigned allele_support_size = p.second[0].size() + p.second[1].size();
                ASSERT(allele_support_size > 0);
                if (allele_support_size > _max_discordant_support) continue;
                // we split the remaining reads in between anchor and next_anchor
                for (int dir = 0; dir < 2; ++dir)
                {
                    for (auto re_cbptr : p.second[dir])
                    {
                        ASSERT(allele_support[(dir + init_c_direction) % 2].count(re_cbptr));
                        allele_support[(dir + init_c_direction) % 2].erase(re_cbptr);
                        auto rp = split_read(
                            re_cbptr,
                            not crt_c_direction? crt_anchor : next_anchor,
                            not crt_c_direction? next_anchor : crt_anchor);
                        // we keep the side of the read that spans crt_anchor
                        auto tmp = get_allele_read_support(crt_anchor, crt_allele, crt_c_direction);
                        (void)tmp;
                        //TODO: remove
                        ASSERT(tmp[dir].count(rp.first) > 0);
                        allele_support[(dir + init_c_direction) % 2].insert(rp.first);
                    }
                    p.second[dir].clear();
                }
            }
            filter_cont(
                next_anchor_support,
                [] (Anchor_Read_Support::const_iterator cit) {
                    return not (cit->second[0].empty() and cit->second[1].empty());
                });
        } // if (next_anchor_support.size() > 1)
        ASSERT(next_anchor_support.size() == 1);
        crt_anchor = next_anchor;
        crt_allele = next_anchor_support.begin()->first;
    }
} // Read_Merger::extend_haploid_support

void Read_Merger::merge_reads(Contig_Entry_BPtr ce_bptr, const RE_OSet& re_oset)
{
    if (re_oset[0].size() + re_oset[1].size() < 2) return;
    LOG("Read_Merger", info) << ptree("begin")
        .put("ce_bptr", ce_bptr.to_int())
        .put("re_oset[0]", re_oset[0])
        .put("re_oset[1]", re_oset[1]);

    // create new name, new re
    static unsigned merge_id = 0;
    ostringstream os;
    os << "merge-" << setfill('0') << setw(9) << merge_id++;
    auto m_re_bptr = Read_Entry_Fact::new_elem(string(os.str()), 0);
    deque< Read_Chunk_BPtr > m_chunk_cont;
    Read_Chunk_BPtr m_rc_bptr;
    set< Contig_Entry_BPtr > m_ce_set;
    for (int dir = 0; dir < 2; ++dir)
    {
        auto m_chunk_inserter = [&] (Read_Chunk_BPtr rc_bptr) {
            not dir? m_chunk_cont.push_back(rc_bptr) : m_chunk_cont.push_front(rc_bptr);
        };
        RC_OSet crt_chunks = get_oriented_chunks(ce_bptr, re_oset);
        ASSERT(crt_chunks[0].size() == re_oset[0].size()
               and crt_chunks[1].size() == re_oset[1].size());
        if (dir == 0)
        {
            // merge chunks in start contig
            m_rc_bptr = merge_contig_chunks(crt_chunks, m_re_bptr);
            m_chunk_inserter(m_rc_bptr);
            m_ce_set.insert(m_rc_bptr->ce_bptr());
        }
        bool crt_dir = dir;
        while (true)
        {
            // advance chunks but save unmappable chunks
            RC_OSet unmappable_chunks;
            tie(crt_chunks, unmappable_chunks) = advance_chunks(crt_chunks, crt_dir);
            if (not (unmappable_chunks[0].empty() and unmappable_chunks[1].empty()))
            {
                // reads may not end with unmappable regions
                if (unmappable_chunks[0].size() + unmappable_chunks[1].size()
                    < crt_chunks[0].size() + crt_chunks[1].size())
                {
                    // inconsistency:
                    // from the set of reads being merged,
                    // some contain an unmappable region but others do not
                    // compute previous ce_bptr
                    // will abort, but first print an informative message
                    int prev_rc_dir = unmappable_chunks[0].empty();
                    auto prev_rc_cbptr = *unmappable_chunks[prev_rc_dir].begin();
                    prev_rc_cbptr = prev_rc_cbptr->re_bptr()->chunk_cont().get_sibling(
                        prev_rc_cbptr, true, not ( (dir + prev_rc_dir + 1) % 2 ));
                    auto prev_ce_cbptr = prev_rc_cbptr->ce_bptr();
                    // compute next ce_bptr
                    auto next_ce_cbptr = not crt_chunks[0].empty()
                        ? (*crt_chunks[0].begin())->ce_bptr()
                        : (*crt_chunks[1].begin())->ce_bptr();
                    // compute reads which have unmappable region
                    RE_OSet unmappable_oset;
                    for (int d = 0; d < 2; ++d)
                    {
                        for (auto rc_bptr : unmappable_chunks[d])
                        {
                            unmappable_oset[d].insert(rc_bptr->re_bptr());
                        }
                    }
                    // compute reads without unmappable region
                    RE_OSet non_unmappable_oset;
                    for (int d = 0; d < 2; ++d)
                    {
                        for (auto rc_bptr : crt_chunks[d])
                        {
                            if (unmappable_oset[d].count(rc_bptr->re_bptr()) == 0)
                            {
                                non_unmappable_oset[d].insert(rc_bptr->re_bptr());
                            }
                        }
                    }
                    LOG("Read_Merger", warning) << ptree("unmappable_inconsistency")
                        .put("prev_ce_cbptr", prev_ce_cbptr.to_int())
                        .put("next_ce_cbptr", next_ce_cbptr.to_int())
                        .put("unmappable_oset[0]", unmappable_oset[0])
                        .put("unmappable_oset[1]", unmappable_oset[1])
                        .put("non_unmappable_oset[0]", non_unmappable_oset[0])
                        .put("non_unmappable_oset[1]", non_unmappable_oset[1]);
                }
                m_rc_bptr = merge_unmappable_chunks(unmappable_chunks, m_re_bptr);
                m_chunk_inserter(m_rc_bptr);
            } // if not unmappable_chunks.empty
            if (crt_chunks[0].empty() and crt_chunks[1].empty())
            {
                break;
            }
            m_rc_bptr = merge_contig_chunks(crt_chunks, m_re_bptr);
            m_chunk_inserter(m_rc_bptr);
            auto p = m_ce_set.insert(m_rc_bptr->ce_bptr());
            if (not p.second)
            {
                LOG("Read_Merger", warning) << ptree("contig_entry_duplication")
                    .put("ce_bptr", m_rc_bptr->ce_bptr().to_int());
            }
        } // while true
    } // for dir

    // now all the new chunks are in m_chunk_cont;
    // compute the m_read offsets and add the chunks to m_re chunk_cont
    Size_Type pos = 0;
    for (auto rc_bptr : m_chunk_cont)
    {
        rc_bptr->r_start() = pos;
        pos += rc_bptr->r_len();
        m_re_bptr->chunk_cont().insert(rc_bptr);
    }
    m_re_bptr->check();

    // finally, destroy the old reads
    for (int d = 0; d < 2; ++d)
    {
        for (auto re_cbptr : re_oset[d])
        {
            _g.remove_read(re_cbptr.unconst());
        }
    }

    LOG("Read_Merger", info) << ptree("end");
}

pair< Read_Entry_CBPtr, Read_Entry_CBPtr >
Read_Merger::split_read(Read_Entry_CBPtr re_cbptr, const Allele_Anchor& l_anchor, const Allele_Anchor& r_anchor)
{
    LOG("Read_Merger", debug) << ptree("begin")
        .put("re_bptr", re_cbptr.to_int())
        .put("re_name", re_cbptr->name())
        .put("l_anchor", l_anchor)
        .put("r_anchor", r_anchor);
    ASSERT(re_cbptr);
    auto re_bptr = re_cbptr.unconst();
    ASSERT(l_anchor.ce_cbptr() == r_anchor.ce_cbptr());
    ASSERT(l_anchor.get_sibling(false) == r_anchor);
    auto ce_bptr = l_anchor.ce_cbptr().unconst();
    auto l_out_pos = l_anchor.is_endpoint()? 0 : l_anchor.mut_cbptr()->rf_start();
    (void)l_out_pos;
    auto l_in_pos = l_anchor.is_endpoint()? 0 : l_anchor.mut_cbptr()->rf_end();
    auto r_in_pos = r_anchor.is_endpoint()? ce_bptr->len() : r_anchor.mut_cbptr()->rf_start();
    auto r_out_pos = r_anchor.is_endpoint()? ce_bptr->len() : r_anchor.mut_cbptr()->rf_end();
    (void)r_out_pos;
    ASSERT(l_in_pos < r_in_pos);
    //auto match_len = r_in_pos - l_in_pos;
    auto rc_bptr = ce_bptr->chunk_cont().search_read(re_bptr).unconst();
    auto dir = rc_bptr->get_rc();
    ASSERT(rc_bptr);
    ASSERT(rc_bptr->get_c_start() <= l_out_pos);
    ASSERT(r_out_pos <= rc_bptr->get_c_end());
    // create new Read_Entry object that will hold the tail
    _g.re_cont().erase(re_bptr);
    auto re_p = Read_Entry::split(rc_bptr, l_in_pos, r_in_pos);
    _g.re_cont().insert(re_p.first);
    _g.re_cont().insert(re_p.second);
    _g.check(set< Read_Entry_CBPtr >{ re_p.first, re_p.second });
    LOG("Read_Merger", debug) << ptree("end")
        .put("head_re", re_p.first.to_int())
        .put("tail_re", re_p.second.to_int());
    return (not dir ? re_p : make_pair(re_p.second, re_p.first));
} // Read_Merger::split_read

Read_Chunk_BPtr Read_Merger::merge_contig_chunks(const RC_OSet& rc_oset, Read_Entry_BPtr m_re_bptr)
{
    LOG("Read_Merger", debug) << ptree("begin")
        .put("rc_oset", rc_oset);
    ASSERT(not (rc_oset[0].empty() and rc_oset[1].empty()));
    auto ce_bptr = not rc_oset[0].empty()
        ? (*rc_oset[0].begin())->ce_bptr()
        : (*rc_oset[0].begin())->ce_bptr();
    auto c_direction = not rc_oset[0].empty()
        ? (*rc_oset[0].begin())->get_rc()
        : (*rc_oset[1].begin())->get_rc();
    ASSERT(all_of(
               rc_oset[0],
               [&] (Read_Chunk_CBPtr rc_cbptr) {
                   return rc_cbptr->ce_bptr() == ce_bptr and rc_cbptr->get_rc() == (0 + c_direction) % 2;
               })
           and all_of(
               rc_oset[1],
               [&] (Read_Chunk_CBPtr rc_cbptr) {
                   return rc_cbptr->ce_bptr() == ce_bptr and rc_cbptr->get_rc() == (1 + c_direction) % 2;
               }));
    // initialize position map
    map< Read_Chunk_CBPtr, Read_Chunk_Pos > pos_map;
    for (int d = 0; d < 2; ++d)
    {
        for (auto rc_cbptr : rc_oset[d])
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
        // find active chunks: the ones whose position is at c_pos
        RC_Set active_rc_set;
        for (const auto& p : pos_map)
        {
            if (p.second.c_pos == c_pos)
            {
                active_rc_set.insert(p.first);
            }
        }
        ASSERT(not active_rc_set.empty());
        // check that either all or none of the reads in the active set have a mutation next
        auto match_len = min_value_of(
            active_rc_set,
            [&] (Read_Chunk_CBPtr rc_cbptr) { return pos_map.at(rc_cbptr).get_match_len(); });
        LOG("Read_Merger", debug) << ptree("before_increment")

            .put("c_pos", c_pos)
            .put("pos_map", pos_map)
            .put("active_rc_set", active_rc_set)
            .put("match_len", match_len);
        // advance active set
        if (match_len > 0)
        {
            // match stretch follows
            br::for_each(
                active_rc_set,
                [&] (Read_Chunk_CBPtr rc_cbptr) { pos_map.at(rc_cbptr).increment(); });
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
        LOG("Read_Merger", debug) << ptree("after_increment")
            .put("c_pos", c_pos)
            .put("pos_map", range_to_ptree(pos_map, [] (const decltype(pos_map)::value_type& p) {
                        return ptree().put("rc_cbptr", p.first).put("pos", p.second);
                    }));
        // remove chunks that reached their end position
        filter_cont(
            pos_map,
            [&] (const decltype(pos_map)::value_type& p) { return p.second != p.first->get_end_pos(); });
        LOG("Read_Merger", debug) << ptree("after_filter")
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
}

Read_Chunk_BPtr Read_Merger::merge_unmappable_chunks(const RC_OSet& unmappable_chunks, Read_Entry_BPtr m_re_bptr)
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
            auto id = seqan::appendRead(store, string(seq));
            seqan::appendAlignedRead(store, id, 0, 0, int(seq.size()));
            seq_v.emplace_back(move(seq));
        }
    }
    seqan::consensusAlignment(store, options);
    ostringstream os;
    os << store.contigStore[0].seq;
    auto m_ce_bptr = Contig_Entry_Fact::new_elem(Seq_Type(os.str()));
    m_ce_bptr->set_unmappable();
    _g.ce_cont().insert(m_ce_bptr);
    auto m_rc_bptr = Read_Chunk_Fact::new_elem(0, m_ce_bptr->len(), 0, m_ce_bptr->len(), false);
    m_rc_bptr->re_bptr() = m_re_bptr;
    m_rc_bptr->ce_bptr() = m_ce_bptr;
    m_ce_bptr->chunk_cont().insert(m_rc_bptr);
    LOG("Read_Merger", debug) << ptree("end")
        .put("unmappable_chunks", seq_v)
        .put("result", m_ce_bptr->seq());
    return m_rc_bptr;
}

pair< RC_OSet, RC_OSet >
Read_Merger::advance_chunks(const RC_OSet& crt_rc_oset, bool direction)
{
    RC_OSet next_rc_oset;
    RC_OSet unmappable_rc_oset;
    for (int d = 0; d < 2; ++d)
    {
        for (auto rc_cbptr : crt_rc_oset[d])
        {
            ASSERT(rc_cbptr);
            rc_cbptr = rc_cbptr->re_bptr()->chunk_cont().get_sibling(
                rc_cbptr, true, not ( (d + direction) % 2 ));
            if (not rc_cbptr) continue;
            if (rc_cbptr->ce_bptr()->is_unmappable())
            {
                // skip one unmappable chunk if necessary
                unmappable_rc_oset[d].insert(rc_cbptr);
                rc_cbptr = rc_cbptr->re_bptr()->chunk_cont().get_sibling(
                    rc_cbptr, true, not ( (d + direction) % 2 ));
                ASSERT(rc_cbptr);
            }
            next_rc_oset[d].insert(rc_cbptr);
        }
    }
    return make_pair(next_rc_oset, unmappable_rc_oset);
}

} // namespace MAC
