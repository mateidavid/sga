#include "Read_Merger.hpp"

#include <seqan/store.h>
#include <seqan/consensus.h>

#include "filter_cont.hpp"

namespace MAC
{

void Read_Merger::operator () () const
{
    LOG("Read_Merger", info) << ptree("begin")
        .put("max_discordant_support", _max_discordant_support)
        .put("merged_weight", _merged_weight);
    for (bool done = false; not done; )
    {
        done = true;
        for (auto ce_bptr : _g.ce_cont() | referenced)
        {
            if (not ce_bptr->is_normal()) continue;
            for (auto mut_bptr : ce_bptr->mut_cont() | referenced)
            {
                Allele_Anchor anchor(mut_bptr);
                auto anchor_support = anchor.read_support();
                for (int al = 0; al < 2; ++al)
                {
                    if (mut_bptr->copy_num(al) != 1) continue;
                    // found a mutation allele with copy number 1
                    // greb all reads supporting it
                    Allele_Specifier allele(al);
                    auto allele_support = move(anchor_support.at(allele));
                    if (allele_support.size() <= 1) continue;
                    LOG("Read_Merger", debug) << ptree("haploid_mutation")
                        .put("anchor", anchor)
                        .put("allele", allele)
                        .put("allele_support", allele_support);
                    extend_haploid_support(anchor, allele, false, allele_support);
                    extend_haploid_support(anchor, allele, true, allele_support);
                    LOG("Read_Merger", debug) << ptree("after_haploid_extension")
                        .put("allele_support", allele_support);
                    if (allele_support.size() <= 1) continue;
                    merge_reads(ce_bptr, allele_support);
                    done = false;
                }
            }
        }
    }
    LOG("Read_Merger", info) << ptree("end");
}

void Read_Merger::extend_haploid_support(
    const Allele_Anchor& init_anchor, const Allele_Specifier& init_allele, bool init_c_direction,
    RE_DSet& allele_support) const
{
    LOG("Read_Merger", debug) << ptree("begin")
        .put("init_anchor", init_anchor)
        .put("init_allele", init_allele)
        .put("init_c_direction", init_c_direction)
        .put("allele_support", allele_support);
    auto crt_anchor = init_anchor;
    auto crt_allele = init_allele;
    auto crt_c_direction = init_c_direction;
    RE_DSet inactive_allele_support;
    while (true)
    {
        auto active_allele_support = allele_support.subtract(allele_support, inactive_allele_support, false);
        LOG("Read_Merger", debug) << ptree("loop_start")
            .put("crt_anchor", crt_anchor)
            .put("crt_allele", crt_allele)
            .put("crt_c_direction", crt_c_direction)
            .put("active_allele_support", active_allele_support);
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
        auto next_anchor_support = next_anchor.read_support(); //get_anchor_read_support(next_anchor, 1, crt_c_direction);
        // remove all reads not currently being tracked
        for (auto& p : next_anchor_support)
        {
            p.second.reverse(crt_c_direction);
            for (int dir = 0; dir < 2; ++dir)
            {
                filter_cont(
                    p.second[dir],
                    [&] (Read_Entry_CBPtr re_cbptr) {
                        return active_allele_support[(dir + init_c_direction) % 2].count(re_cbptr) > 0;
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
            set< Allele_Specifier > next_allele_v;
            for (const auto& p : next_anchor_support)
            {
                auto support_accum = [&] (unsigned s, Read_Entry_CBPtr re_cbptr) {
                    return s + (re_cbptr->name().substr(0, 5) == "merge"? _merged_weight : 1);
                };
                // each previously merged read counts as 5 when computing discordance
                unsigned allele_support_size = accumulate(p.second[0], 0, support_accum);
                allele_support_size += accumulate(p.second[1], 0, support_accum);
                ASSERT(allele_support_size > 0);
                if (allele_support_size > _max_discordant_support)
                {
                    next_allele_v.insert(p.first);
                }
            }
            if (next_allele_v.empty())
            {
                LOG("Read_Merger", info) << ptree("loop_end__no_allele_with_strong_support");
                allele_support[0].clear();
                allele_support[1].clear();
                return;
            }
            else if (next_allele_v.size() > 1)
            {
                LOG("Read_Merger", info) << ptree("loop_end__multiple_alleles_with_strong_support")
                    .put("next_anchor", next_anchor)
                    .put("next_allele_v", next_allele_v);
                allele_support[0].clear();
                allele_support[1].clear();
                return;
            }
            next_allele = *next_allele_v.begin();
            // support at next_anchor is ambiguous,
            // but a single allele has support greater than the discordant threshold
            LOG("Read_Merger", info) << ptree("loop")
                .put("next_allele", next_allele);
            for (auto& p : next_anchor_support)
            {
                if (p.first == next_allele) continue;
                // we split the remaining reads in between anchor and next_anchor
                for (int dir = 0; dir < 2; ++dir)
                {
                    for (auto re_cbptr : p.second[dir])
                    {
                        ASSERT(allele_support[(dir + init_c_direction) % 2].count(re_cbptr)
                               and active_allele_support[(dir + init_c_direction) % 2].count(re_cbptr));
                               and not inactive_allele_support[(dir + init_c_direction) % 2].count(re_cbptr));
                        allele_support[(dir + init_c_direction) % 2].erase(re_cbptr);
                        active_allele_support[(dir + init_c_direction) % 2].erase(re_cbptr);
                        auto rp = _g.split_read(
                            re_cbptr,
                            not crt_c_direction? crt_anchor : next_anchor);
                        Read_Entry_CBPtr new_re_cbptr = (not dir? rp.first : rp.second);
                        // we keep the side of the read that spans crt_anchor
                        auto tmp = crt_anchor.read_support().at(crt_allele);
                        tmp.reverse(crt_c_direction);
                        //TODO: remove
                        ASSERT(tmp[dir].count(new_re_cbptr) > 0);
                        allele_support[(dir + init_c_direction) % 2].insert(new_re_cbptr);
                        active_allele_support[(dir + init_c_direction) % 2].insert(new_re_cbptr);
                    }
                    p.second[dir].clear();
                }
            }
            filter_cont(
                next_anchor_support,
                [] (Anchor_Read_Support::const_iterator cit) {
                    return not cit->second.empty();
                });
        } // if (next_anchor_support.size() > 1)
        ASSERT(next_anchor_support.size() == 1);
        // reads which do not support next anchor become inactive
        active_allele_support.subtract(next_anchor_support.begin()->second, init_c_direction);
        LOG("Read_Merger", debug) << ptree("loop_filter")
            .put("new_inactive_allele_support", active_allele_support);
        inactive_allele_support.add(active_allele_support, false);
        crt_anchor = next_anchor;
        crt_allele = next_anchor_support.begin()->first;
    }
} // Read_Merger::extend_haploid_support

void Read_Merger::merge_reads(Contig_Entry_BPtr ce_bptr, const RE_DSet& re_dset) const
{
    ASSERT(re_dset.size() >= 2);
    LOG("Read_Merger", info) << ptree("begin")
        .put("ce_bptr", ce_bptr.to_int())
        .put("re_dset[0]", re_dset[0])
        .put("re_dset[1]", re_dset[1]);

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
        //RC_DSet crt_chunks = get_oriented_chunks(ce_bptr, re_dset);
        RC_DSet crt_chunks = move(*re_dset.chunks(ce_bptr));
        ASSERT(crt_chunks[0].size() == re_dset[0].size()
               and crt_chunks[1].size() == re_dset[1].size());
        if (dir == 0)
        {
            // merge chunks in start contig
            m_rc_bptr = merge_contig_chunks(crt_chunks, m_re_bptr);
            m_chunk_inserter(m_rc_bptr);
            m_ce_set.insert(m_rc_bptr->ce_bptr());
        }
        while (true)
        {
            LOG("Read_Merger", debug) << ptree("loop__start")
                .put("dir", dir)
                .put("crt_chunks", crt_chunks);
            // advance chunks but save unmappable chunks
            RC_DSet unmappable_chunks;
            tie(crt_chunks, unmappable_chunks) = advance_chunks(crt_chunks, dir);
            LOG("Read_Merger", debug) << ptree("loop__after_advance")
                .put("crt_chunks", crt_chunks)
                .put("unmappable_chunks", unmappable_chunks);
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
                    RE_DSet unmappable_dset;
                    for (int d = 0; d < 2; ++d)
                    {
                        for (auto rc_bptr : unmappable_chunks[d])
                        {
                            unmappable_dset[d].insert(rc_bptr->re_bptr());
                        }
                    }
                    // compute reads without unmappable region
                    RE_DSet non_unmappable_dset;
                    for (int d = 0; d < 2; ++d)
                    {
                        for (auto rc_bptr : crt_chunks[d])
                        {
                            if (unmappable_dset[d].count(rc_bptr->re_bptr()) == 0)
                            {
                                non_unmappable_dset[d].insert(rc_bptr->re_bptr());
                            }
                        }
                    }
                    LOG("Read_Merger", warning) << ptree("unmappable_inconsistency")
                        .put("prev_ce_cbptr", prev_ce_cbptr.to_int())
                        .put("next_ce_cbptr", next_ce_cbptr.to_int())
                        .put("unmappable_dset[0]", unmappable_dset[0])
                        .put("unmappable_dset[1]", unmappable_dset[1])
                        .put("non_unmappable_dset[0]", non_unmappable_dset[0])
                        .put("non_unmappable_dset[1]", non_unmappable_dset[1]);
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
    m_re_bptr->len() = pos;
    _g.re_cont().insert(m_re_bptr);
    _g.check(set< Read_Entry_CBPtr >{ m_re_bptr });

    // finally, destroy the old reads
    for (int d = 0; d < 2; ++d)
    {
        for (auto re_cbptr : re_dset.at(d))
        {
            _g.remove_read(re_cbptr.unconst());
        }
    }
    _g.check(set< Read_Entry_CBPtr >{ m_re_bptr });

    LOG("Read_Merger", info) << ptree("end");
}

Read_Chunk_BPtr Read_Merger::merge_contig_chunks(const RC_DSet& rc_dset, Read_Entry_BPtr m_re_bptr) const
{
    LOG("Read_Merger", debug) << ptree("begin")
        .put("rc_dset", rc_dset);
    ASSERT(not (rc_dset[0].empty() and rc_dset[1].empty()));
    auto ce_bptr = not rc_dset[0].empty()
        ? (*rc_dset[0].begin())->ce_bptr()
        : (*rc_dset[1].begin())->ce_bptr();
    auto c_direction = not rc_dset[0].empty()
        ? (*rc_dset[0].begin())->get_rc()
        : not (*rc_dset[1].begin())->get_rc();
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
        LOG("Read_Merger", debug) << ptree("before_increment")
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
    m_rc_bptr->set_unbreakable(true);
    m_ce_bptr->chunk_cont().insert(m_rc_bptr);
    LOG("Read_Merger", debug) << ptree("end")
        .put("unmappable_chunks", seq_v)
        .put("result", m_ce_bptr->seq());
    return m_rc_bptr;
}

pair< RC_DSet, RC_DSet >
Read_Merger::advance_chunks(const RC_DSet& crt_rc_dset, bool direction) const
{
    RC_DSet next_rc_dset;
    RC_DSet unmappable_rc_dset;
    for (int d = 0; d < 2; ++d)
    {
        for (auto rc_cbptr : crt_rc_dset[d])
        {
            ASSERT(rc_cbptr);
            auto next_rc_cbptr = rc_cbptr->re_bptr()->chunk_cont().get_sibling(
                rc_cbptr, true, not ( (d + direction) % 2 ));
            if (not next_rc_cbptr) continue;
            if (next_rc_cbptr->ce_bptr()->is_unmappable())
            {
                // skip one unmappable chunk if necessary
                unmappable_rc_dset[d].insert(next_rc_cbptr);
                auto next_next_rc_cbptr = next_rc_cbptr->re_bptr()->chunk_cont().get_sibling(
                    next_rc_cbptr, true, not ( (d + direction) % 2 ));
                if (not next_next_rc_cbptr)
                {
                    LOG("Read_Merger", warning) << ptree("terminal_unmappable_chunk")
                        .put("re_name", rc_cbptr->re_bptr()->name());
                    continue;
                }
                next_rc_cbptr = next_next_rc_cbptr;
            }
            next_rc_dset[d].insert(next_rc_cbptr);
        }
    }
    return make_pair(next_rc_dset, unmappable_rc_dset);
}

} // namespace MAC
