//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Contig_Entry.hpp"
#include "Read_Entry.hpp"


namespace MAC
{

void Contig_Entry::cut_mutation(Mutation_BPtr mut_bptr, Size_Type c_offset, Size_Type r_offset)
{
    // Mutation must be in this container
    ASSERT(mut_cont().find(mut_bptr, true));

    // unlink the Mutation from its container
    mut_cont().erase(mut_bptr);

    // cut mutation, save second part
    Mutation_BPtr new_mut_bptr = mut_bptr->cut(c_offset, r_offset).unconst();

    // add both new and old to container
    mut_cont().insert(mut_bptr);
    mut_cont().insert(new_mut_bptr);

    // clone Read_Chunk_Ptr container from original Mutation
    ASSERT(new_mut_bptr->chunk_ptr_cont().empty());
    new_mut_bptr->chunk_ptr_cont().clone_from(mut_bptr->chunk_ptr_cont(), new_mut_bptr);

    // insert new Mutation in the Mutation_Ptr containers of the affected Read_Chunk objects
    {
        auto mca_it = mut_bptr->chunk_ptr_cont().begin();
        auto new_mca_it = new_mut_bptr->chunk_ptr_cont().begin();
        while (mca_it != mut_bptr->chunk_ptr_cont().end())
        {
            ASSERT(new_mca_it != new_mut_bptr->chunk_ptr_cont().end());
            Read_Chunk_BPtr rc_bptr = mca_it->chunk_cbptr().unconst();
            Mutation_Chunk_Adapter_BPtr new_mca_bptr = &*new_mca_it;

            rc_bptr->mut_ptr_cont().insert_after(rc_bptr->mut_ptr_cont().iterator_to(*mca_it), new_mca_bptr);

            ++mca_it;
            ++new_mca_it;
        }
        ASSERT(new_mca_it == new_mut_bptr->chunk_ptr_cont().end());
    }

    check();
}

Contig_Entry_BPtr Contig_Entry::cut(Size_Type c_brk, Mutation_CBPtr mut_left_cbptr, bool strict)
{
    ASSERT(not mut_left_cbptr
           or (mut_left_cbptr->rf_start() == c_brk and mut_left_cbptr->is_ins()));

    LOG("Contig_Entry", debug1) << ptree("cut")
        .put("c_brk", c_brk)
        .put("mut_left_ptr", mut_left_cbptr.to_ptree())
        .put("strict", strict);
    // there is nothing to cut if:
    if (not strict
        and (// cut is at the start, and no insertion goes to the left
            (c_brk == 0 and not mut_left_cbptr)
            // cut is at the end, and:
            or (c_brk == len()
                and (// there are no mutations
                    mut_cont().empty()
                    // or the last one is not an insertion
                    or not mut_cont().rbegin()->is_ins()
                    // or the last insertion is not at the end
                    or mut_cont().rbegin()->rf_start() < len()))))
    {
        return nullptr;
    }
    // split any mutations that span c_pos
    while (true)
    {
        Mutation_BPtr mut_bptr = mut_cont().find_span_pos(c_brk).unconst();
        if (not mut_bptr)
        {
            break;
        }
        ASSERT(mut_bptr->rf_start() < c_brk and c_brk < mut_bptr->rf_end());
        ASSERT(mut_bptr != mut_left_cbptr);
        Size_Type mut_rf_brk = c_brk - mut_bptr->rf_start();
        cut_mutation(mut_bptr, mut_rf_brk, min(mut_rf_brk, mut_bptr->seq_len()));
    }
    // create new contig entry object; set base sequences
    Contig_Entry_BPtr ce_new_bptr = Contig_Entry_Fact::new_elem(Seq_Type(seq().substr(c_brk)));
    seq().resize(c_brk);
    // split Mutation_Cont, save rhs in new Contig_Entry
    ce_new_bptr->mut_cont().splice(mut_cont(), c_brk, mut_left_cbptr);
    // unlink Read_Chunk objects from their RE containers
    auto rc_next_map = chunk_cont().erase_from_re_cont();
    // split Read_Chunk_Cont, save rhs in new Contig_Entry
    ASSERT(ce_new_bptr->chunk_cont().empty());
    auto rc_split_map = ce_new_bptr->chunk_cont().splice(chunk_cont(), c_brk, mut_left_cbptr, strict);
    ASSERT(not chunk_cont().empty());
    ASSERT(chunk_cont().max_end() <= c_brk);
    ASSERT(ce_new_bptr->chunk_cont().empty() or c_brk <= ce_new_bptr->chunk_cont().begin()->get_c_start());
    // rebase all mutations and read chunks from the rhs to the breakpoint
    ce_new_bptr->mut_cont().shift(-int(c_brk));
    ce_new_bptr->chunk_cont().shift(-int(c_brk));
    ce_new_bptr->chunk_cont().set_ce_ptr(ce_new_bptr);
    // fix rc_next_map to incorporate the split information
    // 1. fix second components
    for (auto p : rc_next_map)
    {
        Read_Chunk_CBPtr rc_cbptr = p.first;
        Read_Chunk_CBPtr rc_next_cbptr = p.second;
        if (rc_split_map.count(rc_next_cbptr))
        {
            // can only happen if next chunk is also involved in the split
            ASSERT(rc_next_map.count(rc_next_cbptr));
            // second component changes iff rc_next is mapped to the rc of its ce
            if (rc_next_cbptr->get_rc())
            {
                rc_next_map[rc_cbptr] = rc_split_map[rc_next_cbptr];
            }
        }
    }
    // 2. fix first components
    {
        list< decltype(rc_next_map)::value_type > aux_list;
        for (auto p : rc_next_map)
        {
            Read_Chunk_CBPtr rc_cbptr = p.first;
            Read_Chunk_CBPtr rc_next_cbptr = p.second;
            if (rc_split_map.count(rc_cbptr))
            {
                Read_Chunk_CBPtr rc_new_cbptr = rc_split_map[rc_cbptr];
                ASSERT(rc_next_map.count(rc_new_cbptr) == 0);
                if (not rc_cbptr->get_rc())
                {
                    aux_list.emplace_back(rc_new_cbptr, rc_next_cbptr);
                    rc_next_map[rc_cbptr] = rc_new_cbptr;
                }
                else
                {
                    aux_list.emplace_back(rc_new_cbptr, rc_cbptr);
                }
            }
        }
        rc_next_map.insert(aux_list.begin(), aux_list.end());
    }
    // link back the chunks into their RE containers
    Read_Chunk_CE_Cont::insert_into_re_cont(rc_next_map);
    // remove unused Mutation objects
    mut_cont().drop_unused();
    ce_new_bptr->mut_cont().drop_unused();
    check();
    //TODO: remove
    for (auto rc_bptr : chunk_cont() | referenced)
    {
        rc_bptr->re_bptr()->check();
    }
    if (ce_new_bptr->chunk_cont().empty())
    {
        Contig_Entry_Fact::del_elem(ce_new_bptr);
        return nullptr;
    }
    else
    {
        ce_new_bptr->check();
        return ce_new_bptr;
    }
} // cut

tuple< set< Read_Chunk_CBPtr >, set< Read_Chunk_CBPtr >, set< Read_Chunk_CBPtr > >
Contig_Entry::mut_support(Mutation_CBPtr mut_cbptr) const
{
    set< Read_Chunk_CBPtr > qr_set;
    set< Read_Chunk_CBPtr > full_rf_set;
    set< Read_Chunk_CBPtr > part_rf_set;
    // first, get support for qr allele
    for (auto mca_cbptr : mut_cbptr->chunk_ptr_cont() | referenced)
    {
        qr_set.insert(mca_cbptr->chunk_cbptr());
    }
    // next, get support for rf allele
    // add all chunks fully spanning mutation
    auto iint_res = chunk_cont().iintersect(mut_cbptr->rf_start(), mut_cbptr->rf_end());
    for (auto rc_cbptr : iint_res | referenced)
    {
        // if chunk observes qr allele, skip it
        if (qr_set.count(rc_cbptr))
        {
            continue;
        }
        // if the rf allele is empty, we require >=1bp mapped on either side
        if ((mut_cbptr->rf_len() > 0
             and rc_cbptr->get_c_start() <= mut_cbptr->rf_start()
             and rc_cbptr->get_c_end() >= mut_cbptr->rf_end())
            or (mut_cbptr->rf_len() == 0
                and rc_cbptr->get_c_start() < mut_cbptr->rf_start()
                and rc_cbptr->get_c_end() > mut_cbptr->rf_end()))
        {
            full_rf_set.insert(rc_cbptr);
        }
        else
        {
            part_rf_set.insert(rc_cbptr);
        }
    }
    return make_tuple(qr_set, full_rf_set, part_rf_set);
}

void Contig_Entry::reverse()
{
    // only reverse full contigs
    ASSERT(_seq_offset == 0);
    // reverse the base sequence
    _seq = _seq.revcomp();
    // reverse Mutation objects in their container
    mut_cont().reverse_mutations(len());
    chunk_cont().reverse();
}

map< pair< Contig_Entry_CBPtr, bool >, set< Read_Chunk_CBPtr > >
Contig_Entry::out_chunks_dir(bool c_right, int unmappable_policy, size_t ignore_threshold) const
{
    map< pair< Contig_Entry_CBPtr, bool >, set< Read_Chunk_CBPtr > > res;
    bool c_left = not c_right;
    Size_Type endpoint = (c_left? 0 : len());
    for (auto rc_cbptr : chunk_cont().iintersect(endpoint, endpoint) | referenced)
    {
        Read_Entry_BPtr re_cbptr = rc_cbptr->re_bptr();
        // we might be skipping chunks mapped in different ways; we need read direction to be consistent
        bool r_right = (c_right != rc_cbptr->get_rc());
        Read_Chunk_CBPtr rc_next_cbptr = re_cbptr->chunk_cont().get_sibling(rc_cbptr, true, r_right);
        if (not rc_next_cbptr)
        {
            // rc is terminal, nothing to do
            continue;
        }
        if (unmappable_policy == 3 or unmappable_policy == 4)
        {
            // skip unmappable chunks
            while (rc_next_cbptr and rc_next_cbptr->ce_bptr()->is_unmappable())
            {
                rc_next_cbptr = re_cbptr->chunk_cont().get_sibling(rc_next_cbptr, true, r_right);
            }
            if (not rc_next_cbptr)
            {
                if ( unmappable_policy == 4)
                {
                    res[make_pair(Contig_Entry_CBPtr(nullptr), false)].insert(rc_cbptr);
                }
                continue;
            }
        }
        ASSERT(rc_next_cbptr);
        Contig_Entry_CBPtr ce_next_cbptr = rc_next_cbptr->ce_bptr();
        bool same_orientation = (rc_cbptr->get_rc() == rc_next_cbptr->get_rc());
        if (rc_next_cbptr->ce_bptr()->is_unmappable())
        {
            ASSERT(unmappable_policy != 3 and unmappable_policy != 4);
            switch (unmappable_policy)
            {
            case 0:
                break;
            case 1:
                res[make_pair(ce_next_cbptr, same_orientation)].insert(rc_cbptr);
                break;
            case 2:
                res[make_pair(Contig_Entry_CBPtr(nullptr), false)].insert(rc_cbptr);
                break;
            default:
                ASSERT(false);
            }
        }
        else
        {
            res[make_pair(ce_next_cbptr, same_orientation)].insert(rc_cbptr);
        }
    }
    // remove (contig,orientation) pairs with low support
    auto it = res.begin();
    while (it != res.end())
    {
        if (it->second.size() <= ignore_threshold)
        {
            auto it_next = next(it);
            res.erase(it);
            it = it_next;
        }
        else
        {
            ++it;
        }
    }
    return res;
}

pair< Size_Type, Size_Type >
Contig_Entry::unmappable_neighbour_range(bool c_right, const set< MAC::Read_Chunk_CBPtr >& chunk_cont) const
{
    if (chunk_cont.empty())
    {
        return make_pair(Size_Type(0), Size_Type(0));
    }
    Size_Type min_skip_len = numeric_limits< Size_Type >::max();
    Size_Type max_skip_len = 0;
    for (auto rc_cbptr : chunk_cont)
    {
        ASSERT(rc_cbptr);
        Read_Entry_CBPtr re_cbptr = rc_cbptr->re_bptr();
        ASSERT(c_right or rc_cbptr->get_c_start() == 0);
        ASSERT(not c_right or rc_cbptr->get_c_end() == len());
        Read_Chunk_CBPtr rc_next_cbptr = re_cbptr->chunk_cont().get_sibling(rc_cbptr, false, c_right);
        ASSERT(rc_next_cbptr);
        Size_Type skip_len = 0;
        while (rc_next_cbptr->ce_bptr()->is_unmappable())
        {
            skip_len += rc_next_cbptr->get_r_len();
            rc_next_cbptr = re_cbptr->chunk_cont().get_sibling(rc_next_cbptr, false, c_right);
            ASSERT(rc_next_cbptr);
        }
        min_skip_len = min(min_skip_len, skip_len);
        max_skip_len = max(max_skip_len, skip_len);
    }
    return make_pair(min_skip_len, max_skip_len);
}

#if 0
auto Contig_Entry::neighbours(bool forward, Neighbour_Options opts, size_t ignore_threshold) -> neighbours_type
{
    neighbours_type res;
    Size_Type endpoint = (forward? len() : 0);
    for (auto rc_cbptr : chunk_cont().iintersect(endpoint, endpoint) | referenced)
    {
        Read_Entry_CBPtr re_cbptr = rc_cbptr->re_bptr();
        // we might be skipping chunks mapped in different ways; we need read direction for consistent traversal
        bool r_forward = (forward != rc_cbptr->get_rc());
        Read_Chunk_CBPtr rc_next_cbptr = re_cbptr->chunk_cont().get_sibling(rc_cbptr, true, r_forward);
        if (not rc_next_cbptr)
        {
            // rc is terminal, nothing to do
            continue;
        }
        Contig_Entry_CBPtr ce_next_cbptr = rc_next_cbptr->ce_bptr();
        // skip chunks as dictated by the options
        while (opts.get_opt(ce_next_cbptr->category()) == Neighbour_Options::opt_skip
               or opts.get_opt(ce_next_cbptr->category()) == Neighbour_Options::opt_skip_not_last)
        {
            Read_Chunk_CBPtr rc_next_next_cbptr = re_cbptr->chunk_cont().get_sibling(rc_next_cbptr, true, r_forward);
            if (not rc_next_next_cbptr)
            {
                break;
            }
            rc_next_cbptr = rc_next_next_cbptr;
            ce_next_cbptr = rc_next_cbptr->ce_bptr();
        }
        if (opts.get_opt(ce_next_cbptr->category()) == Neighbour_Options::opt_drop
            or opts.get_opt(ce_next_cbptr->category()) == Neighbour_Options::opt_skip)
        {
            continue;
        }
        ASSERT(opts.get_opt(ce_next_cbptr->category()) == Neighbour_Options::opt_none
               or (opts.get_opt(ce_next_cbptr->category()) == Neighbour_Options::opt_skip_not_last
                   and not re_cbptr->chunk_cont().get_sibling(rc_next_cbptr, true, r_forward)));
        bool same_orientation = (rc_cbptr->get_rc() == rc_next_cbptr->get_rc());
        res[make_tuple(ce_next_cbptr, same_orientation)].push_back(make_tuple(rc_cbptr, rc_next_cbptr));
    }
    // remove (contig,orientation) pairs with low support
    if (ignore_threshold > 0)
    {
        for (auto it = res.begin(); it != res.end(); )
        {
            if (it->second.size() <= ignore_threshold)
            {
                res.erase(it++);
            }
            else
            {
                ++it;
            }
        }
    }
    return res;
}
#endif

tuple< Contig_Entry_CBPtr, bool, set< Read_Chunk_CBPtr > >
Contig_Entry::can_cat_dir(bool c_right) const
{
    ASSERT(not chunk_cont().empty());
    Contig_Entry_CBPtr ce_cbptr = chunk_cont().begin()->ce_bptr();
    ASSERT(ce_cbptr.raw() == this);

    LOG("Contig_Entry", debug1) << ptree("can_cat_dir")
        .put("ce_ptr", ce_cbptr.to_int())
        .put("c_right", c_right);

    auto res = out_chunks_dir(c_right, 1);
    LOG("Contig_Entry", debug1) << ptree("can_cat_dir")
        .put("res_size", res.size());
    if (res.size() != 1)
    {
        return make_tuple(Contig_Entry_CBPtr(nullptr), false, set< Read_Chunk_CBPtr >());
    }
    Contig_Entry_CBPtr ce_next_cbptr;
    bool same_orientation;
    tie(ce_next_cbptr, same_orientation) = res.begin()->first;
    set< Read_Chunk_CBPtr > rc_cbptr_cont = move(res.begin()->second);
    ASSERT(not rc_cbptr_cont.empty());
    LOG("Contig_Entry", debug1) << ptree("can_cat_dir")
        .put("ce_next_ptr", ce_next_cbptr.to_int())
        .put("same_orientation", same_orientation)
        .put("chunk_cont_size", rc_cbptr_cont.size());
    //if (ce_next_cbptr->is_unmappable() != is_unmappable())
    /*
    if (ce_next_cbptr->category() != category())
    {
        return std::make_tuple(Contig_Entry_CBPtr(nullptr), false, vector< Read_Chunk_CBPtr >());
    }
    */
    auto tmp = ce_next_cbptr->out_chunks_dir(c_right != same_orientation, 1);
    LOG("Contig_Entry", debug1) << ptree("can_cat_dir")
        .put("tmp_size", tmp.size());
    if (tmp.size() != 1)
    {
        return make_tuple(Contig_Entry_CBPtr(nullptr), false, set< Read_Chunk_CBPtr >());
    }
    ASSERT(tmp.begin()->first == make_pair(ce_cbptr, same_orientation));
    ASSERT(tmp.begin()->second.size() == rc_cbptr_cont.size());
    return make_tuple(ce_next_cbptr, same_orientation, move(rc_cbptr_cont));
}

void Contig_Entry::cat_c_right(Contig_Entry_BPtr ce_bptr, Contig_Entry_BPtr ce_next_bptr,
                               set< Read_Chunk_CBPtr >& rc_cbptr_cont)
{
    LOG("Contig_Entry", debug1) << ptree("cat_c_right")
        .put("ce_ptr", ce_bptr.to_int())
        .put("ce_next_ptr", ce_next_bptr.to_int());

    ASSERT(ce_bptr->_seq_offset == 0 and ce_next_bptr->_seq_offset == 0);
    ASSERT(ce_next_bptr->is_unlinked());
    // first, shift all mutations and chunks from the second Contig_Entry
    ce_next_bptr->mut_cont().shift(int(ce_bptr->len()));
    ce_next_bptr->chunk_cont().shift(int(ce_bptr->len()));
    // move all chunks and mutations into first Contig_Entry
    ce_bptr->mut_cont().splice(ce_next_bptr->mut_cont());
    ce_bptr->chunk_cont().splice(ce_next_bptr->chunk_cont(), ce_bptr);
    // grab next contig's sequence
    ce_bptr->seq() += ce_next_bptr->seq();
    // merge read chunks previously spanning the break
    for (auto rc_cbptr : rc_cbptr_cont)
    {
        Read_Chunk_BPtr rc_bptr = rc_cbptr.unconst();
        ASSERT(rc_bptr->ce_bptr() == ce_bptr);
        Read_Chunk_BPtr rc_next_bptr = rc_bptr->re_bptr()->chunk_cont().get_sibling(rc_bptr, false, true).unconst();
        ASSERT(rc_next_bptr);
        ASSERT(rc_next_bptr->ce_bptr() == ce_bptr);
        ASSERT(rc_next_bptr->get_c_start() == rc_bptr->get_c_end());
        ASSERT(rc_bptr->get_rc() == rc_next_bptr->get_rc());
        // unlink chunks from RE&CE containers
        rc_bptr->re_bptr()->chunk_cont().erase(rc_bptr);
        rc_bptr->re_bptr()->chunk_cont().erase(rc_next_bptr);
        ce_bptr->chunk_cont().erase(rc_bptr);
        ce_bptr->chunk_cont().erase(rc_next_bptr);
        // cat left&right
        Read_Chunk::cat_c_right(rc_bptr, rc_next_bptr, ce_bptr->mut_cont());
        // insert resulting chunk back in its RE&CE containers
        rc_bptr->re_bptr()->chunk_cont().insert(rc_bptr);
        ce_bptr->chunk_cont().insert(rc_bptr);
    }
    // deallocate rhs contig
    Contig_Entry_Fact::del_elem(ce_next_bptr);
}

Mutation_CBPtr Contig_Entry::swap_mutation_alleles(Contig_Entry_BPtr ce_bptr, Mutation_CBPtr mut_cbptr)
{
    LOG("Contig_Entry", debug) << ptree("swap_mutation_alleles")
        .put("ce_ptr", ce_bptr.to_int())
        .put("mut_ptr", mut_cbptr.to_int());

    // save a read entry that has this mutation, and its coordinates
    ASSERT(not mut_cbptr->chunk_ptr_cont().empty());
    Read_Chunk_CBPtr rc_cbptr = mut_cbptr->chunk_ptr_cont().begin()->chunk_cbptr();
    Read_Entry_CBPtr re_cbptr = rc_cbptr->re_bptr();
    Range_Type c_rg(mut_cbptr->rf_start(), mut_cbptr->rf_end());
    // ask for maximal mapped read range iff mutation is an insertion
    bool maximal = mut_cbptr->is_ins();
    auto r_rg = rc_cbptr->mapped_range(c_rg, true, maximal, maximal);
    ASSERT((c_rg.len() == 0) == mut_cbptr->is_ins());
    ASSERT((r_rg.len() == 0) == mut_cbptr->is_del());

    // check there are no other mutations touching the endpoints of the contig range
    ASSERT(all_of(ce_bptr->mut_cont().begin(), ce_bptr->mut_cont().end(),
                  [&] (Mutation_CBRef mut_other_cbref) {
                      Mutation_CBPtr mut_other_cbptr = &mut_other_cbref;
                      return (mut_other_cbptr == mut_cbptr
                              or mut_other_cbptr->rf_end() < c_rg.start()
                              or mut_other_cbptr->rf_start() > c_rg.end());
                  }));

    // cut contig entry into 3 slices, so that mutation spans the full middle contig
    auto ce_mid_bptr = ce_bptr->cut(c_rg.start(), nullptr, true);
    ASSERT(ce_mid_bptr);
    rc_cbptr = re_cbptr->chunk_cont().get_chunk_with_pos(r_rg.start());
    ASSERT(rc_cbptr->ce_bptr() == ce_mid_bptr);
    ASSERT(rc_cbptr->get_c_start() == 0);
    ASSERT(not rc_cbptr->mut_ptr_cont().empty());
    mut_cbptr = rc_cbptr->mut_ptr_cont().begin()->mut_cbptr();
    ASSERT(mut_cbptr->rf_start() == 0);
    ASSERT(mut_cbptr->rf_len() == c_rg.len());
    auto ce_rhs_bptr = ce_mid_bptr->cut(c_rg.len(), mut_cbptr->is_ins()? mut_cbptr : nullptr, true);
    ASSERT(ce_rhs_bptr);
    ASSERT(ce_mid_bptr->len() == c_rg.len());
    ASSERT(size_one(ce_mid_bptr->mut_cont()));
    ASSERT(mut_cbptr == &*ce_mid_bptr->mut_cont().begin());
    ASSERT(mut_cbptr->rf_start() == 0);
    ASSERT(mut_cbptr->rf_len() == c_rg.len());
    ASSERT(rc_cbptr->ce_bptr() == ce_mid_bptr);
    ASSERT(rc_cbptr->get_c_start() == 0);
    ASSERT(rc_cbptr->get_c_len() == c_rg.len());
    ASSERT(rc_cbptr->get_r_len() == r_rg.len());
    ASSERT(size_one(rc_cbptr->mut_ptr_cont()));

    // create replacement contig entry
    Contig_Entry_BPtr ce_new_bptr = Contig_Entry_Fact::new_elem(Seq_Type(mut_cbptr->seq()));
    // create replacement mutation
    Mutation_BPtr mut_new_bptr = Mutation_Fact::new_elem(0, ce_new_bptr->len(), Seq_Type(ce_mid_bptr->seq()));
    ce_new_bptr->mut_cont().insert(mut_new_bptr);
    auto rc_map = ce_mid_bptr->chunk_cont().erase_from_re_cont();
    // create new chunks
    for (Read_Chunk_CBPtr rc_cbptr : ce_mid_bptr->chunk_cont() | referenced)
    {
        auto it = rc_map.find(rc_cbptr);
        ASSERT(it != rc_map.end());
        auto p = *it;
        ASSERT(p.first == rc_cbptr);
        rc_map.erase(it);
        // remove terminal chunks which end in the middle of a mutation
        if (rc_cbptr->mut_ptr_cont().empty()
            and (rc_cbptr->get_c_start() > 0 or rc_cbptr->get_c_end() < ce_mid_bptr->len()))
        {
            bool r_end = ((rc_cbptr->get_c_start() > 0 and rc_cbptr->get_rc())
                          or (rc_cbptr->get_c_end() < ce_mid_bptr->len() and not rc_cbptr->get_rc()));
            // check the chunk is terminal
            ASSERT(r_end or p.second == &*rc_cbptr->re_bptr()->chunk_cont().begin());
            ASSERT(not r_end or not p.second);
            // trim the Read_Entry
            rc_cbptr->re_bptr()->trim_end(r_end, rc_cbptr->get_seq());
            continue;
        }
        // create new read chunk
        Read_Chunk_BPtr rc_new_bptr = Read_Chunk_Fact::new_elem(
            rc_cbptr->get_r_start(), rc_cbptr->get_r_end(),
            0, ce_new_bptr->len(),
            rc_cbptr->get_rc());
        rc_new_bptr->re_bptr() = rc_cbptr->re_bptr();
        rc_new_bptr->ce_bptr() = ce_new_bptr;
        // if old chunk didn't have the mutation, the new one does
        if (rc_cbptr->mut_ptr_cont().empty())
        {
            Mutation_Chunk_Adapter_BPtr mca_new_bptr = Mutation_Chunk_Adapter_Fact::new_elem(
                mut_new_bptr, rc_new_bptr);
            rc_new_bptr->mut_ptr_cont().push_back(mca_new_bptr);
            mut_new_bptr->chunk_ptr_cont().insert(mca_new_bptr);
        }
        ce_new_bptr->chunk_cont().insert(rc_new_bptr);
        rc_map[rc_new_bptr] = p.second;
    }
    // drop mutation if no chunks use it
    ce_new_bptr->mut_cont().drop_unused();
    // insert new chunks into their re containers
    ce_new_bptr->chunk_cont().insert_into_re_cont(rc_map);
    ASSERT(rc_map.empty());
    // destroy old middle contig
    ce_mid_bptr->chunk_cont().clear_and_dispose();
    ce_mid_bptr->mut_cont().clear_and_dispose();
    // deallocate old middle contig
    Contig_Entry_Fact::del_elem(ce_mid_bptr);
    // concatenate contigs to reconstruct original contig entry
    {
        auto res = ce_new_bptr->can_cat_dir(true);
        ASSERT(get<0>(res) == ce_rhs_bptr);
        ASSERT(get<1>(res));
        Contig_Entry::cat_c_right(ce_new_bptr, ce_rhs_bptr, get<2>(res));
    }
    {
        auto res = ce_bptr->can_cat_dir(true);
        ASSERT(get<0>(res) == ce_new_bptr);
        ASSERT(get<1>(res));
        Contig_Entry::cat_c_right(ce_bptr, ce_new_bptr, get<2>(res));
    }
    auto rg = ce_bptr->mut_cont().iintersect(c_rg.start(), c_rg.start());
    ASSERT(rg.begin() != rg.end());
    ASSERT(++rg.begin() == rg.end());

    return &*rg.begin();
}

void Contig_Entry::check() const
{
#ifndef DISABLE_ASSERTS
    // there are chunks mapped to this contig
    ASSERT(not chunk_cont().empty());
    if (is_unmappable())
    {
        ASSERT(mut_cont().empty());
        ASSERT(size_one(chunk_cont()));
        ASSERT(chunk_cont().begin()->is_unbreakable());
    }
    // check mutation container
    mut_cont().check();
    // mutations:
    for (auto mut_cbptr : mut_cont() | referenced)
    {
        // base coordinates
        ASSERT(mut_cbptr->rf_start() <= seq().size() and mut_cbptr->rf_end() <= seq().size());
        // no empty mutations
        ASSERT(not mut_cbptr->is_empty());
        // read chunk ptr container
        mut_cbptr->chunk_ptr_cont().check();
        // mutations must be observed by chunks
        ASSERT(not mut_cbptr->chunk_ptr_cont().empty());
        // check read_chunk_ptr_cont
        for (auto mca_cbptr : mut_cbptr->chunk_ptr_cont() | referenced)
        {
            // Mutation back pointers
            ASSERT(mca_cbptr->mut_cbptr() == mut_cbptr);
            // Read_Chunk part of current contig
            ASSERT(mca_cbptr->chunk_cbptr()->ce_bptr().raw() == this);
        }
    }
    // check chunk container
    chunk_cont().check();
    // read chunks:
    for (auto rc_bptr : chunk_cont() | referenced)
    {
        // contig entry pointers
        ASSERT(rc_bptr->ce_bptr() and rc_bptr->ce_bptr().raw() == this);
        // contig coordinates
        ASSERT(rc_bptr->get_c_end() <= seq().size());
        // mutation pointers:
        for (auto mca_cbptr : rc_bptr->mut_ptr_cont() | referenced)
        {
            // Read_Chunk back pointers
            ASSERT(mca_cbptr->chunk_cbptr() == rc_bptr);
            // Mutation pointers point inside Mutation container
            ASSERT(mut_cont().find(mca_cbptr->mut_cbptr(), true));
        }
    }
#endif
}

boost::property_tree::ptree Contig_Entry::to_ptree() const
{
    return ptree().put("addr", (void*)this)
                  .put("seq", seq())
                  .put("len", len())
                  .put("tag", tag())
                  .put("is_unmappable", is_unmappable())
                  .put("mut_cont", cont_to_ptree(mut_cont()))
                  .put("chunk_cont", cont_to_ptree(chunk_cont()));
}


} // namespace MAC
