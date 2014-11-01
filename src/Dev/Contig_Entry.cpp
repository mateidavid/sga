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

map< std::tuple< Contig_Entry_CBPtr, bool >, vector< Read_Chunk_CBPtr > >
Contig_Entry::out_chunks_dir(bool c_right, int unmappable_policy, size_t ignore_threshold) const
{
    map< std::tuple< Contig_Entry_CBPtr, bool >, vector< Read_Chunk_CBPtr > > res;
    bool c_left = not c_right;
    Size_Type endpoint = (c_left? 0 : len());
    for (auto rc_cbref : chunk_cont().iintersect(endpoint, endpoint))
    {
        Read_Chunk_CBPtr rc_cbptr = &rc_cbref;
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
                    res[std::make_tuple(Contig_Entry_CBPtr(nullptr), false)].push_back(rc_cbptr);
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
                res[std::make_tuple(ce_next_cbptr, same_orientation)].push_back(rc_cbptr);
                break;
            case 2:
                res[std::make_tuple(Contig_Entry_CBPtr(nullptr), false)].push_back(rc_cbptr);
                break;
            default:
                ASSERT(false);
            }
        }
        else
        {
            res[std::make_tuple(ce_next_cbptr, same_orientation)].push_back(rc_cbptr);
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

tuple< Size_Type, Size_Type >
Contig_Entry::unmappable_neighbour_range(bool c_right, const vector< MAC::Read_Chunk_CBPtr >& chunk_cont) const
{
    if (chunk_cont.empty())
    {
        return std::make_tuple(Size_Type(0), Size_Type(0));
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
    return std::make_tuple(min_skip_len, max_skip_len);
}

auto Contig_Entry::neighbours(bool forward, Neighbour_Options opts, size_t ignore_threshold) -> neighbours_type
{
    neighbours_type res;
    Size_Type endpoint = (forward? len() : 0);
    for (auto rc_cbref : chunk_cont().iintersect(endpoint, endpoint))
    {
        Read_Chunk_CBPtr rc_cbptr = &rc_cbref;
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

tuple< Contig_Entry_CBPtr, bool, vector< Read_Chunk_CBPtr > >
Contig_Entry::can_cat_dir(bool c_right) const
{
    ASSERT(not chunk_cont().empty());
    Contig_Entry_CBPtr ce_cbptr = chunk_cont().begin()->ce_bptr();
    ASSERT(ce_cbptr.raw() == this);

    logger("Contig_Entry", debug1) << ptree("can_cat_dir")
        .put("ce_ptr", ce_cbptr.to_int())
        .put("c_right", c_right);

    auto res = out_chunks_dir(c_right, 1);
    logger("Contig_Entry", debug1) << ptree("can_cat_dir")
        .put("res_size", res.size());
    if (res.size() != 1)
    {
        return make_tuple(Contig_Entry_CBPtr(nullptr), false, vector< Read_Chunk_CBPtr >());
    }
    Contig_Entry_CBPtr ce_next_cbptr;
    bool same_orientation;
    tie(ce_next_cbptr, same_orientation) = res.begin()->first;
    vector< Read_Chunk_CBPtr > rc_cbptr_cont = move(res.begin()->second);
    ASSERT(not rc_cbptr_cont.empty());
    logger("Contig_Entry", debug1) << ptree("can_cat_dir")
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
    logger("Contig_Entry", debug1) << ptree("can_cat_dir")
        .put("tmp_size", tmp.size());
    if (tmp.size() != 1)
    {
        return make_tuple(Contig_Entry_CBPtr(nullptr), false, vector< Read_Chunk_CBPtr >());
    }
    ASSERT(tmp.begin()->first == make_tuple(ce_cbptr, same_orientation));
    ASSERT(tmp.begin()->second.size() == rc_cbptr_cont.size());
    return make_tuple(ce_next_cbptr, same_orientation, move(rc_cbptr_cont));
}

void Contig_Entry::cat_c_right(Contig_Entry_BPtr ce_bptr, Contig_Entry_BPtr ce_next_bptr,
                               vector< MAC::Read_Chunk_CBPtr >& rc_cbptr_cont)
{
    logger("Contig_Entry", debug1) << ptree("cat_c_right")
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

/*
vector< Mutation_CPtr > Contig_Entry::get_separated_het_mutations(
    size_t min_support_report, Size_Type min_separation) const
{
    vector< Mutation_CPtr > res;
    Size_Type last_mut_end = 0;
    for (auto mut_it = _mut_cont.begin(); mut_it != _mut_cont.end(); ++mut_it)
    {
        Mutation_CPtr mut_cptr = &*mut_it;
        auto next_mut_it = mut_it;
        ++next_mut_it;
        if (mut_cptr->rf_start() >= last_mut_end + min_separation
                and mut_cptr->rf_end() + min_separation <= len()
                and (next_mut_it == _mut_cont.end()
                     or mut_cptr->rf_end() + min_separation <= next_mut_it->rf_start()))
        {
            // well separated
            size_t n_chunks_supporting_mut[2] = { 0, 0 };
            for (auto& rc_cptr : _chunk_cptr_cont)
            {
                if (rc_cptr->get_c_start() > mut_cptr->rf_start()
                        or rc_cptr->get_c_end() < mut_cptr->rf_end())
                {
                    // doesn't span mutation
                    continue;
                }
                ++n_chunks_supporting_mut[int(rc_cptr->have_mutation(mut_cptr))];
            }
            if (n_chunks_supporting_mut[0] >= min_support_report
                    and n_chunks_supporting_mut[1] >= min_support_report)
            {
                res.push_back(mut_cptr);
            }
        }
        last_mut_end = std::max(last_mut_end, mut_cptr->rf_end());
    }
    return res;
}

void Contig_Entry::print_separated_het_mutations(
    ostream& os, size_t min_support_report, Size_Type min_separation) const
{
    vector< Mutation_CPtr > v = get_separated_het_mutations(min_support_report, min_separation);
    for (const auto& mut_cptr : v)
    {
        if (mut_cptr->is_del())
        {
            os << _seq_ptr->substr(mut_cptr->rf_start(), mut_cptr->rf_len()) << "\t\t";
        }
        else
        {
            os << mut_cptr->get_seq() << '\t' << _seq_ptr->substr(mut_cptr->rf_start(), mut_cptr->rf_len()) << '\t';
        }
        os << _seq_ptr->substr(mut_cptr->rf_start() - min_separation, min_separation) << '\t'
           << _seq_ptr->substr(mut_cptr->rf_end(), min_separation) << '\n';
    }
}
*/

void Contig_Entry::check() const
{
#ifndef BOOST_DISABLE_ASSERTS
    // there are chunks mapped to this contig
    ASSERT(not chunk_cont().empty());
    if (is_unmappable())
    {
        ASSERT(mut_cont().empty());
        ASSERT(chunk_cont().single_node());
        ASSERT(chunk_cont().begin()->is_unbreakable());
    }
    // check mutation container
    mut_cont().check();
    // mutations:
    for (auto mut_cbref : mut_cont())
    {
        Mutation_CBPtr mut_cbptr = &mut_cbref;
        // base coordinates
        ASSERT(mut_cbptr->rf_start() <= seq().size() and mut_cbptr->rf_end() <= seq().size());
        // no empty mutations
        ASSERT(not mut_cbptr->is_empty());
        // read chunk ptr container
        mut_cbptr->chunk_ptr_cont().check();
        // mutations must be observed by chunks
        ASSERT(not mut_cbptr->chunk_ptr_cont().empty());
        // check read_chunk_ptr_cont
        for (auto mca_cbref : mut_cbptr->chunk_ptr_cont())
        {
            Mutation_Chunk_Adapter_CBPtr mca_cbptr = &mca_cbref;
            // Mutation back pointers
            ASSERT(mca_cbptr->mut_cbptr() == mut_cbptr);
            // Read_Chunk part of current contig
            ASSERT(mca_cbptr->chunk_cbptr()->ce_bptr().raw() == this);
        }
    }
    // check chunk container
    chunk_cont().check();
    // read chunks:
    for (auto rc_cbref : chunk_cont())
    {
        Read_Chunk_CBPtr rc_cbptr = &rc_cbref;
        // contig entry pointers
        ASSERT(rc_cbptr->ce_bptr() and rc_cbptr->ce_bptr().raw() == this);
        // contig coordinates
        ASSERT(rc_cbptr->get_c_end() <= seq().size());
        // mutation pointers:
        for (auto mca_cbref : rc_cbptr->mut_ptr_cont())
        {
            Mutation_Chunk_Adapter_CBPtr mca_cbptr = &mca_cbref;
            // Read_Chunk back pointers
            ASSERT(mca_cbptr->chunk_cbptr() == rc_cbptr);
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
