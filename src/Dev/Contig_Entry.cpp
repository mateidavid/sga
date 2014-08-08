//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Contig_Entry.hpp"
#include "Mutation.hpp"
#include "Read_Chunk.hpp"
#include "Read_Entry.hpp"
#include "print_seq.hpp"
#include "../Util/Util.h"

using namespace std;


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
    ASSERT(new_mut_bptr->chunk_ptr_cont().size() == 0);
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

    ASSERT(check());
}

void Contig_Entry::reverse()
{
    // only reverse full contigs
    ASSERT(_seq_offset == 0);
    // reverse the base sequence
    _seq = reverseComplement(_seq);
    // reverse Mutation objects in their container
    mut_cont().reverse_mutations(get_len());
    for (auto rc_bref : chunk_cont())
    {
        Read_Chunk_BPtr rc_bptr = &rc_bref;
        rc_bptr->mut_ptr_cont().reverse();
        rc_bptr->reverse();
    }
}

std::tuple< size_t, size_t, size_t, size_t >
Contig_Entry::get_out_degrees(int unmappable_policy, size_t ignore_threshold) const
{
    auto neighbours_left_cont = out_chunks_dir(false, unmappable_policy, ignore_threshold);
    auto neighbours_right_cont = out_chunks_dir(true, unmappable_policy, ignore_threshold);
    size_t total_left = 0;
    for (const auto& t : neighbours_left_cont)
    {
        total_left += t.second.size();
    }
    size_t total_right = 0;
    for (const auto& t : neighbours_right_cont)
    {
        total_left += t.second.size();
    }
    return std::make_tuple(total_left, neighbours_left_cont.size(),
                           total_right, neighbours_right_cont.size());
}

map< std::tuple< Contig_Entry_CBPtr, bool >, vector< Read_Chunk_CBPtr > >
Contig_Entry::out_chunks_dir(bool c_right, int unmappable_policy, size_t ignore_threshold) const
{
    map< std::tuple< Contig_Entry_CBPtr, bool >, vector< Read_Chunk_CBPtr > > res;
    bool c_left = not c_right;
    Size_Type endpoint = (c_left? 0 : get_len());
    for (auto rc_cbref : chunk_cont().iintersect(endpoint, endpoint))
    {
        Read_Chunk_CBPtr rc_cbptr = &rc_cbref;
        Read_Entry_BPtr re_cbptr = rc_cbptr->re_bptr();
        Read_Chunk_CBPtr rc_next_cbptr = re_cbptr->chunk_cont().get_sibling(rc_cbptr, false, c_right);
        if (not rc_next_cbptr)
        {
            continue;
        }
        if (unmappable_policy == 3)
        {
            // skip unmappable chunks
            while (rc_next_cbptr and rc_next_cbptr->is_unmappable())
            {
                rc_next_cbptr = re_cbptr->chunk_cont().get_sibling(rc_next_cbptr, false, c_right);
            }
            if (not rc_next_cbptr)
            {
                continue;
            }
        }
        ASSERT(rc_next_cbptr);
        Contig_Entry_CBPtr ce_next_cbptr = rc_next_cbptr->ce_bptr();
        bool same_orientation = (rc_cbptr->get_rc() == rc_next_cbptr->get_rc());
        if (rc_next_cbptr->is_unmappable())
        {
            ASSERT(unmappable_policy != 3);
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
        ASSERT(not c_right or rc_cbptr->get_c_end() == get_len());
        Read_Chunk_CBPtr rc_next_cbptr = re_cbptr->chunk_cont().get_sibling(rc_cbptr, false, c_right);
        ASSERT(rc_next_cbptr);
        Size_Type skip_len = 0;
        while (rc_next_cbptr->is_unmappable())
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

std::tuple< Contig_Entry_CBPtr, bool, vector< Read_Chunk_CBPtr > >
Contig_Entry::can_cat_dir(bool c_right) const
{
    // ignore unmappable out chunks iff this ce is not unmappable
    auto res = out_chunks_dir(c_right, not is_unmappable()? 0 : 1);
    if (res.size() != 1)
    {
        return std::make_tuple(Contig_Entry_CBPtr(nullptr), false, vector< Read_Chunk_CBPtr >());
    }
    ASSERT(not res.begin()->second.empty());
    Contig_Entry_CBPtr ce_next_cbptr;
    bool same_orientation;
    vector< Read_Chunk_CBPtr > rc_cbptr_cont = std::move(res.begin()->second);
    std::tie(ce_next_cbptr, same_orientation) = res.begin()->first;
    if (ce_next_cbptr->is_unmappable() != is_unmappable())
    {
        return std::make_tuple(Contig_Entry_CBPtr(nullptr), false, vector< Read_Chunk_CBPtr >());
    }
    auto tmp = ce_next_cbptr->out_chunks_dir(c_right != same_orientation, not is_unmappable()? 0 : 1);
    if (tmp.size() != 1)
    {
        return std::make_tuple(Contig_Entry_CBPtr(nullptr), false, vector< Read_Chunk_CBPtr >());
    }
    ASSERT(not tmp.begin()->second.empty());
    ASSERT(tmp.begin()->first == std::make_tuple(Contig_Entry_CBPtr(rc_cbptr_cont.front()->ce_bptr()), same_orientation));
    ASSERT(tmp.begin()->second.size() == rc_cbptr_cont.size());
    return std::make_tuple(ce_next_cbptr, same_orientation, std::move(rc_cbptr_cont));
}

void Contig_Entry::cat_c_right(Contig_Entry_BPtr ce_bptr, Contig_Entry_BPtr ce_next_bptr,
                               vector< MAC::Read_Chunk_CBPtr >& rc_cbptr_cont)
{
    ASSERT(ce_bptr->_seq_offset == 0 and ce_next_bptr->_seq_offset == 0);
    ASSERT(ce_next_bptr->is_unlinked());
    // first, shift all mutations and chunks from the second Contig_Entry
    ce_next_bptr->mut_cont().shift(int(ce_bptr->get_len()));
    ce_next_bptr->chunk_cont().shift(int(ce_bptr->get_len()));
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
map< std::tuple< Contig_Entry_CBPtr, bool >, std::tuple< unsigned int, Size_Type, Size_Type > >
Contig_Entry::neighbour_stats(bool c_right, int unmappable_policy, bool trim_results) const
{
    map< std::tuple< Contig_Entry_CBPtr, bool >, std::tuple< unsigned int, Size_Type, Size_Type > > res;
    //unsigned int support = 0;
    auto out_chunks_cont = out_chunks_dir(c_right, unmappable_policy);
    for (auto rc_cptr_it = _chunk_cptr_cont.begin(); rc_cptr_it != _chunk_cptr_cont.end(); ++rc_cptr_it)
    {
        Read_Chunk_CPtr rc_cptr = *rc_cptr_it;
        if ((not c_right and not rc_cptr->get_c_start() == 0)
                or (c_right and not rc_cptr->get_c_end() == get_len()))
        {
            continue;
        }
        Read_Chunk_CPtr rc_next_cptr = rc_cptr->get_re_ptr()->get_sibling(rc_cptr, c_right != rc_cptr->get_rc());
        Size_Type skipped_len = 0;
        if (rc_next_cptr != NULL and skip_unmappable and rc_next_cptr->is_unmappable())
        {
            skipped_len = rc_next_cptr->get_r_len();
            rc_next_cptr = rc_next_cptr->get_re_ptr()->get_sibling(rc_next_cptr, c_right != rc_cptr->get_rc());
        }
        if (rc_next_cptr == NULL or rc_next_cptr->is_unmappable())
        {
            continue;
        }
        //++support;
        std::tuple< const Contig_Entry*, bool > t = std::make_tuple(rc_next_cptr->get_ce_ptr(), rc_cptr->get_rc() != rc_next_cptr->get_rc());
        if (res->count(t) == 0)
        {
            (*res)[t] = std::make_tuple(1, skipped_len, skipped_len);
        }
        else
        {
            unsigned int tmp_support;
            Size_Type min_skipped_len;
            Size_Type max_skipped_len;
            std::tie(tmp_support, min_skipped_len, max_skipped_len) = (*res)[t];
            (*res)[t] = std::make_tuple(tmp_support + 1, min(min_skipped_len, skipped_len), max(max_skipped_len, skipped_len));
        }
    }
    // remove neighbours with single support
    if (skip_unmappable and trim_results)
    {
        for (auto it = res->begin(); it != res->end(); ++it)
        {
            if (get<0>(it->second) == 1)
            {
                / *
                const Contig_Entry* ce_next_cptr;
                bool flip;
                std::tie(ce_next_cptr, flip) = it->first;
                auto neighbours_next_sptr = ce_next_cptr->get_neighbours(dir == flip, true, false);
                ASSERT(neighbours_next_sptr->count(std::make_tuple(this, flip)) == 1);
                unsigned int neighbour_support = 0;
                for (auto it2 = neighbours_next_sptr->begin(); it2 != neighbours_next_sptr->end(); ++it2)
                {
                    neighbour_support += it2->second;
                }
                * /
                res->erase(it);
            }
        }
    }
    return res;
}
*/

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
        if (mut_cptr->get_start() >= last_mut_end + min_separation
                and mut_cptr->get_end() + min_separation <= get_len()
                and (next_mut_it == _mut_cont.end()
                     or mut_cptr->get_end() + min_separation <= next_mut_it->get_start()))
        {
            // well separated
            size_t n_chunks_supporting_mut[2] = { 0, 0 };
            for (auto& rc_cptr : _chunk_cptr_cont)
            {
                if (rc_cptr->get_c_start() > mut_cptr->get_start()
                        or rc_cptr->get_c_end() < mut_cptr->get_end())
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
        last_mut_end = std::max(last_mut_end, mut_cptr->get_end());
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
            os << _seq_ptr->substr(mut_cptr->get_start(), mut_cptr->get_len()) << "\t\t";
        }
        else
        {
            os << mut_cptr->get_seq() << '\t' << _seq_ptr->substr(mut_cptr->get_start(), mut_cptr->get_len()) << '\t';
        }
        os << _seq_ptr->substr(mut_cptr->get_start() - min_separation, min_separation) << '\t'
           << _seq_ptr->substr(mut_cptr->get_end(), min_separation) << '\n';
    }
}
*/

bool Contig_Entry::check() const
{
    // check there are chunks mapped to this contig
    ASSERT(_chunk_cont.size() > 0);
    // mutations:
    for (const auto& mut_cbref : mut_cont())
    {
        Mutation_CBPtr mut_cbptr = &mut_cbref;
        // check base coordinates
        ASSERT(mut_cbptr->get_start() <= seq().size() and mut_cbptr->get_end() <= seq().size());
        // check no empty mutations
        ASSERT(not mut_cbptr->is_empty());
        // check read_chunk_ptr_cont
        for (const auto& mca_cbref : mut_cbptr->chunk_ptr_cont())
        {
            Mutation_Chunk_Adapter_CBPtr mca_cbptr = &mca_cbref;
            // check mutation back pointers
            ASSERT(mca_cbptr->mut_cbptr() == mut_cbptr);
            // check read chunk part of current contig
            ASSERT(mca_cbptr->chunk_cbptr()->ce_bptr().raw() == this);
        }
    }
    // read chunks:
    auto ce_bptr = bptr_to();
    for (const auto& rc_cbref : chunk_cont())
    {
        Read_Chunk_CBPtr rc_cbptr = &rc_cbref;
        // check contig entry pointers
        ASSERT(rc_cbptr->ce_bptr() == ce_bptr);
        // check contig coordinates
        ASSERT(rc_cbptr->get_c_end() <= seq().size());
        // mutation pointers:
        for (const auto& mca_cbref : rc_cbptr->mut_ptr_cont())
        {
            Mutation_Chunk_Adapter_CBPtr mca_cbptr = &mca_cbref;
            // check back Read_Chunk pointers
            ASSERT(mca_cbptr->chunk_cbptr() == rc_cbptr);
            // check Mutation pointers point inside Mutation container
            ASSERT(mut_cont().find(mca_cbptr->mut_cbptr(), true));
        }
    }
    return true;
}

/*
bool Contig_Entry::check_colour(bool dir) const
{
    if (is_unmappable())
    {
        return true;
    }
    auto neighbours_sptr = get_neighbours(dir);
    if (neighbours_sptr->size() == 1)
    {
        ASSERT((get_colour() & (dir == false? 0x4 : 0x8)) == 0);
    }
    else
    {
        ASSERT((get_colour() & (dir == false? 0x4 : 0x8)) != 0);
    }
    return true;
}
*/

/*
ostream& operator << (ostream& os, const Contig_Entry& rhs)
{
    os << indent::tab << "(Contig_Entry &=" << (void*)&rhs
       << indent::inc << indent::nl << "seq=\"" << rhs.seq() << "\",len=" << rhs.get_len()
       << ",col=" << rhs.colour() << ",is_unmappable=" << (int)rhs.is_unmappable()
       << indent::nl << "mut_cont:"
       << indent::inc << '\n';
    print_seq(os, rhs._mut_cont, indent::nl, indent::tab, '\n');
    os << indent::dec << indent::tab << "chunk_cont:"
       << indent::inc << '\n';
    print_seq(os, rhs._chunk_cont, indent::nl, indent::tab, '\n');
    os << indent::dec << indent::dec << indent::tab << ")" << '\n';
    return os;
}
*/

boost::property_tree::ptree Contig_Entry::to_ptree() const
{
    return ptree().put("addr", (void*)this)
                  .put("seq", seq())
                  .put("seq_len", get_len())
                  .put("col", colour())
                  .put("is_unmappable", is_unmappable())
                  .put("mut_cont", cont_to_ptree(mut_cont()))
                  .put("chunk_cont", cont_to_ptree(chunk_cont()));
}

} // namespace MAC
