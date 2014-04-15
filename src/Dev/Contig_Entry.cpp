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

            rc_bptr->mut_ptr_cont().insert_after(rc_bptr->mut_ptr_cont().iterator_to(&*mca_it), new_mca_bptr);

            ++mca_it;
            ++new_mca_it;
        }
        ASSERT(new_mca_it == new_mut_bptr->chunk_ptr_cont().end());
    }

    ASSERT(check());
}

/*
vector< Read_Chunk_CPtr > Contig_Entry::get_chunks_with_mutation(const Mutation* mut_cptr) const
{
    vector< Read_Chunk_CPtr > res;
    for (auto it = _chunk_cptr_cont.begin(); it != _chunk_cptr_cont.end(); ++it)
    {
        if ((*it)->have_mutation(mut_cptr))
            res.push_back(*it);
    }
    return res;
}

map< const Mutation*, const Mutation* > Contig_Entry::acquire_second_half_mutations(
    const Contig_Entry* ce_cptr, Size_Type c_pos, const Mutation* mut_left_cptr)
{
    map< const Mutation*, const Mutation* > res;

    ASSERT(mut_left_cptr == NULL or (mut_left_cptr->get_start() == c_pos and mut_left_cptr->is_ins()));

    Mutation_Cont::iterator mut_old_it = ce_cptr->_mut_cont.lower_bound(c_pos);
    while (mut_old_it != ce_cptr->_mut_cont.end())
    {
        if (&(*mut_old_it) != mut_left_cptr)
        {
            Mutation_Cont::iterator mut_new_it;
            bool success;
            if (mut_old_it->have_seq())
                tie(mut_new_it, success) = _mut_cont.insert(
                                               Mutation(mut_old_it->get_start() - c_pos, mut_old_it->get_len(), mut_old_it->get_seq()));
            else
                tie(mut_new_it, success) = _mut_cont.insert(
                                               Mutation(mut_old_it->get_start() - c_pos, mut_old_it->get_len(), mut_old_it->get_seq_len()));
            ASSERT(success);
            res[&(*mut_old_it)] = &(*mut_new_it);
        }
        ++mut_old_it;
    }
    return res;
}

void Contig_Entry::drop_mutations(const map< const Mutation*, const Mutation* >& mut_cptr_map)
{
    for (auto it = mut_cptr_map.begin(); it != mut_cptr_map.end(); ++it)
    {
        _mut_cont.erase(_mut_cont.iterator_to(*(it->first)));
    }
}
*/

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

/*
tuple< size_t, size_t, size_t, size_t > Contig_Entry::get_out_degrees() const
{
    size_t cnt_0 = 0;
    size_t cnt_1 = 0;
    set< const Contig_Entry* > set_0;
    set< const Contig_Entry* > set_1;
    for (auto rc_cptr_it = _chunk_cptr_cont.begin(); rc_cptr_it != _chunk_cptr_cont.end(); ++rc_cptr_it)
    {
        Read_Chunk_CPtr rc_cptr = *rc_cptr_it;
        if (rc_cptr->get_c_start() == 0
                and (not rc_cptr->get_rc()? rc_cptr->get_r_start() > 0 : rc_cptr->get_r_end() < rc_cptr->get_re_ptr()->get_len()))
        {
            ++cnt_0;
            Read_Chunk_CPtr rc_next_cptr = rc_cptr->get_re_ptr()->get_sibling(rc_cptr, not rc_cptr->get_rc()? false : true);
            ASSERT(rc_next_cptr != NULL);
            set_0.insert(rc_next_cptr->get_ce_ptr());
        }
        if (rc_cptr->get_c_end() == get_len()
                and (not rc_cptr->get_rc()? rc_cptr->get_r_end() < rc_cptr->get_re_ptr()->get_len() : rc_cptr->get_r_start() > 0))
        {
            ++cnt_1;
            Read_Chunk_CPtr rc_next_cptr = rc_cptr->get_re_ptr()->get_sibling(rc_cptr, not rc_cptr->get_rc()? true : false);
            ASSERT(rc_next_cptr != NULL);
            set_1.insert(rc_next_cptr->get_ce_ptr());
        }
    }
    return std::make_tuple(cnt_0, set_0.size(), cnt_1, set_1.size());
}

shared_ptr< vector< Read_Chunk_CPtr > > Contig_Entry::get_chunks_spanning_pos(Size_Type start, Size_Type end) const
{
    shared_ptr< vector< Read_Chunk_CPtr > > res(new vector< Read_Chunk_CPtr >());
    for (auto rc_cptr_it = _chunk_cptr_cont.begin(); rc_cptr_it != _chunk_cptr_cont.end(); ++rc_cptr_it)
    {
        Read_Chunk_CPtr rc_cptr = *rc_cptr_it;
        if (rc_cptr->get_c_start() <= start and rc_cptr->get_c_end() >= end)
        {
            res->push_back(rc_cptr);
        }
    }
    return res;
}
*/

map< std::tuple< Contig_Entry_CBPtr, bool >, vector< Read_Chunk_CBPtr > >
Contig_Entry::out_chunks_dir(bool c_right, bool skip_next_unmappable) const
{
    map< std::tuple< Contig_Entry_CBPtr, bool >, vector< Read_Chunk_CBPtr > > res;
    bool c_left = not c_right;
    Size_Type endpoint = (c_left? 0 : get_len());
    for (auto rc_cbref : chunk_cont().interval_intersect(endpoint, endpoint))
    {
        Read_Chunk_CBPtr rc_cbptr = &rc_cbref;
        Read_Chunk_CBPtr rc_next_cbptr = rc_cbptr->re_bptr()->chunk_cont().get_sibling(rc_cbptr, false, c_right);
        if (not rc_next_cbptr or (skip_next_unmappable and rc_next_cbptr->is_unmappable()))
        {
            continue;
        }
        Contig_Entry_CBPtr ce_next_cbptr = rc_next_cbptr->ce_bptr();
        bool same_orientation = (rc_cbptr->get_rc() == rc_next_cbptr->get_rc());
        res[std::make_tuple(ce_next_cbptr, same_orientation)].push_back(rc_cbptr);
    }
    return res;
}

std::tuple< Contig_Entry_CBPtr, bool, vector< Read_Chunk_CBPtr > >
Contig_Entry::can_cat_dir(bool c_right) const
{
    auto res = out_chunks_dir(c_right, false);
    if (res.size() != 1)
    {
        return std::make_tuple(Contig_Entry_CBPtr(nullptr), false, vector< Read_Chunk_CBPtr >());
    }
    ASSERT(not res.begin()->second.empty());
    Contig_Entry_CBPtr ce_next_cbptr;
    bool same_orientation;
    vector< Read_Chunk_CBPtr > rc_cbptr_cont = std::move(res.begin()->second);
    std::tie(ce_next_cbptr, same_orientation) = res.begin()->first;
    auto tmp = ce_next_cbptr->out_chunks_dir(c_right != same_orientation);
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
shared_ptr< map< std::tuple< const Contig_Entry*, bool >, std::tuple< unsigned int, Size_Type, Size_Type > > >
Contig_Entry::get_neighbours(bool dir, bool skip_unmappable, bool trim_results) const
{
    shared_ptr< map< std::tuple< const Contig_Entry*, bool >, std::tuple< unsigned int, Size_Type, Size_Type > > >
    res(new map< std::tuple< const Contig_Entry*, bool >, std::tuple< unsigned int, Size_Type, Size_Type > >());
    //unsigned int support = 0;
    for (auto rc_cptr_it = _chunk_cptr_cont.begin(); rc_cptr_it != _chunk_cptr_cont.end(); ++rc_cptr_it)
    {
        Read_Chunk_CPtr rc_cptr = *rc_cptr_it;
        if ((not dir and not rc_cptr->get_c_start() == 0)
                or (dir and not rc_cptr->get_c_end() == get_len()))
        {
            continue;
        }
        Read_Chunk_CPtr rc_next_cptr = rc_cptr->get_re_ptr()->get_sibling(rc_cptr, dir != rc_cptr->get_rc());
        Size_Type skipped_len = 0;
        if (rc_next_cptr != NULL and skip_unmappable and rc_next_cptr->is_unmappable())
        {
            skipped_len = rc_next_cptr->get_r_len();
            rc_next_cptr = rc_next_cptr->get_re_ptr()->get_sibling(rc_next_cptr, dir != rc_cptr->get_rc());
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
    for (const auto& mut_cbref : _mut_cont)
    {
        const Mutation& mut = mut_cbref.raw();
        // check base coordinates
        ASSERT(mut.get_start() <= _seq.size() and mut.get_end() <= _seq.size());
        // check no empty mutations
        ASSERT(not mut.is_empty());
    }
    // read chunks:
    auto ce_bptr = bptr_to();
    for (const auto& rc_cbref : _chunk_cont)
    {
        Read_Chunk_CBPtr rc_cbptr = &rc_cbref;
        // check contig entry pointers
        ASSERT(rc_cbptr->ce_bptr() == ce_bptr);
        // check contig coordinates
        ASSERT(rc_cbptr->get_c_end() <= _seq.size());
        // mutation pointers:
        for (const auto& mca_cbref : rc_cbptr->mut_ptr_cont())
        {
            Mutation_Chunk_Adapter_CBPtr mca_cbptr = &mca_cbref;
            // check back Read_Chunk pointers
            ASSERT(mca_cbptr->chunk_cbptr() == rc_cbptr);
            // check Mutation pointers point inside Mutation container
            ASSERT(_mut_cont.find(mca_cbptr->mut_cbptr(), true));
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

ostream& operator << (ostream& os, const Contig_Entry& rhs)
{
    os << indent::tab << "(Contig_Entry &=" << (void*)&rhs
       << indent::inc << indent::nl << "seq=\"" << rhs.seq() << "\",len=" << rhs.get_len()
       << ",col=" << rhs.colour() << ",is_unmappable=" << (int)rhs.is_unmappable()
       << indent::nl << "mut_cont:"
       << indent::inc << '\n';
    print_seq(os, rhs._mut_cont, indent::nl, indent::tab, '\n');
    os << indent::dec << indent::tab << "chunk_cptr_cont:"
       << indent::inc << '\n';
    print_seq(os, rhs._chunk_cont, indent::nl, indent::tab, '\n');
    os << indent::dec << indent::dec << indent::tab << ")" << '\n';
    return os;
}

} // namespace MAC
