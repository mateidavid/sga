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
    /*
    void Contig_Entry::remove_chunk(Read_Chunk_CPtr rc_cptr)
    {
        Read_Chunk_CPtr_Cont::iterator it, it_end;
        for (auto it = _chunk_cptr_cont.begin(); it != _chunk_cptr_cont.end(); ++it)
        {
            if (*it == rc_cptr)
            {
                _chunk_cptr_cont.erase(it);
                return;
            }
        }
        ASSERT(false);
    }

    void Contig_Entry::remove_chunks(const set< Read_Chunk_CPtr >& rc_cptr_set)
    {
        size_t i = 0;
        while (i < _chunk_cptr_cont.size())
        {
            if (rc_cptr_set.count(_chunk_cptr_cont[i]) > 0)
                _chunk_cptr_cont.erase(_chunk_cptr_cont.begin() + i);
            else
                ++i;
        }
    }

    const Mutation* Contig_Entry::cut_mutation(const Mutation* mut_cptr, Size_Type c_offset, Size_Type r_offset)
    {
        // construct modifier that cuts existing mutation, saving remaining part
        Mutation m_new;
        auto mut_modifier = [&] (Mutation& m) { m_new = m.cut(c_offset, r_offset); };

        // apply it
        modify_element<Mutation_Cont>(_mut_cont, mut_cptr, mut_modifier);

        // insert remaining part
        Mutation_Cont::iterator it_new;
        bool success;
        tie(it_new, success) = _mut_cont.insert(m_new);
        ASSERT(success);

        return &(*it_new);
    }

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

    vector< const Mutation* > Contig_Entry::get_mutations_spanning_pos(Size_Type c_pos) const
    {
        vector< const Mutation* > res;
        for (auto it = _mut_cont.begin(); it != _mut_cont.end() and it->get_start() < c_pos; ++it)
        {
            if (c_pos < it->get_end())
                res.push_back(&(*it));
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

    void Contig_Entry::drop_unused_mutations()
    {
        set< Mutation_CPtr > used_mut;
        for (auto rc_cptr_it = _chunk_cptr_cont.begin(); rc_cptr_it != _chunk_cptr_cont.end(); ++rc_cptr_it)
        {
            Read_Chunk_CPtr rc_cptr = *rc_cptr_it;
            for (auto mut_cptr_it = rc_cptr->get_mut_ptr_cont().begin(); mut_cptr_it != rc_cptr->get_mut_ptr_cont().end(); ++mut_cptr_it)
            {
                Mutation_CPtr mut_cptr = *mut_cptr_it;
                used_mut.insert(mut_cptr);
            }
        }
        for (auto mut_it = _mut_cont.begin(); mut_it != _mut_cont.end(); ++mut_it)
        {
            if (used_mut.count(&*mut_it) == 0)
            {
                _mut_cont.erase(mut_it);
            }
        }
    }

    void Contig_Entry::reverse(const Read_Chunk::ext_mod_type& rc_reverse_mod)
    {
        // only reverse full contigs
        ASSERT(_seq_offset == 0);

        // reverse the string
        std::shared_ptr< Seq_Type > old_seq_ptr = _seq_ptr;
        _seq_ptr.reset(new string(reverseComplement(*old_seq_ptr)));

        // save mutation pointers
        vector< Mutation_CPtr > mut_cptr_cont;
        for (auto mut_it = _mut_cont.begin(); mut_it != _mut_cont.end(); ++mut_it)
        {
            mut_cptr_cont.push_back(&*mut_it);
        }

        // using saved pointers, change mutations in place (incudes re-keying)
        for (auto mut_cptr_it = mut_cptr_cont.begin(); mut_cptr_it != mut_cptr_cont.end(); ++mut_cptr_it)
        {
            Mutation_CPtr mut_cptr = *mut_cptr_it;
            modify_element<Mutation_Cont>(_mut_cont, mut_cptr, [&] (Mutation& mut) { mut.reverse(get_len()); });
        }

        // apply external modifier to reverse read chunks in-place
        for (auto rc_cptr_it = _chunk_cptr_cont.begin(); rc_cptr_it != _chunk_cptr_cont.end(); ++rc_cptr_it)
        {
            Read_Chunk_CPtr rc_cptr = *rc_cptr_it;
            rc_reverse_mod(rc_cptr);
        }
    }

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

    Read_Chunk_CPtr Contig_Entry::get_next_chunk(bool dir, Read_Chunk_CPtr rc_cptr) const
    {
        return rc_cptr->get_re_ptr()->get_sibling(rc_cptr, dir != rc_cptr->get_rc());
    }

    shared_ptr< vector< Read_Chunk_CPtr > > Contig_Entry::get_chunks_out(bool dir, bool skip_next_unmappable) const
    {
        shared_ptr< vector< Read_Chunk_CPtr > > res(new vector< Read_Chunk_CPtr >());
        for (auto rc_cptr_it = _chunk_cptr_cont.begin(); rc_cptr_it != _chunk_cptr_cont.end(); ++rc_cptr_it)
        {
            Read_Chunk_CPtr rc_cptr = *rc_cptr_it;
            if ((not dir and not rc_cptr->get_c_start() == 0)
                or (dir and not rc_cptr->get_c_end() == get_len()))
            {
                continue;
            }
            Read_Chunk_CPtr rc_next_cptr = rc_cptr->get_re_ptr()->get_sibling(rc_cptr, dir != rc_cptr->get_rc());
            if (rc_next_cptr == NULL or (skip_next_unmappable and rc_next_cptr->is_unmappable()))
            {
                continue;
            }
            res->push_back(rc_cptr);
        }
        return res;
    }

    shared_ptr< vector< Read_Chunk_CPtr > > Contig_Entry::is_mergeable_one_way(bool dir) const
    {
        auto chunks_out_cont_sptr = get_chunks_out(dir, false);
        if (chunks_out_cont_sptr->size() == 0)
            return NULL;
        // check all chunks leaving go to the same (different) contig, in the same orientation
        const Contig_Entry* candidate_ce_cptr = NULL;
        bool candidate_orientation;
        for (auto rc_cptr_it = chunks_out_cont_sptr->begin(); rc_cptr_it != chunks_out_cont_sptr->end(); ++rc_cptr_it)
        {
            Read_Chunk_CPtr rc_cptr = *rc_cptr_it;
            Read_Chunk_CPtr rc_next_cptr = get_next_chunk(dir, rc_cptr);
            ASSERT(rc_next_cptr != NULL);
            //ASSERT(not rc_next_cptr->is_unmappable());
            const Contig_Entry* tmp_ce_cptr = rc_next_cptr->get_ce_ptr();
            bool tmp_orientation = (rc_cptr->get_rc() == rc_next_cptr->get_rc());
            if (tmp_ce_cptr == this) // multiple chunks of the same read mapped across this boundary: unmergeable
                return NULL;
            if (candidate_ce_cptr == NULL)
            {
                candidate_ce_cptr = tmp_ce_cptr;
                candidate_orientation = tmp_orientation;
            }
            else
            {
                if (tmp_ce_cptr != candidate_ce_cptr or tmp_orientation != candidate_orientation)
                    return NULL;
            }
        }
        ASSERT(candidate_ce_cptr != NULL);
        return chunks_out_cont_sptr;
    }

    shared_ptr< vector< Read_Chunk_CPtr > > Contig_Entry::is_mergeable(bool dir) const
    {
        auto chunks_out_cont_sptr = is_mergeable_one_way(dir);
        if (not chunks_out_cont_sptr)
            return NULL;
        Read_Chunk_CPtr rc_cptr = *chunks_out_cont_sptr->begin();
        Read_Chunk_CPtr rc_next_cptr = get_next_chunk(dir, rc_cptr);
        const Contig_Entry* ce_next_cptr = rc_next_cptr->get_ce_ptr();
        bool same_orientation = (rc_cptr->get_rc() == rc_next_cptr->get_rc());
        auto chunks_in_cont_sptr = ce_next_cptr->is_mergeable_one_way(dir != same_orientation);
        if (not chunks_in_cont_sptr)
            return NULL;
        ASSERT(chunks_out_cont_sptr->size() == chunks_in_cont_sptr->size());
        return chunks_out_cont_sptr;
    }

    void Contig_Entry::merge_forward(const Contig_Entry* ce_next_cptr, const Read_Chunk::ext_mod_with_map_type& rc_rebase_mod)
    {
        ASSERT(_seq_offset == 0 and ce_next_cptr->_seq_offset == 0);
        Mutation_Trans_Cont mut_map;
        // rebase mutations: copy them into this contig, and in the translation table
        for (auto mut_it = ce_next_cptr->_mut_cont.begin(); mut_it != ce_next_cptr->_mut_cont.end(); ++mut_it)
        {
            Mutation_Trans mut_trans;
            mut_trans.old_mut_cptr = &*mut_it;
            Mutation m(*mut_it);
            m.add_base_prefix(get_len());
            Mutation_Cont::iterator new_mut_it;
            bool success;
            std::tie(new_mut_it, success) = _mut_cont.insert(m);
            ASSERT(success);
            mut_trans.new_mut_cptr = &*new_mut_it;
            Mutation_Trans_Cont::iterator it;
            std::tie(it, success) = mut_map.insert(mut_trans);
            ASSERT(success);
        }
        // rebase chunks using external modifier
        for (auto rc_cptr_it = ce_next_cptr->_chunk_cptr_cont.begin(); rc_cptr_it != ce_next_cptr->_chunk_cptr_cont.end(); ++rc_cptr_it)
        {
            Read_Chunk_CPtr rc_cptr = *rc_cptr_it;
            rc_rebase_mod(rc_cptr, mut_map);
            _chunk_cptr_cont.push_back(rc_cptr);
        }
        // lastly, merge the next contig's sequence into this one
        *_seq_ptr += *ce_next_cptr->_seq_ptr;
    }

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
        // check base sequence exists
        ASSERT(_seq_ptr);
        // check there are chunks mapped to this contig
        ASSERT(_chunk_cont.size() > 0);
        // mutations:
        for (const auto& mut_cbref : _mut_cont)
        {
            const Mutation& mut = mut_cbref.raw();
            // check base coordinates
            ASSERT(mut.get_start() <= _seq_ptr->size() and mut.get_end() <= _seq_ptr->size());
            // check no empty mutations
            ASSERT(not mut.is_empty());
        }
        // read chunks:
        auto ce_bptr = bptr_to();
        for (const auto& rc_cbref : _chunk_cont)
        {
            const Read_Chunk& rc = rc_cbref.raw();
            // check contig entry pointers
            ASSERT(rc.ce_bptr() == ce_bptr);
            // check contig coordinates
            ASSERT(rc.get_c_end() <= _seq_ptr->size());
            // mutation pointers:
            for (const auto& mca_cbref : rc.mut_ptr_cont())
            {
                Mutation_CBPtr mut_cbptr = mca_cbref.raw().mut_cbptr();
                // check they point inside mutation container
                ASSERT(_mut_cont.find(mut_cbptr) != _mut_cont.end());
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
           << indent::inc << indent::nl << "seq=\"" << rhs.get_seq() << "\",len=" << rhs.get_len()
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
}
