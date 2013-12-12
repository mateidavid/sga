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

using namespace std;


namespace MAC
{
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
        assert(false);
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
        assert(success);

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

        assert(mut_left_cptr == NULL or (mut_left_cptr->get_start() == c_pos and mut_left_cptr->is_ins()));

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
                assert(success);
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
                assert(rc_next_cptr != NULL);
                set_0.insert(rc_next_cptr->get_ce_ptr());
            }
            if (rc_cptr->get_c_end() == get_len()
                and (not rc_cptr->get_rc()? rc_cptr->get_r_end() < rc_cptr->get_re_ptr()->get_len() : rc_cptr->get_r_start() > 0))
            {
                ++cnt_1;
                Read_Chunk_CPtr rc_next_cptr = rc_cptr->get_re_ptr()->get_sibling(rc_cptr, not rc_cptr->get_rc()? true : false);
                assert(rc_next_cptr != NULL);
                set_1.insert(rc_next_cptr->get_ce_ptr());
            }
        }
        return std::make_tuple(cnt_0, set_0.size(), cnt_1, set_1.size());
    }

    bool Contig_Entry::check() const
    {
        // check base sequence exists
        assert(_seq_ptr);
        // mutations:
        for (auto mut_it = _mut_cont.begin(); mut_it != _mut_cont.end(); ++mut_it)
        {
            // check base coordinates
            assert(mut_it->get_start() <= _seq_ptr->size() and mut_it->get_end() <= _seq_ptr->size());
        }
        // read chunks:
        for (auto rc_cptr_it = _chunk_cptr_cont.begin(); rc_cptr_it != _chunk_cptr_cont.end(); ++rc_cptr_it)
        {
            // check contig pointers
            assert((*rc_cptr_it)->get_ce_ptr() == this);
            // check contig coordinates
            assert((*rc_cptr_it)->get_c_end() <= _seq_ptr->size());
            // mutation pointers:
            for (auto mut_cptr_it = (*rc_cptr_it)->get_mut_ptr_cont().begin(); mut_cptr_it != (*rc_cptr_it)->get_mut_ptr_cont().end(); ++mut_cptr_it)
            {
                // check they point inside mutation container
                auto mut_it = _mut_cont.begin();
                while (mut_it != _mut_cont.end() and *mut_cptr_it != &(*mut_it))
                {
                    ++mut_it;
                }
                assert(mut_it != _mut_cont.end());
            }
        }
        return true;
    }

    std::ostream& operator << (std::ostream& os, const Contig_Entry& rhs)
    {
        os << indent::tab << "(Contig_Entry &=" << (void*)&rhs
           << indent::inc << indent::nl << "seq=\"" << rhs.get_seq() << "\""
           << indent::nl << "mut_cont:"
           << indent::inc << '\n';
        print_seq(os, rhs._mut_cont, indent::nl, indent::tab, '\n');
        os << indent::dec << indent::tab << "chunk_cptr_cont:"
           << indent::inc << '\n';
        print_seq(os, rhs._chunk_cptr_cont, indent::nl, indent::tab, '\n');
        os << indent::dec << indent::dec << indent::tab << ")" << '\n';
        return os;
    }
}
