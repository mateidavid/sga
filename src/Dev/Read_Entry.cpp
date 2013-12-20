//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Read_Entry.hpp"
#include "Contig_Entry.hpp"

#include "print_seq.hpp"
#include "indent.hpp"

using namespace std;


namespace MAC
{
    const Read_Chunk* Read_Entry::get_chunk_with_pos(Size_Type r_pos) const
    {
        if (r_pos >= get_len())
            return NULL;
        Read_Chunk_Cont::const_iterator it = _chunk_cont.lower_bound(r_pos);
        if (it != _chunk_cont.end() and it->get_r_start() == r_pos)
        {
            return &(*it);
        }
        else
        {
            assert(it != _chunk_cont.begin());
            return &(*(--it));
        }
    }

    bool Read_Entry::is_terminal(bool check_start) const
    {
        Read_Chunk_CPtr rc_cptr = (check_start? &*_chunk_cont.begin() : &*_chunk_cont.rbegin());
        return ((check_start and not rc_cptr->get_rc() and rc_cptr->get_c_start() == 0)
                or (check_start and rc_cptr->get_rc() and rc_cptr->get_c_end() == rc_cptr->get_ce_ptr()->get_len())
                or (not check_start and not rc_cptr->get_rc() and rc_cptr->get_c_end() == rc_cptr->get_ce_ptr()->get_len())
                or (not check_start and rc_cptr->get_rc() and rc_cptr->get_c_start() == 0));
    }

    Seq_Type Read_Entry::get_seq() const
    {
        Seq_Type res;
        for (auto it = _chunk_cont.begin(); it != _chunk_cont.end(); ++it)
        {
            res += it->get_seq();
        }
        return res;
    }

    bool Read_Entry::check() const
    {
        // name is set
        assert(_name_ptr);
        // name not empty
        assert(_name_ptr->size() > 0);
        for (auto rc_it = _chunk_cont.begin(); rc_it != _chunk_cont.end(); ++rc_it)
        {
            Read_Chunk_CPtr rc_cptr = &*rc_it;
            // re_ptr
            assert(rc_cptr->get_re_ptr() == this);
            // no empty chunks
            assert(rc_cptr->get_r_len() > 0);
            // read start&end bounds
            if (rc_it == _chunk_cont.begin())
            {
                assert(rc_cptr->get_r_start() == 0);
            }
            auto rc_next_it = rc_it;
            rc_next_it++;
            if (rc_next_it == _chunk_cont.end())
            {
                assert(rc_cptr->get_r_end() == _len);
            }
            if (rc_next_it != _chunk_cont.end())
            {
                assert(rc_cptr->get_r_end() == rc_next_it->get_r_start());
            }
            // contigs coordinates
            assert(rc_cptr->get_c_start() <= rc_cptr->get_c_end());
            assert(rc_cptr->get_c_start() <= rc_cptr->get_ce_ptr()->get_seq_offset() + rc_cptr->get_ce_ptr()->get_len());
            assert(rc_cptr->get_c_end() <= rc_cptr->get_ce_ptr()->get_seq_offset() + rc_cptr->get_ce_ptr()->get_len());
            // mapped length
            Size_Type c_len = rc_cptr->get_c_end() - rc_cptr->get_c_start();
            Size_Type r_len = rc_cptr->get_r_end() - rc_cptr->get_r_start();
            long long delta = 0;
            for (size_t i = 0; i < rc_cptr->get_mut_ptr_cont().size(); ++i)
            {
                // no empty mutations
                assert(not rc_cptr->get_mut_ptr_cont()[i]->is_empty());
                // mutations must be in contig order
                assert(i == 0 or rc_cptr->get_mut_ptr_cont()[i - 1]->get_end() <= rc_cptr->get_mut_ptr_cont()[i]->get_start());
#ifndef ALLOW_CONSECUTIVE_MUTATIONS
                assert(i == 0 or rc_cptr->get_mut_ptr_cont()[i - 1]->get_end() < rc_cptr->get_mut_ptr_cont()[i]->get_start());
#endif
                delta += (long long)rc_cptr->get_mut_ptr_cont()[i]->get_seq_len() - (long long)rc_cptr->get_mut_ptr_cont()[i]->get_len();
            }
            assert((long long)c_len + delta == (long long)r_len);
#ifndef ALLOW_PART_MAPPED_INNER_CHUNKS
            // chunks must end on contig breaks except for first and last
            assert(rc_it == _chunk_cont.begin()
                   or (not rc_cptr->get_rc()?
                       rc_cptr->get_c_start() == 0
                       : rc_cptr->get_c_end() == rc_cptr->get_ce_ptr()->get_seq_offset() + rc_cptr->get_ce_ptr()->get_len()));
            assert(rc_next_it == _chunk_cont.end()
                   or (not rc_cptr->get_rc()?
                       rc_cptr->get_c_end() == rc_cptr->get_ce_ptr()->get_seq_offset() + rc_cptr->get_ce_ptr()->get_len()
                       : rc_cptr->get_c_start() == 0));
#endif
        }
        return true;
    }

    ostream& operator << (ostream& os, const Read_Entry& rhs)
    {
        os << indent::tab << "(Read_Entry &=" << (void*)&rhs
           << indent::inc << indent::nl << "name=\"" << rhs.get_name() << "\",len=" << rhs.get_len()
           << indent::nl << "chunk_cont:\n"
           << indent::inc;
        print_seq(os, rhs._chunk_cont, "");
        os << indent::dec << indent::dec << indent::tab << ")\n";
        return os;
    }
}
