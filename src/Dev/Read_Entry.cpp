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

    void Read_Entry::check() const
    {
        // name is set
        assert(_name_ptr);
        // name not empty
        assert(_name_ptr->size() > 0);
        for (auto rc_it = _chunk_cont.begin(); rc_it != _chunk_cont.end(); ++rc_it)
        {
            // re_ptr
            assert(rc_it->get_re_ptr() == this);
            // no empty chunks
            assert(rc_it->get_r_len() > 0);
            // read start&end bounds
            if (rc_it == _chunk_cont.begin())
            {
                assert(rc_it->get_r_start() == 0);
            }
            auto rc_next_it = rc_it;
            rc_next_it++;
            if (rc_next_it == _chunk_cont.end())
            {
                assert(rc_it->get_r_end() == _len);
            }
            if (rc_next_it != _chunk_cont.end())
            {
                assert(rc_it->get_r_end() == rc_next_it->get_r_start());
            }
        }
    }

    ostream& operator << (ostream& os, const Read_Entry& rhs)
    {
        os << "(Read_Entry &=" << (void*)&rhs
           << indent::inc << indent::nl
           << "name=" << rhs.get_name()
           << indent::nl << "chunk_cont="
           << indent::inc;
        print_seq(os, rhs._chunk_cont, indent::nl, indent::nl);
        os << indent::dec << indent::dec << indent::nl << ")";
        return os;
    }
}
