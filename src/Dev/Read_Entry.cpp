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

bool Read_Entry::is_terminal(bool check_start) const
{
    Read_Chunk_BPtr rc_bptr = (check_start? &*_chunk_cont.begin() : &*_chunk_cont.rbegin());
    return ((check_start and not rc_bptr->get_rc() and rc_bptr->get_c_start() == 0)
            or (check_start and rc_bptr->get_rc() and rc_bptr->get_c_end() == rc_bptr->ce_bptr()->get_len())
            or (not check_start and not rc_bptr->get_rc() and rc_bptr->get_c_end() == rc_bptr->ce_bptr()->get_len())
            or (not check_start and rc_bptr->get_rc() and rc_bptr->get_c_start() == 0));
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
    ASSERT(_name_ptr);
    // name not empty
    ASSERT(_name_ptr->size() > 0);
    for (auto rc_cit = _chunk_cont.begin(); rc_cit != _chunk_cont.end(); ++rc_cit)
    {
        Read_Chunk_CBPtr rc_cbptr = &*rc_cit;
        // re_ptr
        ASSERT(rc_cbptr->re_bptr() == this->bptr_to());
        // read start&end bounds
        if (rc_cit == _chunk_cont.begin())
        {
            ASSERT(rc_cbptr->get_r_start() == 0);
        }
        auto rc_next_cit = rc_cit;
        rc_next_cit++;
        if (rc_next_cit == _chunk_cont.end())
        {
            ASSERT(rc_cbptr->get_r_end() == _len);
        }
         if (rc_next_cit != _chunk_cont.end())
         {
             ASSERT(rc_cbptr->get_r_end() == rc_next_cit->get_r_start());
         }
         ASSERT(rc_cbptr->check());
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
