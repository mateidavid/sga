//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Read_Entry.hpp"
#include "Contig_Entry.hpp"

using namespace std;


namespace MAC
{

bool Read_Entry::is_terminal(bool check_start) const
{
    // get the first/last chunk
    Read_Chunk_CBPtr rc_cbptr = (check_start? &*_chunk_cont.begin() : &*_chunk_cont.rbegin());
    return ((check_start and not rc_cbptr->get_rc() and rc_cbptr->get_c_start() == 0)
            or (check_start and rc_cbptr->get_rc() and rc_cbptr->get_c_end() == rc_cbptr->ce_bptr()->get_len())
            or (not check_start and not rc_cbptr->get_rc() and rc_cbptr->get_c_end() == rc_cbptr->ce_bptr()->get_len())
            or (not check_start and rc_cbptr->get_rc() and rc_cbptr->get_c_start() == 0));
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
    // name not empty
    ASSERT(not _name.empty());
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

boost::property_tree::ptree Read_Entry::to_ptree() const
{
    return ptree().put("name", get_name())
                  .put("len", get_len())
                  .put("chunk_cont", cont_to_ptree(chunk_cont()));
}

} // namespace MAC
