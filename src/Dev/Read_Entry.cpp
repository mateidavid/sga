//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Read_Entry.hpp"
#include "Contig_Entry.hpp"


namespace MAC
{

bool Read_Entry::is_terminal(bool check_start) const
{
    // get the first/last chunk
    Read_Chunk_CBPtr rc_cbptr = (check_start? &*_chunk_cont.begin() : &*_chunk_cont.rbegin());
    return ((check_start and not rc_cbptr->get_rc() and rc_cbptr->get_c_start() == 0)
            or (check_start and rc_cbptr->get_rc() and rc_cbptr->get_c_end() == rc_cbptr->ce_bptr()->len())
            or (not check_start and not rc_cbptr->get_rc() and rc_cbptr->get_c_end() == rc_cbptr->ce_bptr()->len())
            or (not check_start and rc_cbptr->get_rc() and rc_cbptr->get_c_start() == 0));
}

Seq_Type Read_Entry::get_seq() const
{
    Seq_Type res;
    for (auto rc_cbptr : chunk_cont() | referenced)
    {
        res += rc_cbptr->get_seq();
    }
    return res;
}

void Read_Entry::check() const
{
#ifndef DISABLE_ASSERTS
    // name not empty
    ASSERT(not name().empty());
    // check integrity of Read_Chunk container
    chunk_cont().check();
    // check Read_Chunk contiguity
    for (auto rc_cit = chunk_cont().begin(); rc_cit != chunk_cont().end(); ++rc_cit)
    {
        Read_Chunk_CBPtr rc_cbptr = &*rc_cit;
        // re_bptr points here
        ASSERT(rc_cbptr->re_bptr() and rc_cbptr->re_bptr().raw() == this);
        // read start&end bounds
        if (rc_cit == chunk_cont().begin())
        {
            ASSERT(rc_cbptr->get_r_start() == start());
        }
        auto rc_next_cit = rc_cit;
        rc_next_cit++;
        if (rc_next_cit == chunk_cont().end())
        {
            ASSERT(rc_cbptr->get_r_end() == end());
        }
        else
        {
            ASSERT(rc_cbptr->get_r_end() == rc_next_cit->get_r_start());
        }
    }
    // check individual chunks
    for (auto rc_cbptr : chunk_cont() | referenced)
    {
        rc_cbptr->check();
    }
#endif
}

boost::property_tree::ptree Read_Entry::to_ptree() const
{
    return ptree().put("name", name())
                  .put("len", len())
                  .put("chunk_cont", cont_to_ptree(chunk_cont()));
}

void Read_Entry::save_strings(ostream& os, size_t& n_strings, size_t& n_bytes) const
{
    os.write(_name.c_str(), _name.size() + 1);
    ++n_strings;
    n_bytes += _name.size() + 1;
}

void Read_Entry::init_strings()
{
    new (&_name) string();
}

void Read_Entry::load_strings(istream& is, size_t& n_strings, size_t& n_bytes)
{
    getline(is, _name, '\0');
    ++n_strings;
    n_bytes += _name.size() + 1;
}

} // namespace MAC
