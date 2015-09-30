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

pair< Read_Entry_BPtr, Read_Entry_BPtr >
Read_Entry::split(Read_Chunk_BPtr rc_bptr, Size_Type c_overlap_start, Size_Type c_overlap_end)
{
    LOG("Read_Entry", debug) << ptree("begin")
        .put("rc_bptr", rc_bptr.to_int())
        .put("c_overlap_start", c_overlap_start)
        .put("c_overlap_end", c_overlap_end);
    ASSERT(c_overlap_start <= c_overlap_end);
    auto re_bptr = rc_bptr->re_bptr();
    auto ce_bptr = rc_bptr->ce_bptr();
    ASSERT(re_bptr->is_unlinked());
    ASSERT(rc_bptr->get_c_start() <= c_overlap_start);
    ASSERT(c_overlap_end <= rc_bptr->get_c_end());
    // create new Read_Entry for tail, and copy following chunks
    auto tail_re_bptr = Read_Entry_Fact::new_elem();
    auto rc_it = re_bptr->chunk_cont().iterator_to(*rc_bptr);
    ASSERT(rc_it != re_bptr->chunk_cont().end());
    tail_re_bptr->chunk_cont().splice(re_bptr->chunk_cont(), next(rc_it), tail_re_bptr);
    // erase chunk from all containers
    re_bptr->chunk_cont().erase(rc_bptr);
    ce_bptr->chunk_cont().erase(rc_bptr);
    // cut chunk
    auto rc_p = Read_Chunk::split(rc_bptr, c_overlap_end, nullptr, true);
    // extend rhs to include overlap
    Size_Type match_len = c_overlap_end - c_overlap_start;
    rc_p.second->c_start() -= match_len;
    rc_p.second->c_len() += match_len;
    if (not rc_p.second->get_rc())
    {
        rc_p.second->r_start() -= match_len;
    }
    rc_p.second->r_len() += match_len;
    ASSERT(rc_p.first->get_c_end() == c_overlap_end);
    ASSERT(rc_p.second->get_c_start() == c_overlap_start);
    // place new chunks in their respective RE containers
    Size_Type re_end_pos;
    Size_Type tail_re_start_pos;
    if (not rc_p.first->get_rc())
    {
        re_bptr->chunk_cont().insert(rc_p.first);
        rc_p.second->re_bptr() = tail_re_bptr;
        tail_re_bptr->chunk_cont().insert(rc_p.second);
        re_end_pos = rc_p.first->get_r_end();
        tail_re_start_pos = rc_p.second->get_r_start();
    }
    else
    {
        re_bptr->chunk_cont().insert(rc_p.second);
        rc_p.first->re_bptr() = tail_re_bptr;
        tail_re_bptr->chunk_cont().insert(rc_p.first);
        re_end_pos = rc_p.second->get_r_end();
        tail_re_start_pos = rc_p.first->get_r_start();
    }
    LOG("Read_Entry", debug) << ptree("pos")
        .put("re_end_pos", re_end_pos)
        .put("tail_re_start_pos", tail_re_start_pos);
    // place chunks in the CE container
    ce_bptr->chunk_cont().insert(rc_p.first);
    ce_bptr->chunk_cont().insert(rc_p.second);
    // adjust RE names
    LOG("Read_Entry", info) << ptree("overlapping_read_break")
        .put("re_name", re_bptr->name())
        .put("brk_start", tail_re_start_pos)
        .put("brk_end", re_end_pos);
    auto append_coordinates = [] (string& dest, Size_Type first, Size_Type last) {
        ostringstream tmp;
        tmp << ":(" << first << "," << last << ")";
        dest += tmp.str();
    };
    tail_re_bptr->name() = re_bptr->name();
    append_coordinates(re_bptr->name(), re_bptr->start(), re_end_pos);
    append_coordinates(tail_re_bptr->name(), tail_re_start_pos, re_bptr->end());
    // adjust start, len
    tail_re_bptr->chunk_cont().shift(-long(tail_re_start_pos));
    tail_re_bptr->start() = 0;
    tail_re_bptr->len() = re_bptr->end() - tail_re_start_pos;
    re_bptr->chunk_cont().shift(-long(re_bptr->start()));
    re_bptr->len() = re_end_pos - re_bptr->start();
    re_bptr->start() = 0;
    re_bptr->check();
    tail_re_bptr->check();
    LOG("Read_Entry", debug) << ptree("end");
    return make_pair(re_bptr, tail_re_bptr);
} // Read_Entry::split

pair< Read_Entry_BPtr, Read_Entry_BPtr >
Read_Entry::split(Read_Chunk_CBPtr rc_cbptr)
{
    LOG("Read_Entry", debug) << ptree("begin")
        .put("re_bptr", rc_cbptr->re_bptr().to_int())
        .put("rc_bptr", rc_cbptr.to_int());
    Read_Entry_BPtr re_bptr = rc_cbptr->re_bptr().unconst();
    ASSERT(re_bptr->is_unlinked());
    Read_Chunk_CBPtr next_rc_cbptr = re_bptr->chunk_cont().get_sibling(rc_cbptr, true, true);
    ASSERT(next_rc_cbptr);
    // create tail re, copy chunks
    Read_Entry_BPtr tail_re_bptr = Read_Entry_Fact::new_elem();
    auto rc_it = re_bptr->chunk_cont().iterator_to(*rc_cbptr.unconst());
    ASSERT(rc_it != re_bptr->chunk_cont().end());
    tail_re_bptr->chunk_cont().splice(re_bptr->chunk_cont(), next(rc_it), tail_re_bptr);
    // fix names and trim positions
    Size_Type tail_re_start_pos = rc_cbptr->get_r_end();
    LOG("Read_Entry", info) << ptree("non_overlapping_read_break")
        .put("re_name", re_bptr->name())
        .put("brk", tail_re_start_pos);
    auto append_coordinates = [] (string& dest, Size_Type first, Size_Type last) {
        ostringstream tmp;
        tmp << ":(" << first << "," << last << ")";
        dest += tmp.str();
    };
    tail_re_bptr->name() = re_bptr->name();
    append_coordinates(re_bptr->name(), re_bptr->start(), tail_re_start_pos);
    append_coordinates(tail_re_bptr->name(), tail_re_start_pos, re_bptr->end());
    tail_re_bptr->chunk_cont().shift(-long(tail_re_start_pos));
    tail_re_bptr->start() = 0;
    tail_re_bptr->len() = re_bptr->end() - tail_re_start_pos;
    re_bptr->chunk_cont().shift(-long(re_bptr->start()));
    re_bptr->len() = tail_re_start_pos - re_bptr->start();
    re_bptr->start() = 0;
    re_bptr->check();
    tail_re_bptr->check();
    LOG("Read_Entry", debug) << ptree("end");
    return make_pair(re_bptr, tail_re_bptr);
} // Read_Entry::split

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

ptree Read_Entry::to_ptree() const
{
    return ptree().put("name", name())
                  .put("len", len())
                  .put("chunk_cont", chunk_cont());
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
