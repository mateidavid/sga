//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Mutation.hpp"
#include "Contig_Entry.hpp"


namespace MAC
{

Mutation_CBPtr Mutation::cut(Size_Type rf_offset, Size_Type qr_offset)
{
    Mutation_BPtr res;
    ASSERT(rf_offset <= _len);
    ASSERT(qr_offset <= _seq_len);
    if (have_seq())
    {
        res = Mutation_Fact::new_elem(_start + rf_offset, _len - rf_offset,
                                      Seq_Type(_seq.substr(qr_offset)));
        _len = rf_offset;
        _seq_len = qr_offset;
        _seq.resize(qr_offset);
    }
    else
    {
        res = Mutation_Fact::new_elem(_start + rf_offset, _len - rf_offset,
                                      _seq_len - qr_offset);
        _len = rf_offset;
        _seq_len = qr_offset;
    }
    return res;
}

void Mutation::simplify(const Seq_Proxy_Type& rf_seq)
{
    ASSERT(rf_seq.size() == _len);
    if (not have_seq())
    {
        return;
    }
    Seq_Proxy_Type qr_seq = seq();
    while (_len > 0 and _seq_len > 0 and rf_seq[_len - 1] == qr_seq[_seq_len - 1])
    {
        --_len;
        --_seq_len;
        //_seq.resize(_seq.size() - 1);
        qr_seq = qr_seq.substr(0, _seq_len);
    }
    size_t i = 0;
    while (i < _len and i < _seq_len and rf_seq[i] == _seq[i])
    {
        ++i;
    }
    if (i > 0)
    {
        _start += i;
        _len -= i;
        _seq_len -= i;
        //_seq = _seq.substr(i);
        qr_seq = qr_seq.substr(i);
    }
    if (qr_seq.size() != _seq.size())
    {
        _seq = qr_seq;
    }
}

void Mutation::to_stream(ostream& os, Mutation_CBPtr mut_cbptr, Contig_Entry_CBPtr ce_cbptr)
{
    ASSERT(mut_cbptr);
    ASSERT(ce_cbptr);
    set< Read_Chunk_CBPtr > qr_set;
    set< Read_Chunk_CBPtr > full_rf_set;
    tie(qr_set, full_rf_set, ignore) = ce_cbptr->mut_support(mut_cbptr);
    os << "mut " << setw(5) << left << mut_cbptr.to_int()
       << setw(5) << right << mut_cbptr->rf_start()
       << setw(5) << right << mut_cbptr->rf_len()
       << setw(5) << right << mut_cbptr->seq_len()
       << setw(5) << right << (int)mut_cbptr->uniq_flank(0)
       << setw(5) << right << (int)mut_cbptr->uniq_flank(1)
       << setw(5) << right << (int)mut_cbptr->copy_num(0)
       << setw(5) << right << (int)mut_cbptr->copy_num(1);
    Seq_Type flank_left(".");
    Seq_Type flank_right(".");
    if (mut_cbptr->uniq_flank(0) > 0 or mut_cbptr->uniq_flank(1) > 0)
    {
        auto flank_len = max(mut_cbptr->uniq_flank(0), mut_cbptr->uniq_flank(1));
        flank_left = ce_cbptr->substr(mut_cbptr->rf_start() - flank_len, flank_len);
        flank_right = ce_cbptr->substr(mut_cbptr->rf_end(), flank_len);
    }
    os << " " << flank_left << " (" << ce_cbptr->substr(mut_cbptr->rf_start(), mut_cbptr->rf_len())
       << "|" << mut_cbptr->seq() << ") " << flank_right << " ";

    os << setw(5) << right << full_rf_set.size()
       << setw(5) << right << qr_set.size();
    os << "   ";
    bool first = true;
    for (auto rc_cbptr : full_rf_set)
    {
        os << (not first? "," : "") << rc_cbptr.to_int();
        first = false;
    }
    os << "   ";
    first = true;
    for (auto rc_cbptr : qr_set)
    {
        os << (not first? "," : "") << rc_cbptr.to_int();
        first = false;
    }
}

void Mutation::save_strings(ostream& os, size_t& n_strings, size_t& n_bytes) const
{
    os.write(_seq.c_str(), _seq.size() + 1);
    ++n_strings;
    n_bytes += _seq.size() + 1;
}

void Mutation::init_strings()
{
    new (&_seq) string();
}

void Mutation::load_strings(istream& is, size_t& n_strings, size_t& n_bytes)
{
    getline(is, _seq, '\0');
    ++n_strings;
    n_bytes += _seq.size() + 1;
}

} // namespace MAC
