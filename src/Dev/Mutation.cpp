//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Mutation.hpp"
#include "Contig_Entry.hpp"


namespace MAC
{

Mutation_CBPtr Mutation::cut(Size_Type base_offset, Size_Type alt_offset)
{
    Mutation_BPtr res;
    ASSERT(base_offset <= _len);
    ASSERT(alt_offset <= _seq_len);
    if (have_seq())
    {
        res = Mutation_Fact::new_elem(_start + base_offset, _len - base_offset, Seq_Type(_seq.substr(alt_offset)));
        _len = base_offset;
        _seq_len = alt_offset;
        _seq.erase(alt_offset);
    }
    else
    {
        res = Mutation_Fact::new_elem(_start + base_offset, _len - base_offset, _seq_len - alt_offset);
        _len = base_offset;
        _seq_len = alt_offset;
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
    auto p = ce_cbptr->mut_support(mut_cbptr);
    os << "mut " << setw(5) << left << mut_cbptr.to_int()
       << setw(5) << right << mut_cbptr->rf_start()
       << setw(5) << right << mut_cbptr->rf_len()
       << setw(5) << right << mut_cbptr->seq_len()
       << setw(5) << right << p.first.size()
       << setw(5) << right << p.second.size();
}


} // namespace MAC
