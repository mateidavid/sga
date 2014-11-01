//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Mutation.hpp"


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

bool operator == (const Mutation& lhs, const Mutation& rhs)
{
    return (lhs._start == rhs._start
            and lhs._len == rhs._len
            and lhs._seq_len == rhs._seq_len
            and lhs.have_seq() == rhs.have_seq()
            and (not lhs.have_seq() or lhs._seq == rhs._seq));
}

boost::property_tree::ptree Mutation::to_ptree() const
{
    return ptree().put("addr", (void*)this)
                  .put("rf_start", rf_start())
                  .put("rf_len", rf_len())
                  .put("seq_len", seq_len())
                  .put("seq", seq());
}


} // namespace MAC
