//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Mutation.hpp"

#include <iostream>
#include "../Util/Util.h"
#include "indent.hpp"

using namespace std;


namespace MAC
{

Mutation_CBPtr Mutation::cut(Size_Type base_offset, Size_Type alt_offset)
{
    Mutation_BPtr res;
    ASSERT(base_offset <= _len);
    ASSERT(alt_offset <= _seq_len);
    if (have_seq())
    {
        res = Mutation_Fact::new_elem(_start + base_offset, _len - base_offset, _seq.substr(alt_offset));
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

void Mutation::simplify(const Seq_Type& rf)
{
    ASSERT(rf.size() == _len);
    if (not have_seq())
    {
        return;
    }
    while (_len > 0 and _seq_len > 0 and rf[_len - 1] == _seq[_seq_len - 1])
    {
        --_len;
        --_seq_len;
        _seq.resize(_seq.size() - 1);
    }
    size_t i = 0;
    while (i < _len and i < _seq_len and rf[i] == _seq[i])
    {
        ++i;
    }
    if (i > 0)
    {
        _start += i;
        _len -= i;
        _seq_len -= i;
        _seq = _seq.substr(i);
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
                  .put("start", get_start())
                  .put("len", get_len())
                  .put("seq_len", get_seq_len())
                  .put("seq", get_seq());
}


} // namespace MAC
