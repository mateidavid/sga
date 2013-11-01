//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Mutation.hpp"

#include <iostream>
#include "indent.hpp"

using namespace std;


namespace MAC
{
    Mutation Mutation::cut(Size_Type base_offset, Size_Type alt_offset)
    {
        Mutation res;

        assert(base_offset <= _len);
        assert(alt_offset <= _seq_len);

        if (have_seq())
        {
            res = Mutation(_start + base_offset, _len - base_offset, _seq.substr(alt_offset));
            _len = base_offset;
            _seq_len = alt_offset;
            _seq.erase(alt_offset);
        }
        else
        {
            res = Mutation(_start + base_offset, _len - base_offset, _seq_len - alt_offset);
            _len = base_offset;
            _seq_len = alt_offset;
        }
        return res;
    }

    ostream& operator << (ostream& os, const Mutation& rhs)
    {
        os << "(Mutation &=" << (void*)&rhs
           << indent::inc << indent::nl << "start=" << (size_t)rhs.get_start()
           << ",len=" << (size_t)rhs.get_len()
           << ",seq_len=" << rhs.get_seq_len()
           << ",seq=" << rhs.get_seq()
           << indent::dec << indent::nl << ")";
        return os;
    }
}
