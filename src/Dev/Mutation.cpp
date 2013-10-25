//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Mutation.hpp"

#include <iostream>

using namespace std;


namespace MAC
{
    ostream& operator << (ostream& os, const Mutation& rhs)
    {
        os << "(start=" << (size_t)rhs.get_start() << ",len=" << (size_t)rhs.get_len() << ",seq=" << rhs.get_seq() << ")";
        return os;
    }
}
