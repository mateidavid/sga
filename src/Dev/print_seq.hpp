//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __PRINT_SEQ_HPP
#define __PRINT_SEQ_HPP

#include <iostream>

#include "indent.hpp"


template <class T, class D1, class D2 = const char*, class D3 = const char*>
void print_seq(std::ostream& os, const T& t, D1 delim, D2 delim_start = "", D3 delim_end = "")
{
    bool empty = true;
    for (auto it = t.begin(); it != t.end(); ++it)
    {
        if (it == t.begin())
            os << delim_start;
        else
            os << delim;
        os << *it;
        empty = false;
    }
    if (not empty)
        os << delim_end;
}

template <class T, class D1, class D2 = const char*, class D3 = const char*>
void print_ptr_seq(std::ostream& os, const T& t, D1 delim, D2 delim_start = "", D3 delim_end = "")
{
    bool empty = true;
    for (auto it = t.begin(); it != t.end(); ++it)
    {
        if (it == t.begin())
            os << delim_start;
        else
            os << delim;
        os << **it;
        empty = false;
    }
    if (not empty)
        os << delim_end;
}


#endif
