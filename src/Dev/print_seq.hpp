//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __PRINT_SEQ_HPP
#define __PRINT_SEQ_HPP

#include <iostream>


template <class T>
void print_seq(std::ostream& os, const T& t, const char *delim = ",")
{
    for (auto it = t.begin(); it != t.end(); ++it)
    {
        if (it != t.begin()) os << delim;
        os << *it;
    }
}

template <class T>
void print_ptr_seq(std::ostream& os, const T& t, const char *delim = ",")
{
    for (auto it = t.begin(); it != t.end(); ++it)
    {
        if (it != t.begin()) os << delim;
        os << *(*it);
    }
}


#endif
