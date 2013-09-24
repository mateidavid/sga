//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

// Generate kmers out of a given string

#include "kmer_gen.hpp"

#include <cassert>

using namespace std;


size_t Kmer_gen::get_next()
{
    assert((end == start and start == 0)
           or (end == start + k and start != string::npos)
           or (end == start and start == string::npos));

    if (start == string::npos)
        // already done
        return start;

    if (end == start + k)
        // discard previously generated kmer
        ++start;

    // include next base
    ++end;
    while (end <= sq_p->size())
    {
        // consider string [start, ..., end-1]
        if ((qv_p != NULL
             and not qv_p->empty()
             and (*qv_p)[end - 1] - phred_offset < min_qv)
            or (allowed_chars_p != NULL
                and not allowed_chars_p->empty()
                and allowed_chars_p->find((*sq_p)[end - 1]) == string::npos))
        {
            // bad qv or bad base: skip all kmers containing this pos
            start = end;
            end = start + 1;
            continue;
        }
        assert(start < end);
        assert(end <= start + k);
        if (end == start + k)
            return start;
        ++end;
    }
    start = string::npos;
    end = string::npos;
    return start;
}
