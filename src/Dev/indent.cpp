//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
// Inspired by:
// http://stackoverflow.com/questions/13842676/c-custom-output-stream-with-indentation
//-----------------------------------------------

#include "indent.hpp"

using namespace std;


namespace indent
{
    /** Helper function to get a storage index in a stream. */
    int get_indent_index()
    {
        /* ios_base::xalloc allocates indices for custom-storage locations. These indices are valid for all streams */
        static int index = ios_base::xalloc();
        return index;
    }

    ios_base& inc(ios_base& stream)
    {
        /* The iword(index) function gives a reference to the index-th custom storage location as a integer */
        stream.iword(get_indent_index())++;
        return stream;
    }

    ios_base& dec(ios_base& stream)
    {
        /* The iword(index) function gives a reference to the index-th custom storage location as a integer */
        if (stream.iword(get_indent_index()) > 0)
            stream.iword(get_indent_index())--;
        return stream;
    }

    ostream& tab(ostream& os)
    {
        int indent = os.iword(get_indent_index());
        while (indent)
        {
            os << '\t';
            indent--;
        }
        return os;
    }

    ostream& nl(ostream& os)
    {
        os << '\n' << tab;
        return os;
    }

    ostream& endl(ostream& os)
    {
        os << std::endl << tab;
        return os;
    }
}
