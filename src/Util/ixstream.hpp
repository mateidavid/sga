//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

//
// ixstream - istream/ifstream wrapper that works with text&gzip transparently
//


#ifndef __IXSTREAM_HPP
#define __IXSTREAM_HPP


#include <cstdlib>
#include <iostream>
#include <fstream>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>


class ixstream : public boost::iostreams::filtering_istream
{
public:
    ixstream() {}
    ixstream(const char* name) { open(name); }
    ixstream(const std::string& name) { open(name.c_str()); }
    ixstream(std::istream& is) { open(is); }

    void open(const char* name)
    {
        if (strncmp(name, "-", 2) == 0)
        {
            open(std::cin);
        }
        else
        {
            p_file = std::unique_ptr<std::ifstream>(new std::ifstream(name));
            if (!*p_file)
            {
                std::cerr << "error opening file [" << name << "]\n";
                exit(EXIT_FAILURE);
            }
            open(*p_file);
        }
    }

    void open(const std::string& name) { open(name.c_str()); }

    void open(std::istream& is)
    {
        int c = is.peek();
        if (is.good() and c == 31)
        {
            push(boost::iostreams::gzip_decompressor());
        }
        is.clear();
        push(is);
    }

private:
    // if necessary, store ifstream object, and auto-delete it when done
    std::unique_ptr<std::ifstream> p_file;

    // no copy constructor & copy-assignment operator:
    // declared private & *undefined*
    ixstream(const ixstream&);
    ixstream& operator = (const ixstream&);
};


#endif
