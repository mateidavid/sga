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
    ixstream(const std::string& name) { open(name); }
    ixstream(std::istream& is) { open(is); }

    ixstream(const ixstream&) = delete;
    ixstream(ixstream&&) = delete;
    ixstream& operator = (const ixstream&) = delete;
    ixstream& operator = (ixstream&&) = delete;

    void open(const std::string& name)
    {
        if (name == "-")
        {
            open(std::cin);
        }
        else
        {
            p_file = std::unique_ptr<std::ifstream>(new std::ifstream(name.c_str()));
            if (!*p_file)
            {
                std::cerr << "error opening file [" << name << "]\n";
                exit(EXIT_FAILURE);
            }
            open(*p_file);
        }
    }

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
};


#endif
