//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __GLOBAL_HPP
#define __GLOBAL_HPP

#include <string>


struct global
{
    static std::string& program_name()
    {
        static std::string _program_name;
        return _program_name;
    }
    static std::string& assert_message()
    {
        static std::string _assert_message;
        return _assert_message;
    }
};


#endif
