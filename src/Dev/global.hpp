//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __GLOBAL_HPP
#define __GLOBAL_HPP

#include <iostream>
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
    static void assertion_failed(const std::string& expr, const std::string& msg,
                                 const std::string& function, const std::string& file, long line)
    {
        std::cerr << program_name() << ": "
                  << file << ":" << line << ": "
                  << function << ": "
                  << "Assertion '" << expr << "' failed";
        if (not msg.empty())
        {
            std::cerr << ": " << msg;
        }
        if (not assert_message().empty())
        {
            std::cerr << ": [" << assert_message() << "]";
        }
        std::cerr << std::endl;
        abort();
    }
}; // struct global

#undef ASSERT
#undef ASSERT_MSG

#if defined(DISABLE_ASSERTS)

#define ASSERT(expr) ((void)0)
#define ASSERT_MSG(expr, msg) ((void)0)

#else

#include <boost/config.hpp> // for BOOST_LIKELY
#include <boost/current_function.hpp>

#define ASSERT_MSG(expr, msg) (BOOST_LIKELY(!!(expr))? ((void)0): global::assertion_failed(#expr, msg, BOOST_CURRENT_FUNCTION, __FILE__, __LINE__))
#define ASSERT(expr) ASSERT_MSG(expr, "")

#endif

#endif
