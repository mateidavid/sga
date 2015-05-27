//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __GLOBAL_ASSERT_HPP
#define __GLOBAL_ASSERT_HPP

#include <iostream>
#include <string>

struct global_assert
{
    static std::string& prog_name()
    {
        static std::string _prog_name;
        return _prog_name;
    }
    static std::string& global_msg()
    {
        static std::string _global_msg;
        return _global_msg;
    }

    static void assertion_failed(const std::string& expr, const std::string& msg,
                                 const std::string& function, const std::string& file, long line)
    {
        std::cerr << prog_name() << ": "
                  << file << ":" << line << ": "
                  << function << ": "
                  << "Assertion '" << expr << "' failed";
        if (not msg.empty())
        {
            std::cerr << ": " << msg;
        }
        if (not global_msg().empty())
        {
            std::cerr << ": [" << global_msg() << "]";
        }
        std::cerr << std::endl;
        abort();
    }
}; // struct global_assert

#undef ASSERT
#undef ASSERT_MSG

#if defined(DISABLE_ASSERTS)

#define ASSERT(expr) ((void)0)
#define ASSERT_MSG(expr, msg) ((void)0)

#else

#include <boost/config.hpp> // for BOOST_LIKELY
#include <boost/current_function.hpp>

#define ASSERT_MSG(expr, msg) (BOOST_LIKELY(!!(expr))? ((void)0): global_assert::assertion_failed(#expr, msg, BOOST_CURRENT_FUNCTION, __FILE__, __LINE__))
#define ASSERT(expr) ASSERT_MSG(expr, "")

#endif

#endif
