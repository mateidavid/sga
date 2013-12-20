#include "common.hpp"

#include <iostream>
#include <cstdlib>

#include "globals.hpp"


namespace boost
{
    void assertion_failed(char const * expr, char const * function, char const * file, long line)
    {
        BOOST_ASSERT_MSG_OSTREAM << (global::prog_name != NULL? global::prog_name : "?") << ": "
            << file << ':' << line << ": "
            << function << ": "
            << "Assertion '" << expr << "' failed";
        if (global::assert_message != NULL)
        {
            BOOST_ASSERT_MSG_OSTREAM << ": [" << global::assert_message << ']';
        }
        BOOST_ASSERT_MSG_OSTREAM << std::endl;
        std::abort();
    }
    void assertion_failed_msg(char const * expr, char const * msg, char const * function, char const * file, long line)
    {
        BOOST_ASSERT_MSG_OSTREAM << (global::prog_name != NULL? global::prog_name : "?") << ": "
            << file << ':' << line << ": "
            << function << ": "
            << "Assertion '" << expr << "' failed: " << msg;
        if (global::assert_message != NULL)
        {
            BOOST_ASSERT_MSG_OSTREAM << ": [" << global::assert_message << ']';
        }
        BOOST_ASSERT_MSG_OSTREAM << std::endl;
        std::abort();
    }
}
