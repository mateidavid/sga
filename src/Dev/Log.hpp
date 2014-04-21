//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __LOG_HPP
#define __LOG_HPP

#include <string>
#include <map>
#include <sstream>
#include <iostream>


namespace detail
{

class Log
{
public:
    enum Level
    {
        error = 0,
        warning,
        info,
        debug,
        debug1,
        debug2
    };

    /** Constructor: initialize buffer. */
    Log(const std::string& facility, Level msg_level)
    {
        _buffer << facility << "." << int(msg_level) << ": ";
    }

    /** Destructor: dump buffer to output. */
    ~Log()
    {
        _buffer << std::endl;
        std::cerr.write(_buffer.str().c_str(), _buffer.str().size());
    }

    /** Accumulate output in internal buffer. */
    template < typename T >
    Log& operator << (const T& obj)
    {
        _buffer << obj;
        return *this;
    }

    /** Get and set the log level of the given facility. */
    static Level& level(const std::string& facility)
    {
        if (_level.count(facility) == 0)
        {
            _level[facility] = error;
        }
        return _level[facility];
    }

private:
    std::ostringstream _buffer;

    static std::map< std::string, Level > _level;
}; // class Log

} // namespace detail

#define log_l(msg_level) \
    if (detail::Log::msg_level > detail::Log::level(LOG_FACILITY)) ; \
    else detail::Log(LOG_FACILITY, detail::Log::msg_level)

#define log_i(int_msg_level) \
    if (static_cast< detail::Log::Level >(int_msg_level) > detail::Log::level(LOG_FACILITY)) ; \
    else detail::Log(LOG_FACILITY, static_cast< detail::Log::Level >(int_msg_level))


#endif
