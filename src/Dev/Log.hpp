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
    : public std::ostringstream
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
        *this << "= " << facility << "." << int(msg_level) << ": ";
    }

    /** Destructor: dump buffer to output. */
    ~Log()
    {
        //_buffer << std::endl;
        std::cerr.write(this->str().c_str(), this->str().size());
    }

    /** Get and set the log level of the given facility. */
    static Level& level(const std::string& facility)
    {
        if (_level.count(facility) == 0)
        {
            _level[facility] = _default_level;
        }
        return _level[facility];
    }

    static Level& default_level() { return _default_level; }

private:
    static std::map< std::string, Level > _level;
    static Level _default_level;
}; // class Log

} // namespace detail

#define log_l(msg_level) \
    if (::detail::Log::msg_level > ::detail::Log::level(LOG_FACILITY)) ; \
    else ::detail::Log(LOG_FACILITY, ::detail::Log::msg_level)

#define log_i(int_msg_level) \
    if (static_cast< ::detail::Log::Level >(int_msg_level) > ::detail::Log::level(LOG_FACILITY)) ; \
    else ::detail::Log(LOG_FACILITY, static_cast< ::detail::Log::Level >(int_msg_level))


#endif
