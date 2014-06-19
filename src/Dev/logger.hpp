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


namespace logger
{

enum level
{
    error = 0,
    warning,
    info,
    debug,
    debug1,
    debug2
};

class log
    : public std::ostringstream
{
public:

    /** Constructor: initialize buffer. */
    log(const std::string& facility, level msg_level)
    {
        *this << "= " << facility << "." << int(msg_level) << ": ";
    }

    /** Destructor: dump buffer to output. */
    ~log()
    {
        //_buffer << std::endl;
        std::cerr.write(this->str().c_str(), this->str().size());
    }

    /** Get and set the log level of the given facility. */
    static level& facility_level(const std::string& facility)
    {
        if (_facility_level.count(facility) == 0)
        {
            _facility_level[facility] = _default_level;
        }
        return _facility_level[facility];
    }

    static level& default_level() { return _default_level; }

private:
    static std::map< std::string, level > _facility_level;
    static level _default_level;
}; // class logger

} // namespace log

#define log_l(msg_level) \
    if (logger::msg_level > logger::log::facility_level(LOG_FACILITY)) ; \
    else logger::log(LOG_FACILITY, logger::msg_level)

#define log_i(int_msg_level) \
    if (static_cast< logger::level >(int_msg_level) > logger::log::facility_level(LOG_FACILITY)) ; \
    else logger::log(LOG_FACILITY, static_cast< logger::level >(int_msg_level))


#endif
