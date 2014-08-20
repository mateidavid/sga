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
#include <mutex>


namespace level_wrapper
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

} // namespace level_wrapper

class Logger
    : public std::ostringstream
{
public:
    static level_wrapper::level get_level(level_wrapper::level l) { return l; }
    static level_wrapper::level get_level(int i) { return static_cast< level_wrapper::level >(i); }
    static level_wrapper::level get_level(const std::string& s) { return level_from_string(s); }

    /** Constructor: initialize buffer. */
    Logger(const std::string& facility, level_wrapper::level msg_level)
    {
        *this << "= " << facility << "." << int(msg_level) << ": ";
    }
    /** Destructor: dump buffer to output. */
    ~Logger()
    {
        //_buffer << std::endl;
        std::clog.write(this->str().c_str(), this->str().size());
    }
    /** Produce l-value for output chaining. */
    std::ostream& l_value() { return *this; }

    static level_wrapper::level level_from_string(const std::string& s)
    {
        std::istringstream iss(s + "\n");
        int tmp_int = -1;
        iss >> tmp_int;
        if (iss.good())
        {
            return level_wrapper::level(tmp_int);
        }
        else
        {
            if (s == "error") return level_wrapper::error;
            else if (s == "warning") return level_wrapper::warning;
            else if (s == "info") return level_wrapper::info;
            else if (s == "debug") return level_wrapper::debug;
            else if (s == "debug1") return level_wrapper::debug1;
            else if (s == "debug2") return level_wrapper::debug2;
            else
            {
                std::cerr << "could not parse log level: " << s << "\n";
                std::exit(1);
            }
        }
    }

    static level_wrapper::level& last_level() {
        static thread_local level_wrapper::level _last_level = level_wrapper::error;
        return _last_level;
    }

    static level_wrapper::level get_default_level()
    {
        return default_level();
    }
    static void set_default_level(level_wrapper::level l)
    {
        mutex().lock();
        default_level() = l;
        mutex().unlock();
    }

    static level_wrapper::level get_facility_level(const std::string& facility)
    {
        return facility_level().count(facility) > 0? facility_level()[facility] : get_default_level();
    }
    static void set_facility_level(const std::string& facility, level_wrapper::level l)
    {
        mutex().lock();
        facility_level()[facility] = l;
        mutex().unlock();
    }

private:
    static level_wrapper::level& default_level()
    {
        static level_wrapper::level _default_level = level_wrapper::error;
        return _default_level;
    }

    static std::map< std::string, level_wrapper::level >& facility_level()
    {
        static std::map< std::string, level_wrapper::level > _facility_level;
        return _facility_level;
    }

    static std::mutex& mutex()
    {
        static std::mutex m;
        return m;
    }
}; // class logger


/**
 * log() macro
 *
 * log(level_spec) << message
 *   `level_spec` : integer, string, or logger level
 * 
 * log(facility, level_spec) << message
 *   `facility` : string
 *   `level_spec` : integer, string, or logger level
 * 
 * Log to `facility` (2nd form) or LOG_FACILITY (1st form) at logger level `level_spec`.
 */

#define log_2(facility, level_spec) \
    { using namespace level_wrapper; Logger::last_level() = Logger::get_level(level_spec); } \
    if (Logger::last_level() > Logger::get_facility_level(facility)) ; \
    else Logger(facility, Logger::last_level()).l_value()

#define log_1(level_spec) \
    log_2(LOG_FACILITY, level_spec)

#define _NARGS_AUX(_0, _1, _2, _3, _4, _5, _6, _7, _8, _9, ...) _9
#define _NARGS(...) _NARGS_AUX(__VA_ARGS__, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 0)

#define _log_aux_1(level_spec) log_1(level_spec)
#define _log_aux_2(facility, level_spec) log_2(facility, level_spec)

// we need 2-level indirection in order to trigger expansion after token pasting
// http://stackoverflow.com/questions/1597007/creating-c-macro-with-and-line-token-concatenation-with-positioning-macr
// http://stackoverflow.com/a/11763196/717706
#define __log_aux(N, ...) _log_aux_ ## N (__VA_ARGS__)
#define _log_aux(N, ...) __log_aux(N, __VA_ARGS__)

#define logger(...) _log_aux(_NARGS(__VA_ARGS__), __VA_ARGS__)


#endif
