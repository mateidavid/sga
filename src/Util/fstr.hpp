#ifndef __FSTR_HPP
#define __FSTR_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>


/** Wrapper for an std::fstream object.
 *
 * The main purpose is to perform a check on open(), and if necessary, to abort the program
 * with an informative error message.
 */
class fstr
    : public std::fstream
{
public:
    fstr() {}
    fstr(const std::string& filename, std::ios_base::openmode mode = std::ios_base::in)
    {
        open(filename, mode);
    }

    fstr(const fstr&) = delete;
    fstr(fstr&&) = default;
    fstr& operator = (const fstr&) = delete;
    fstr& operator = (fstr&&) = default;

    void open(const std::string& filename, std::ios_base::openmode mode = std::ios_base::in)
    {
        std::fstream::open(filename, mode);
        if (this->fail()
            or ((mode & std::ios_base::in) and (this->peek(), this->fail())))
        {
            std::ostringstream oss;
            oss << "error in fstr::open(filename=\"" << filename << "\",mode=" << mode_to_string(mode) << ")";
            perror(oss.str().c_str());
            abort();
        }
        this->clear();
    }

    static std::string mode_to_string(std::ios_base::openmode mode)
    {
        static const std::vector< std::pair < std::ios_base::openmode, std::string > > v =
            {
                { std::ios_base::app, "app" },
                { std::ios_base::ate, "ate" },
                { std::ios_base::binary, "binary" },
                { std::ios_base::in, "in" },
                { std::ios_base::out, "out" },
                { std::ios_base::trunc, "trunc" }
            };
        std::string res;
        for (const auto& p : v)
        {
            if (mode & p.first)
            {
                res += (not res.empty()? "|" : "") + p.second;
            }
        }
        return res;
    }
};


#endif
