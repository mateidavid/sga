#ifndef __VARIABLES_MAP_STRINGFIER_HPP
#define __VARIABLES_MAP_STRINGFIER_HPP

#include <boost/program_options.hpp>

#include "shortcuts.hpp"
#include "any_converter.hpp"
#include "ptree.hpp"


struct variables_map_converter
{
    // Options in the variable map are bundled into this type.
    struct option
    {
        const std::string& name;
        const boost::program_options::variable_value& value;
        const boost::any& converted_value;
    }; // struct option

    // Type of function to apply to each option in turn.
    using option_handler = std::function< void(const option&) >;

    // Apply handler to every option in the variables map,
    // converting values using the given any_converter.
    static void apply(const boost::program_options::variables_map& vm,
                      const any_converter& ac,
                      option_handler f)
    {
        for (const auto& p : vm)
        {
            option o = { p.first, p.second, ac.convert(p.second.value()) };
            f(o);
        }
    }

    // Dump options to stream.
    // PRE: ac must convert boost::any objects to string
    static void to_stream(const boost::program_options::variables_map& vm,
                          const any_converter& ac,
                          std::ostream& os)
    {
        auto f = [&os] (const option& o) {
            os << o.name << ": ";
            if (o.value.empty())
            {
                os << "(empty)";
            }
            else
            {
                if (o.value.defaulted())
                {
                    os << "(default) ";
                }
                if (o.converted_value.type() != typeid(std::string))
                {
                    std::cerr << "error: converted value is not a string" << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                os << boost::any_cast< const std::string& >(o.converted_value);
            }
            os << std::endl;
        };
        apply(vm, ac, f);
    }

    // Dump options to ptree.
    // PRE: ac must convert boost::any objects to either string or ptree
    static boost::property_tree::ptree to_ptree(const boost::program_options::variables_map& vm,
                                                const any_converter& ac)
    {
        ptree res;
        auto f = [&res] (const option& o) {
            if (o.value.empty())
            {
                res.put(o.name, "(empty)");
            }
            else
            {
                if (o.converted_value.type() == typeid(std::string))
                {
                    res.put(o.name, boost::any_cast< const std::string& >(o.converted_value));
                }
                else if (o.converted_value.type() == typeid(boost::property_tree::ptree))
                {
                    res.put(o.name, boost::any_cast< const boost::property_tree::ptree& >(o.converted_value));
                }
                else
                {
                    std::cerr << "error: converted value is neither a string nor a ptree" << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            }
        };
        apply(vm, ac, f);
        return res;
    }
}; // struct variables_map_converter


#endif
