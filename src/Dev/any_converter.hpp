#ifndef __ANY_CONVERTER_HPP
#define __ANY_CONVERTER_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <typeindex>
#include <functional>
#include <boost/any.hpp>


namespace detail
{

template < typename Dest_Type, typename Src_Type >
struct static_cast_converter;
template < typename Src_Type >
struct string_converter;

} // namespace detail

// Collection of converters that can be applied of boost::any objects.
// Every converter takes a boost::any as argument and produces a boost::any result.
// The class chooses the converter to apply based on the type_index of the argument.
class any_converter
{
public:
    // Add converter for type Src_Type.
    template < typename Src_Type >
    any_converter& add_converter(std::function< boost::any(const boost::any&) > f)
    {
        std::type_index ti = typeid(Src_Type);
        if (not f)
        {
            std::cerr << "error: type index [" << ti.name() << "]: empty converter" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (_converter_map.count(ti) > 0)
        {
            std::clog << "warning: type index [" << ti.name() << "]: already in the converter map" << std::endl;
        }
        _converter_map[ti] = std::move(f);
        return *this;
    }

    // Add regular function pointer as converter.
    template < typename Dest_Type, typename Src_Type >
    any_converter& add_converter(Dest_Type (*f)(const Src_Type&))
    {
        if (not f)
        {
            std::cerr << "error: null function pointer" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        return this->add_converter< Src_Type >([f] (const boost::any& v) {
                return boost::any(f(boost::any_cast< const Src_Type& >(v)));
            });
    }

    // Add static_cast converter.
    template < typename Dest_Type, typename Src_Type >
    any_converter& add_static_cast_converter()
    {
        return this->add_converter< Src_Type >(detail::static_cast_converter< Dest_Type, Src_Type >());
    }

    // Add string converter.
    template < typename Src_Type >
    any_converter& add_string_converter()
    {
        return this->add_converter< Src_Type >(detail::string_converter< Src_Type >());
    }

    // Convert boost::any argument using current converter table.
    boost::any convert(const boost::any& v) const
    {
        if (v.empty())
        {
            std::cerr << "error: v is empty" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        std::type_index ti = v.type();
        if (_converter_map.count(ti) == 0)
        {
            std::cerr << "error: type index [" << ti.name() << "]: not in the converter map" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        return _converter_map.at(ti)(v);
    }

private:
    std::map< std::type_index, std::function< boost::any(const boost::any&) > > _converter_map;
};  // class any_converter

namespace detail
{

// Default converter via static_cast.
template < typename Dest_Type, typename Src_Type >
struct static_cast_converter
{
    boost::any operator () (const boost::any& v) const
    {
        return static_cast< Dest_Type >(boost::any_cast< const Src_Type& >(v));
    }
};

// Default converter to string.
// NOTE: use operator << (ostream&, const Src_Type&)
template < typename Src_Type >
struct string_converter
{
    boost::any operator () (const boost::any& v) const
    {
        std::ostringstream os;
        os << boost::any_cast< const Src_Type& >(v);
        return os.str();
    }
};

} // namespace detail


#endif
