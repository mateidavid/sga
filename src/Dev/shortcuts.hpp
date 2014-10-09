//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __SHORTCUTS_HPP
#define __SHORTCUTS_HPP

#include <type_traits>
#include <boost/type_traits/intrinsics.hpp>
#ifndef BOOST_IS_CONVERTIBLE
#define BOOST_IS_CONVERTIBLE(T, U) (std::is_convertible<T, U>::value)
#endif
#include <boost/utility.hpp>


#define DELETE_COPY_CTOR(_type) \
    _type(const _type&) = delete;
#define DELETE_MOVE_CTOR(_type) \
    _type(_type&&) = delete;
#define DELETE_DEF_CTOR(_type) \
    _type() = delete;
#define DELETE_COPY_ASOP(_type) \
    _type& operator = (const _type&) = delete;
#define DELETE_MOVE_ASOP(_type) \
    _type& operator = (_type&&) = delete;

#define DEFAULT_COPY_CTOR(_type) \
    _type(const _type&) = default;
#define DEFAULT_MOVE_CTOR(_type) \
    _type(_type&&) = default;
#define DEFAULT_DEF_CTOR(_type) \
    _type() = default;
#define DEFAULT_COPY_ASOP(_type) \
    _type& operator = (const _type&) = default;
#define DEFAULT_MOVE_ASOP(_type) \
    _type& operator = (_type&&) = default;

#define USING_ITERATORS(_type) \
    using typename _type::value_type; \
    using typename _type::iterator; \
    using typename _type::const_iterator; \
    using _type::begin; \
    using _type::end; \
    using _type::rbegin; \
    using _type::rend; \
    using _type::cbegin; \
    using _type::cend; \
    using _type::crbegin; \
    using _type::crend;

#define USING_ITERATOR_TO(_type) \
    using _type::iterator_to;

#define USING_SIZE_MEMBERS(_type) \
    using _type::size; \
    using _type::empty;

#define USING_STD_CONT(_type) \
    USING_ITERATORS(_type) \
    USING_SIZE_MEMBERS(_type)

#define USING_INTRUSIVE_CONT(_type) \
    USING_ITERATORS(_type) \
    USING_ITERATOR_TO(_type) \
    USING_SIZE_MEMBERS(_type)

#define T_ENABLE_IF(_cond) \
    bool _unused = true, typename boost::enable_if_c< _unused and (_cond), int >::type = 0

#define GETTER(_type, _pub_id, _pri_id) \
    const _type& _pub_id() const { return _pri_id; } \
    _type& _pub_id() { return _pri_id; }

#define BOOL_CONVERSION(_type, _access, _expr)  \
    private: \
        typedef void (_type::*_bool_type)() const; \
        void _bool_type_func() const {} \
    _access: \
        operator _bool_type () const { return (_expr)? &_type::_bool_type_func : nullptr; }


#endif
