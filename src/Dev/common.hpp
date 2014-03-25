//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __COMMON_HPP
#define __COMMON_HPP

#include <iostream>


#define BOOST_ENABLE_ASSERT_HANDLER
#define BOOST_ASSERT_MSG_OSTREAM std::cerr
#include <boost/assert.hpp>
#define ASSERT BOOST_ASSERT
#define ASSERT_MSG BOOST_ASSERT_MSG

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


#endif
