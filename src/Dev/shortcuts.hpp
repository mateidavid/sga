//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __SHORTCUTS_HPP
#define __SHORTCUTS_HPP


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
    using _type::size; \
    using _type::iterator; \
    using _type::const_iterator; \
    using _type::begin; \
    using _type::end; \
    using _type::rbegin; \
    using _type::rend; \
    using _type::cbegin; \
    using _type::cend; \
    using _type::crbegin; \
    using _type::crend;


#endif
