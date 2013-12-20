//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __GLOBALS_HPP
#define __GLOBALS_HPP

#ifdef __DEFINE_GLOBALS
#define GLOBAL_INIT(_type, _id, _init_val) extern _type _id; _type _id = _init_val
#define GLOBAL(_type, _id) extern _type _id; _type _id
#else
#define GLOBAL_INIT(_type, _id, _init_val) extern _type _id
#define GLOBAL(_type, _id) extern _type _id
#endif

#define CONCAT(...) __VA_ARGS__


namespace global
{
    GLOBAL(const char *, prog_name);
    GLOBAL(const char *, assert_message);
}


#endif
