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

#include <string>


namespace global
{
    GLOBAL(std::string, program_name);
    GLOBAL(std::string, assert_message);
    GLOBAL_INIT(bool, merge_contigs_at_each_step, false);
    GLOBAL_INIT(bool, merge_contigs_at_end, false);
    GLOBAL_INIT(bool, print_graph, false);
    GLOBAL_INIT(bool, print_graph_each_step, false);
    GLOBAL_INIT(bool, progress_graph_op, false);
    GLOBAL_INIT(int, progress_count, 0);
    GLOBAL_INIT(size_t, unmap_trigger_len, 9);
}


#endif
