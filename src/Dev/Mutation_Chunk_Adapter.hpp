//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __MUTATION_CHUNK_ADAPTER_HPP
#define __MUTATION_CHUNK_ADAPTER_HPP

#include "MAC_forward.hpp"


namespace MAC
{

namespace detail
{

struct Mutation_Ptr_List_Node_Traits;
struct Read_Chunk_Ptr_List_Node_Traits;

};

class Mutation_Chunk_Adapter
{
    // only constructed by Factory object
private:
    friend class Factory< Mutation_Chunk_Adapter >;

    // disallow copy or move
private:
    DEFAULT_DEF_CTOR(Mutation_Chunk_Adapter)
    DELETE_COPY_CTOR(Mutation_Chunk_Adapter)
    DELETE_MOVE_CTOR(Mutation_Chunk_Adapter)
public:
    DELETE_COPY_ASOP(Mutation_Chunk_Adapter)
    DELETE_MOVE_ASOP(Mutation_Chunk_Adapter)

    /** Constructor from Mutation and Read_Chunk pointers. */
private:
    Mutation_Chunk_Adapter(Mutation_CBPtr mut_cbptr, Read_Chunk_CBPtr chunk_cbptr)
        : _mut_cbptr(mut_cbptr), _chunk_cbptr(chunk_cbptr) {}

    /** Check it is unlinked before destruction. */
    ~Mutation_Chunk_Adapter() { ASSERT(is_unlinked()); }

    /** Getters. */
public:
    const Mutation_CBPtr& mut_cbptr() const { return _mut_cbptr; }
    Mutation_CBPtr& mut_cbptr() { return _mut_cbptr; }
    const Read_Chunk_CBPtr& chunk_cbptr() const { return _chunk_cbptr; }
    Read_Chunk_CBPtr& chunk_cbptr() { return _chunk_cbptr; }

    friend std::ostream& operator << (std::ostream& os, const Mutation_Chunk_Adapter&);

private:
    Mutation_CBPtr _mut_cbptr;
    Read_Chunk_CBPtr _chunk_cbptr;

    /** Hooks for storage in intrusive list inside Read_Chunk objects. */
    friend struct detail::Mutation_Ptr_List_Node_Traits;
    Mutation_Chunk_Adapter_BPtr _mut_ptr_previous;
    Mutation_Chunk_Adapter_BPtr _mut_ptr_next;

    /** Hooks for storage in intrusive list inside Mutation objects. */
    friend struct detail::Read_Chunk_Ptr_List_Node_Traits;
    Mutation_Chunk_Adapter_BPtr _chunk_ptr_previous;
    Mutation_Chunk_Adapter_BPtr _chunk_ptr_next;

    bool is_unlinked() const { return not(_mut_ptr_previous or _mut_ptr_next or _chunk_ptr_previous or _chunk_ptr_next); }
};

} // namespace MAC


#endif
