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

class Mutation_Chunk_Adapter
{
    // only constructed by Factory object
private:
    friend class bounded::Factory< Mutation_Chunk_Adapter >;

    /// Trivial default constructor
    DEFAULT_DEF_CTOR(Mutation_Chunk_Adapter);

    // disallow copy or move
    DELETE_COPY_CTOR(Mutation_Chunk_Adapter);
    DELETE_MOVE_CTOR(Mutation_Chunk_Adapter);
    DELETE_COPY_ASOP(Mutation_Chunk_Adapter);
    DELETE_MOVE_ASOP(Mutation_Chunk_Adapter);

    /** Constructor from Mutation and Read_Chunk pointers. */
    Mutation_Chunk_Adapter(Mutation_CBPtr mut_cbptr, Read_Chunk_CBPtr chunk_cbptr)
        : _mut_cbptr(mut_cbptr), _chunk_cbptr(chunk_cbptr) {}

    /** Check it is unlinked before destruction. */
    ~Mutation_Chunk_Adapter() { ASSERT(is_unlinked()); }

    /** Getters. */
public:
    GETTER(Mutation_CBPtr, mut_cbptr, _mut_cbptr)
    GETTER(Read_Chunk_CBPtr, chunk_cbptr, _chunk_cbptr)

    //friend std::ostream& operator << (std::ostream& os, const Mutation_Chunk_Adapter&);
    boost::property_tree::ptree to_ptree() const
    {
        return ptree().put("mut_ptr", mut_cbptr().to_ptree())
                      .put("chunk_ptr", chunk_cbptr().to_ptree());
    }

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
