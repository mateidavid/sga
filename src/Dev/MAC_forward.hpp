//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __MAC_FORWARD_HPP
#define __MAC_FORWARD_HPP

#include <string>
#include <functional>

#include "common.hpp"
#include "factory.hpp"


/** Multi-Allelic Contig namespace. */
namespace MAC
{

/** Type for absolute and relative offsets in read and contig sequences. */
typedef size_t Size_Type;
/** Type holding a read sequence. */
typedef std::string Seq_Type;
/** Type holding a single base pair. */
typedef Seq_Type::value_type Symb_Type;

class Mutation;
class Mutation_Ptr_Node;
class Read_Chunk;
class Read_Chunk_Ptr_Node;
class Read_Entry;
class Contig_Entry;

typedef Factory< Mutation > Mutation_Fact;
typedef Mutation_Fact::ptr_type Mutation_BPtr;
typedef Mutation_Fact::const_ptr_type Mutation_CBPtr;
typedef Mutation_Fact::ref_type Mutation_BRef;
typedef Mutation_Fact::const_ref_type Mutation_CBRef;

typedef Factory< Mutation_Ptr_Node > Mutation_Ptr_Node_Fact;
typedef Mutation_Ptr_Node_Fact::ptr_type Mutation_Ptr_Node_BPtr;
typedef Mutation_Ptr_Node_Fact::const_ptr_type Mutation_Ptr_Node_CBPtr;
typedef Mutation_Ptr_Node_Fact::ref_type Mutation_Ptr_Node_BRef;
typedef Mutation_Ptr_Node_Fact::const_ref_type Mutation_Ptr_Node_CBRef;

typedef Factory< Read_Chunk > Read_Chunk_Fact;
typedef Read_Chunk_Fact::ptr_type Read_Chunk_BPtr;
typedef Read_Chunk_Fact::const_ptr_type Read_Chunk_CBPtr;
typedef Read_Chunk_Fact::ref_type Read_Chunk_BRef;
typedef Read_Chunk_Fact::const_ref_type Read_Chunk_CBRef;

typedef Factory< Read_Chunk_Ptr_Node > Read_Chunk_Ptr_Node_Fact;
typedef Read_Chunk_Ptr_Node_Fact::ptr_type Read_Chunk_Ptr_Node_BPtr;
typedef Read_Chunk_Ptr_Node_Fact::const_ptr_type Read_Chunk_Ptr_Node_CBPtr;
typedef Read_Chunk_Ptr_Node_Fact::ref_type Read_Chunk_Ptr_Node_BRef;
typedef Read_Chunk_Ptr_Node_Fact::const_ref_type Read_Chunk_Ptr_Node_CBRef;

typedef Factory< Read_Entry > Read_Entry_Fact;
typedef Read_Entry_Fact::ptr_type Read_Entry_BPtr;
typedef Read_Entry_Fact::const_ptr_type Read_Entry_CBPtr;
typedef Read_Entry_Fact::ref_type Read_Entry_BRef;
typedef Read_Entry_Fact::const_ref_type Read_Entry_CBRef;

typedef Factory< Contig_Entry > Contig_Entry_Fact;
typedef Contig_Entry_Fact::ptr_type Contig_Entry_BPtr;
typedef Contig_Entry_Fact::const_ptr_type Contig_Entry_CBPtr;
typedef Contig_Entry_Fact::ref_type Contig_Entry_BRef;
typedef Contig_Entry_Fact::const_ref_type Contig_Entry_CBRef;

}


/** General-purpose modifier for boost containers.
 * @param cont A container.
 * @param e_it Iterator to element inside container.
 * @param modifier Modifier to apply to element.
 */
template <class Container>
inline void modify_element(Container& cont, typename Container::iterator e_it, std::function<void(typename Container::value_type&)> modifier)
{
    ASSERT(e_it != cont.end());
    bool success = cont.modify(e_it, modifier);
    ASSERT(success);
}

/** General-purpose modifier for boost containers.
 * @param cont A container.
 * @param e Constant reference to element inside container.
 * @param modifier Modifier to apply to element.
 */
template <class Container>
inline void modify_element(Container& cont, const typename Container::value_type& e, std::function<void(typename Container::value_type&)> modifier)
{
    modify_element<Container>(cont, cont.iterator_to(e), modifier);
}

/** General-purpose modifier for boost containers.
 * @param cont A container.
 * @param e_cptr Constant pointer to element inside container.
 * @param modifier Modifier to apply to element.
 */
template <class Container>
inline void modify_element(Container& cont, const typename Container::value_type* e_cptr, std::function<void(typename Container::value_type&)> modifier)
{
    ASSERT(e_cptr != NULL);
    modify_element<Container>(cont, cont.iterator_to(*e_cptr), modifier);
}

template <class LHS_T, class RHS_T>
bool operator != (const LHS_T& lhs, const RHS_T& rhs)
{
    return !(lhs == rhs);
}


#endif
