//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __MAC_FORWARD_HPP
#define __MAC_FORWARD_HPP

#include <string>
#include <functional>
#include <cassert>


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
    class Read_Chunk;
    class Read_Entry;
    class Contig_Entry;
}


/** General-purpose modifier for boost containers.
 * @param cont A container.
 * @param e_it Iterator to element inside container.
 * @param modifier Modifier to apply to element.
 */
template <class Container>
inline void modify_element(Container& cont, typename Container::iterator e_it, std::function<void(typename Container::value_type&)> modifier)
{
    assert(e_it != cont.end());
    bool success = cont.modify(e_it, modifier);
    assert(success);
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
    assert(e_cptr != NULL);
    modify_element<Container>(cont, cont.iterator_to(*e_cptr), modifier);
}


#endif
