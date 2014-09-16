//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __MUTATION_HPP
#define __MUTATION_HPP

#include <iostream>

#include "MAC_forward.hpp"
#include "Read_Chunk_Ptr_Cont.hpp"
#include "Cigar.hpp"
#include "../Util/Util.h"


namespace MAC
{

namespace detail {

struct Mutation_ITree_Node_Traits;

}

/** Holds information about a mutation from a base sequence.
 *
 * The class holds the start and span of the base sequence region affected by a mutation,
 * as well as the alternate sequence.
 */
class Mutation
{
private:
    // Can only be created by Factory object
    friend class bounded::Factory< Mutation >;

    /** Default constructor. */
    Mutation()
        : _start(0), _len(0), _seq_len(0) {}

    /** Constructor.
     * @param start Start of the mutation, i.e., length of base sequence prior to mutation.
     * @param len Length of the base sequence affected by the mutation.
     * @param seq Alternate sequence.
     */
    Mutation(Size_Type start, Size_Type len, const Seq_Type& seq = Seq_Type())
        : _seq(seq), _start(start), _len(len), _seq_len(seq.size()) {}

    /** Constructor.
     * @param start Start of the mutation, i.e., length of base sequence prior to mutation.
     * @param len Length of the base sequence affected by the mutation.
     * @param seq_len Length of alternate sequence.
     */
    Mutation(Size_Type start, Size_Type len, Size_Type seq_len)
        : _start(start), _len(len), _seq_len(seq_len) {}

    /** No copy constructor. */
    DELETE_COPY_CTOR(Mutation)
    /** Move constructor. */
    Mutation(Mutation&& rhs) : Mutation() { *this = std::move(rhs); }
public:
    /** No copy assignment. */
    DELETE_COPY_ASOP(Mutation)
    /** Move assignment: allow move only when unlinked. */
    Mutation& operator = (Mutation&& rhs)
    {
        if (this != &rhs)
        {
            ASSERT(is_unlinked() and rhs.is_unlinked());
            _seq = std::move(rhs._seq);
            _start = std::move(rhs._start);
            _len = std::move(rhs._len);
            _seq_len = std::move(rhs._seq_len);
            _chunk_ptr_cont = std::move(rhs._chunk_ptr_cont);
        }
        return *this;
    }

private:
    /** Destructor. */
    ~Mutation()
    {
        ASSERT(chunk_ptr_cont().empty());
        ASSERT(is_unlinked());
    }

public:
    /** @name Getters */
    /**@{*/
    Size_Type get_start() const { return _start; }
    Size_Type get_len() const { return _len; }
    Size_Type get_end() const { return _start + _len; }
    Size_Type get_seq_len() const { return _seq_len; }
    bool have_seq() const { return _seq.size() == _seq_len; }
    const Seq_Type& get_seq() const { return _seq; }
    Seq_Type& seq() { return _seq; }
    /**@}*/

    /** Find bounded pointer to this object.
     * Pre: Must be linked.
     */
    /*
    Mutation_CBPtr bptr_to() const
    {
        ASSERT(not is_unlinked());
        if (_parent->_l_child.raw() == this)
        {
            return _parent->_l_child;
        }
        if (_parent->_r_child.raw() == this)
        {
            return _parent->_r_child;
        }
        if (_parent->_parent.raw() == this)
        {
            return _parent->_parent;
        }
        ASSERT(false);
        return nullptr;
    }
    Mutation_BPtr bptr_to()
    {
        return static_cast< Mutation_BPtr >(const_cast< const Mutation* >(this)->bptr_to());
    }
    */

    /** @name Basic queries */
    /**@{*/
    bool is_ins() const { return _len == 0 and _seq_len > 0; }
    bool is_snp() const { return _len == 1 and _len == _seq_len; }
    bool is_del() const { return _len > 0 and _seq_len == 0; }
    bool is_empty() const { return _len == 0 and _seq_len == 0; }
    /**@}*/

    /** Merge with given Mutation.
     * Pre: Mutations must be unlinked, adjacent on rf, and appear in the same read chunks.
     * @param rhs Next Mutation.
     */
    /*
    void merge(Mutation&& rhs)
    {
        ASSERT(is_unlinked() and rhs.is_unlinked());
        ASSERT(_chunk_ptr_cont.size() == rhs._chunk_ptr_cont.size());
        if (is_empty())
        {
            *this = std::move(rhs);
        }
        else
        {
            ASSERT(get_end() == rhs.get_start());
            ASSERT(have_seq() == rhs.have_seq());
            _len += rhs._len;
            _seq_len += rhs._seq_len;
            _seq += rhs._seq;
        }
    }
    */

    /** Extend Mutation.
     * Pre: This Mutation is unlinked.
     * Pre: Mutation contains its alternate sequence.
     * @param start Reference position.
     * @param extra_len Extra reference length.
     * @param extra_seq Extra alternate sequence.
     */
    void extend(Size_Type start, Size_Type extra_len, const Seq_Type& extra_seq)
    {
        ASSERT(is_unlinked());
        ASSERT(have_seq());
        if (is_empty())
        {
            _start = start;
        }
        ASSERT(get_end() == start);
        _len += extra_len;
        _seq += extra_seq;
        _seq_len += extra_seq.size();
    }
    /** Extend Mutation; if empty, copy the given Mutation.
     * Pre: This Mutation is unlinked.
     * Pre: Both Mutations contain alternate sequences.
     * @param extra_mut_cbptr Extra mutation.
     */
    void extend(Mutation_CBPtr extra_mut_cbptr)
    {
        ASSERT(is_unlinked());
        ASSERT(have_seq() and extra_mut_cbptr->have_seq());
        if (is_empty())
        {
            _start = extra_mut_cbptr->get_start();
        }
        ASSERT(get_end() == extra_mut_cbptr->get_start());
        _len += extra_mut_cbptr->get_len();
        _seq += extra_mut_cbptr->get_seq();
        _seq_len += extra_mut_cbptr->get_seq_len();
    }

    /** Cut mutation at given offsets, allocate new Mutation to keep leftover.
     * @param base_offset Base offset, 0-based.
     * @param alt_offset Alternate sequence offset, 0-based.
     * @return The second part of the Mutation that was cut.
     */
    Mutation_CBPtr cut(Size_Type base_offset, Size_Type alt_offset);

    /** Simplify Mutation by dropping the ends of rf and qr if they match.
     * @param rf Reference sequence spanned by the mutation.
     */
    void simplify(const Seq_Type& rf);

    /** Reverse the mutation.
     * @param c_len The contig length.
     */
    void reverse(Size_Type c_len)
    {
        _start = c_len - (_start + _len);
        if (have_seq())
        {
            _seq = reverseComplement(_seq);
        }
    }

    /** Shift Mutation.
     * @param delta Signed integer value to add to start point.
     */
    template < typename delta_type >
    void shift(delta_type delta)
    {
        ASSERT(delta_type(_start) + delta >= 0);
        _start = Size_Type(delta_type(_start) + delta);
    }

    const Read_Chunk_Ptr_Cont& chunk_ptr_cont() const { return _chunk_ptr_cont; }
    Read_Chunk_Ptr_Cont& chunk_ptr_cont() { return _chunk_ptr_cont; }

    friend bool operator == (const Mutation&, const Mutation&);
    //friend std::ostream& operator << (std::ostream&, const Mutation&);
    boost::property_tree::ptree to_ptree() const;

private:
    /** Hooks for storage in intrusive interval trees inside Contig_Entry objects. */
    friend struct detail::Mutation_ITree_Node_Traits;
    bool is_unlinked() const { return not(_parent or _l_child or _r_child); }

    Seq_Type _seq;
    Size_Type _start;
    Size_Type _col_n_max_end;
    Read_Chunk_Ptr_Cont _chunk_ptr_cont;
    Mutation_BPtr _parent;
    Mutation_BPtr _l_child;
    Mutation_BPtr _r_child;
    uint32_t _len;
    uint32_t _seq_len;
}; // class Mutation

} // namespace MAC


#endif
