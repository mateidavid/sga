//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __MUTATION_HPP
#define __MUTATION_HPP

#include <iostream>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/intrusive/itree.hpp>

#include "MAC_forward.hpp"
#include "Read_Chunk_forward.hpp"
#include "Cigar.hpp"
#include "../Util/Util.h"


namespace MAC
{

struct Mutation_ITree_Node_Traits;

/** Holds information about a mutation from a base sequence.
 *
 * The class holds the start and span of the base sequence region affected by a mutation,
 * as well as the alternate sequence.
 */
class Mutation
{
public:
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
    Mutation(const Mutation&) = delete;
    /** No copy assignment. */
    Mutation& operator = (const Mutation&) = delete;
    /** Move constructor. */
    Mutation(Mutation&& rhs) : Mutation() { *this = std::move(rhs); }
    /** Move assignment. */
    Mutation& operator = (Mutation&& rhs)
    {
        if (this != &rhs)
        {
            ASSERT(is_unlinked() and rhs.is_unlinked());
            _seq = std::move(rhs._seq);
            _start = std::move(rhs._start);
            _len = std::move(rhs._len);
            _seq_len = std::move(rhs._seq_len);
            _rcpn_cont = std::move(rhs._rcpn_cont);
        }
        return *this;
    }
    /** Destructor. */
    ~Mutation() { ASSERT(is_unlinked()); }

    /** @name Getters */
    /**@{*/
    Size_Type get_start() const { return _start; }
    Size_Type get_len() const { return _len; }
    Size_Type get_end() const { return _start + _len; }
    Size_Type get_seq_len() const { return _seq_len; }
    bool have_seq() const { return _seq.size() == _seq_len; }
    const Seq_Type& get_seq() const { return _seq; }
    /**@}*/

    /** @name Basic queries */
    /**@{*/
    bool is_ins() const { return _len == 0 and _seq_len > 0; }
    bool is_snp() const { return _len == 1 and _len == _seq_len; }
    bool is_del() const { return _len > 0 and _seq_len == 0; }
    bool is_empty() const { return _len == 0 and _seq_len == 0; }
    /**@}*/

    /** Merge with given Mutation.
     * Pre: Mutations must be adjacent on rf.
     * @param rhs Next Mutation.
     */
    void merge(const Mutation& rhs)
    {
        if (is_empty())
        {
            _seq = rhs._seq;
            _start = rhs._start;
            _len = rhs._len;
            _seq_len = rhs._seq_len;
        }
        else
        {
            ASSERT(get_end() == rhs.get_start());
            ASSERT(have_seq() == rhs.have_seq());
            _len += rhs._len;
            _seq_len += rhs._seq_len;
            _seq += rhs._seq;
        }
        //TODO: merge chunk ptr lists
    }

    /** Cut mutation at given offsets.
     * @param base_offset Base offset, 0-based.
     * @param alt_offset Alternate sequence offset, 0-based.
     * @return The part of the original mutation that was cut from this object.
     */
    Mutation cut(Size_Type base_offset, Size_Type alt_offset);

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
            _seq = reverseComplement(_seq);
    }

    /** Change Mutation base by adding a prefix of the given length. */
    void add_base_prefix(Size_Type len) { _start += len; }

    const Read_Chunk_Ptr_Node_Cont& chunk_ptr_cont() const { return _rcpn_cont; }
    Read_Chunk_Ptr_Node_Cont& chunk_ptr_cont() { return _rcpn_cont; }

    friend bool operator == (const Mutation&, const Mutation&);
    friend std::ostream& operator << (std::ostream&, const Mutation&);

private:
    Seq_Type _seq;
    Size_Type _start;
    Size_Type _len;
    Size_Type _seq_len;

    Read_Chunk_Ptr_Node_Cont _rcpn_cont;

    /** Hooks for storage in intrusive interval trees inside Contig_Entry objects. */
    friend struct Mutation_ITree_Node_Traits;
    Mutation_BPtr _parent;
    Mutation_BPtr _l_child;
    Mutation_BPtr _r_child;
    Size_Type _max_end;
    int _col;
    bool is_unlinked() const { return _parent or _l_child or _r_child; }
};

struct Mutation_ITree_Node_Traits
{
    typedef Holder< Mutation > node;
    typedef Factory< Mutation > fact_type;
    typedef typename fact_type::ptr_type node_ptr;
    typedef typename fact_type::const_ptr_type const_node_ptr;
    typedef int color;
    typedef Size_Type key_type;

    static node_ptr get_parent(const_node_ptr n) { return n->_parent; }
    static void set_parent(node_ptr n, node_ptr ptr) { n->_parent = ptr; }
    static node_ptr get_left(const_node_ptr n) { return n->_l_child; }
    static void set_left(node_ptr n, node_ptr ptr) { n->_l_child = ptr; }
    static node_ptr get_right(const_node_ptr n) { return n->_r_child; }
    static void set_right(node_ptr n, node_ptr ptr) { n->_r_child = ptr; }
    static color get_color(const_node_ptr n) { return n->_col; }
    static void set_color(node_ptr n, color c) { n->_col = c ; }
    static color black() { return 0; }
    static color red() { return 1; }
    static key_type get_max_end(const_node_ptr n) { return n->_max_end; }
    static void set_max_end(node_ptr n, key_type k) { n->_max_end = k ; }
};

struct Mutation_ITree_Value_Traits
{
    typedef Mutation value_type;
    typedef Mutation_ITree_Node_Traits node_traits;
    typedef typename node_traits::key_type key_type;
    typedef typename node_traits::node_ptr node_ptr;
    typedef typename node_traits::const_node_ptr const_node_ptr;
    typedef node_ptr pointer;
    typedef const_node_ptr const_pointer;
    typedef typename node_traits::fact_type::ref_type reference;
    typedef typename node_traits::fact_type::const_ref_type const_reference;

    static const boost::intrusive::link_mode_type link_mode = boost::intrusive::normal_link;

    static node_ptr to_node_ptr (reference value) { return &value; }
    static const_node_ptr to_node_ptr (const_reference value) { return &value; }
    static pointer to_value_ptr(node_ptr n) { return n; }
    static const_pointer to_value_ptr(const_node_ptr n) { return n; }
    static key_type get_start(const_pointer n) { return n->get_start(); }
    static key_type get_start(const value_type* n) { return n->get_start(); }
    static key_type get_end(const_pointer n) { return n->get_end(); }
    static key_type get_end(const value_type* n) { return n->get_end(); }
};

typedef boost::intrusive::itree< Mutation_ITree_Value_Traits > Mutation_Cont;

struct Mutation_Ptr_Node_Set_Node_Traits;

class Mutation_Ptr_Node
{
public:
    Mutation_Ptr_Node(Mutation_BPtr m_bptr) : _m_bptr(m_bptr) {}

    Mutation_BPtr get() const { return _m_bptr; }

private:
    const Mutation_BPtr _m_bptr;

    friend struct Mutation_Ptr_Node_Set_Node_Traits;
    Mutation_Ptr_Node_BPtr _parent;
    Mutation_Ptr_Node_BPtr _l_child;
    Mutation_Ptr_Node_BPtr _r_child;
    int _col;
};

struct Mutation_Ptr_Node_Set_Node_Traits
{
    typedef Holder< Mutation_Ptr_Node > node;
    typedef Mutation_Ptr_Node_Fact fact_type;
    typedef fact_type::ptr_type node_ptr;
    typedef fact_type::const_ptr_type const_node_ptr;
    typedef int color;

    static node_ptr get_parent(const_node_ptr n) { return n->_parent; }
    static void set_parent(node_ptr n, node_ptr ptr) { n->_parent = ptr; }
    static node_ptr get_left(const_node_ptr n) { return n->_l_child; }
    static void set_left(node_ptr n, node_ptr ptr) { n->_l_child = ptr; }
    static node_ptr get_right(const_node_ptr n) { return n->_r_child; }
    static void set_right(node_ptr n, node_ptr ptr) { n->_r_child = ptr; }
    static color get_color(const_node_ptr n) { return n->_col; }
    static void set_color(node_ptr n, color c) { n->_col = c ; }
    static color black() { return 0; }
    static color red() { return 1; }
};

struct Mutation_Ptr_Node_Set_Value_Traits
{
    typedef Mutation_Ptr_Node value_type;
    typedef Mutation_Ptr_Node_Set_Node_Traits node_traits;
    typedef node_traits::node_ptr node_ptr;
    typedef node_traits::const_node_ptr const_node_ptr;
    typedef node_ptr pointer;
    typedef const_node_ptr const_pointer;
    typedef node_traits::fact_type::ref_type reference;
    typedef node_traits::fact_type::const_ref_type const_reference;

    static const boost::intrusive::link_mode_type link_mode = boost::intrusive::normal_link;

    static node_ptr to_node_ptr (reference value) { return &value; }
    static const_node_ptr to_node_ptr (const_reference value) { return &value; }
    static pointer to_value_ptr(node_ptr n) { return n; }
    static const_pointer to_value_ptr(const_node_ptr n) { return n; }
};

typedef boost::intrusive::set<
    Mutation_Ptr_Node,
    boost::intrusive::value_traits< Mutation_Ptr_Node_Set_Value_Traits >
> Mutation_Ptr_Node_Cont;

/** Create a set of mutations to a reference string based on a cigar object.
 * Pre: Cigar contains no 'M' operations (use disambiguate() first).
 * Post: Adjacent non-match operations are merged.
 * @param cigar Cigar object describing the match.
 * @param qr Query string; optional: if not given, Mutation objects contain alternate string lengths only.
 * @return Container of Mutation objects.
 */
//std::shared_ptr< Mutation_Cont > make_mutations_from_cigar(const Cigar& cigar, const std::string& qr = std::string());

/** Add Mutation to container, use existing Mutation if it already exists.
 * @param mut_cont Mutation container.
 * @param mut Mutation to add.
 * @return Pointer to Mutation inside container.
 */
//Mutation_CPtr add_mut_to_cont(Mutation_Cont& mut_cont, const Mutation& mut);

}


#endif
