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
private:
    // Can only be created by Factory object
    friend class Factory< Mutation >;

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
            _rcpn_cont = std::move(rhs._rcpn_cont);
        }
        return *this;
    }

private:
    /** Destructor. */
    ~Mutation() { ASSERT(is_unlinked()); }

public:
    /** @name Getters */
    /**@{*/
    Size_Type get_start() const { return _start; }
    Size_Type get_len() const { return _len; }
    Size_Type get_end() const { return _start + _len; }
    Size_Type get_seq_len() const { return _seq_len; }
    bool have_seq() const { return _seq.size() == _seq_len; }
    const Seq_Type& get_seq() const { return _seq; }
    /**@}*/

    /** Find bounded pointer to this object.
     * Pre: Must be linked.
     */
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
    void merge(Mutation&& rhs)
    {
        ASSERT(is_unlinked() and rhs.is_unlinked());
        ASSERT(_rcpn_cont.size() == rhs._rcpn_cont.size());
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
    bool is_unlinked() const { return not(_parent or _l_child or _r_child); }
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

class Mutation_Cont
 : private boost::intrusive::itree< Mutation_ITree_Value_Traits >
{
private:
    typedef boost::intrusive::itree< Mutation_ITree_Value_Traits > Base;
public:
    // allow move only
    DEFAULT_DEF_CTOR(Mutation_Cont)
    DELETE_COPY_CTOR(Mutation_Cont)
    DEFAULT_MOVE_CTOR(Mutation_Cont)
    DELETE_COPY_ASOP(Mutation_Cont)
    DEFAULT_MOVE_ASOP(Mutation_Cont)

    using Base::iterator;
    using Base::const_iterator;
    using Base::size;
    using Base::begin;
    using Base::end;
    using Base::iterator_to;

    // check it is empty before deallocating
    ~Mutation_Cont() { ASSERT(size() == 0); }

    /** Create a Mutation container using mutations from a cigar string.
     * Pre: Cigar contains no 'M' operations (use disambiguate() first).
     * Post: Adjacent non-match operations are merged.
     * @param cigar Cigar object describing the match.
     * @param qr Query string; optional: if not given, Mutation objects contain alternate string lengths only.
     */
    Mutation_Cont(const Cigar& cigar, const std::string& qr = std::string());

    /** Insert Mutation in container. */
    iterator insert(Mutation_BPtr mut_bptr) { return Base::insert(*mut_bptr); }

    /** Get iterator to Mutation inside container. */
    const_iterator iterator_to(Mutation_CBPtr mut_cbptr) const { return Base::iterator_to(*mut_cbptr); }
    iterator iterator_to(Mutation_BPtr mut_bptr) { return Base::iterator_to(*mut_bptr); }

    /** Find an equivalent Mutation in container. */
    const_iterator find(Mutation_CBPtr mut_cbptr) const
    {
        const_iterator it;
        const_iterator it_end;
        for (std::tie(it, it_end) = this->equal_range(*mut_cbptr); it != it_end; ++it)
        {
            if (*mut_cbptr == *it)
            {
                return it;
            }
        }
        return end();
    }
    iterator find(Mutation_BPtr mut_bptr)
    {
        return iterator_to(static_cast< Mutation_BRef >(*const_this(this)->find(mut_bptr)));
    }

    /** Erase Mutation from container. */
    void erase(Mutation_BPtr mut_bptr) { Base::erase(*mut_bptr); }

    /** Add Mutation to container; if an equivalent one already exists, use that one.
     * Note: Does not deallocate new Mutation when reusing old one.
     * @param mut_bptr Pointer to Mutation to add.
     * @return Pointer to Mutation inside container.
     */
    //Mutation_BPtr add_mut(Mutation_BPtr mut_bptr);
};


struct Mutation_Ptr_Node_Set_Node_Traits;

class Mutation_Ptr_Node
{
private:
    // Can only be created by Factory object
    friend class Factory< Mutation_Ptr_Node >;

    Mutation_Ptr_Node(Mutation_BPtr m_bptr = nullptr) : _m_bptr(m_bptr) {}

    // no copy or move
    DELETE_COPY_CTOR(Mutation_Ptr_Node)
    DELETE_MOVE_CTOR(Mutation_Ptr_Node)
public:
    DELETE_COPY_ASOP(Mutation_Ptr_Node)
    DELETE_MOVE_ASOP(Mutation_Ptr_Node)

    //Mutation_BPtr get() const { return _m_bptr; }
    operator Mutation_BPtr () const { return _m_bptr; }

    /** Find bounded pointer to this object.
     * Pre: Must be linked.
     */
    Mutation_Ptr_Node_CBPtr bptr_to() const
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
    Mutation_Ptr_Node_BPtr bptr_to()
    {
        return static_cast< Mutation_Ptr_Node_BPtr >(const_this(this)->bptr_to());
    }

private:
    const Mutation_BPtr _m_bptr;

    /** Hooks for storage in intrusive set inside Read_Chunk objects. */
    friend struct Mutation_Ptr_Node_Set_Node_Traits;
    Mutation_Ptr_Node_BPtr _parent;
    Mutation_Ptr_Node_BPtr _l_child;
    Mutation_Ptr_Node_BPtr _r_child;
    int _col;
    bool is_unlinked() const { return not( _parent or _l_child or _r_child); }
};

/** Comparator for storage in tree. */
struct Mutation_Ptr_Node_Compare
{
    bool operator () (const Mutation_Ptr_Node& lhs, const Mutation_Ptr_Node& rhs) const
    {
        Mutation_BPtr lhs_ptr(lhs);
        Mutation_BPtr rhs_ptr(rhs);
        ASSERT(lhs_ptr and rhs_ptr);
        return lhs_ptr->get_start() < rhs_ptr->get_start();
    }
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

class Mutation_Ptr_Cont
    : private boost::intrusive::set<
                Mutation_Ptr_Node,
                boost::intrusive::compare< Mutation_Ptr_Node_Compare >,
                boost::intrusive::value_traits< Mutation_Ptr_Node_Set_Value_Traits >
                >
{
private:
    typedef boost::intrusive::set<
                Mutation_Ptr_Node,
                boost::intrusive::compare< Mutation_Ptr_Node_Compare >,
                boost::intrusive::value_traits< Mutation_Ptr_Node_Set_Value_Traits >
                > Base;
public:
    // allow move only
    DEFAULT_DEF_CTOR(Mutation_Ptr_Cont)
    DELETE_COPY_CTOR(Mutation_Ptr_Cont)
    DEFAULT_MOVE_CTOR(Mutation_Ptr_Cont)
    DELETE_COPY_ASOP(Mutation_Ptr_Cont)
    DEFAULT_MOVE_ASOP(Mutation_Ptr_Cont)

    using Base::iterator;
    using Base::const_iterator;
    using Base::size;
    using Base::begin;
    using Base::end;

    // check it is empty when deallocating
    ~Mutation_Ptr_Cont() { ASSERT(size() == 0); }

    /** Construct a Mutation_Ptr_Cont from a Mutation_Cont.
     * Pre: Mutations may not touch in reference.
     */
    Mutation_Ptr_Cont(const Mutation_Cont& mut_cont)
    {
        Mutation_BPtr last_bptr = nullptr;
        for (const auto& mut_bref : mut_cont)
        {
            Mutation_BPtr crt_bptr(&mut_bref); // de-const
            if (last_bptr)
            {
                ASSERT(last_bptr->get_end() < crt_bptr->get_start());
            }
            auto tmp = insert(crt_bptr);
            ASSERT(tmp.second);
        }
    }

    /** Add a mutation pointer to this container.
     * Creates Mutation_Ptr_Node holding the given pointer, and inserts it in the container.
     */
    std::pair< Base::iterator, bool > insert(Mutation_BPtr mut_bptr)
    {
        Mutation_Ptr_Node_BPtr mp_bptr = Mutation_Ptr_Node_Fact::new_elem(mut_bptr);
        return Base::insert(*mp_bptr);
    }

    /** Find a given mutation pointer in this container. */
    Base::const_iterator find(Mutation_BPtr mut_bptr) const
    {
        for (const auto& mpn : *static_cast< const Base* >(this))
        {
            if (Mutation_BPtr(static_cast< const Mutation_Ptr_Node& >(mpn)) == mut_bptr)
            {
                return Base::iterator_to(mpn);
            }
        }
        return end();
    }
    Base::iterator find(Mutation_BPtr mut_bptr)
    {
        auto cit = const_this(this)->find(mut_bptr);
        return cit != end() ? Base::iterator_to(static_cast< Mutation_Ptr_Node_BRef >(*cit)) : end();
    }

    /** Add new mutation if old mutation exists. */
    void cond_insert(Mutation_BPtr old_mut_bptr, Mutation_BPtr new_mut_bptr)
    {
        if (find(old_mut_bptr) != end())
        {
            insert(new_mut_bptr);
        }
    }
};

}


#endif
