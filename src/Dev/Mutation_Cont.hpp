//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __MUTATION_CONT_HPP
#define __MUTATION_CONT_HPP

#include "shortcuts.hpp"
#include <boost/intrusive/itree.hpp>

#include "Mutation.hpp"

namespace MAC
{

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

    USING_ITERATORS(Base)
    //using Base::iterator_to;

    // check it is empty before deallocating
    ~Mutation_Cont() { ASSERT(size() == 0); }

    /** Create a Mutation container using mutations from a cigar string.
     * Pre: Cigar contains no 'M' operations (use disambiguate() first).
     * Post: Adjacent non-match operations are merged.
     * @param cigar Cigar object describing the match.
     * @param qr Query string; optional: if not given, Mutation objects contain alternate string lengths only.
     */
    Mutation_Cont(const Cigar& cigar, const std::string& qr = std::string());

    /** Get iterator to Mutation inside container. */
    const_iterator iterator_to(Mutation_CBPtr mut_cbptr) const { return Base::iterator_to(*mut_cbptr); }
    iterator iterator_to(Mutation_BPtr mut_bptr) { return Base::iterator_to(*mut_bptr); }

    /** Insert Mutation in container. */
    iterator insert(Mutation_BPtr mut_bptr) { return Base::insert(*mut_bptr); }

    /** Find an equivalent Mutation in container.
     * @param find_exact Bool; if true, look for that Mutation only; if false, look for equivalent Mutations as well.
     */
    const_iterator find(Mutation_CBPtr mut_cbptr, bool find_exact = false) const
    {
        const_iterator it;
        const_iterator it_end;
        for (std::tie(it, it_end) = Base::equal_range(*mut_cbptr); it != it_end; ++it)
        {
            if (&*it == mut_cbptr or (not find_exact and *mut_cbptr == *it))
            {
                return it;
            }
        }
        return end();
    }

    /** Erase Mutation from container. */
    void erase(Mutation_BPtr mut_bptr) { Base::erase(*mut_bptr); }

    /** Add Mutation to container; if an equivalent one already exists, use that one.
     * Note: Does not deallocate new Mutation when reusing old one.
     * @param mut_bptr Pointer to Mutation to add.
     * @return Pointer to Mutation inside container.
     */
    //Mutation_BPtr add_mut(Mutation_BPtr mut_bptr);

    /** Search for a Mutation that completely spans a given contig position.
     * @param c_pos Contig position, 0-based.
     * @return Pointer to Mutation, or NULL if none exists.
     */
    Mutation_CBPtr find_span_pos(Size_Type c_pos) const;

    /** Extract second half mutations, place them in new container.
     * Pre: No Mutation objects cross the cut.
     * @param c_brk Position of the cut.
     * @param mut_left_cbptr Insertion at c_pos to remain on the left of the cut, if any.
     * @return New Mutation_Cont.
     */
    Mutation_Cont split(Size_Type c_brk, Mutation_CBPtr mut_left_cbptr);

    /** Shift all Mutations in this container.
     * @param delta Signed integer value to add to all start points.
     */
    void shift(int delta)
    {
        for (auto mut_bref : *this)
        {
            Mutation_BPtr mut_bptr = &mut_bref;
            mut_bptr->shift(delta);
        }
        Base::implement_shift(delta);
    }

    /** Drop unused Mutation objects. */
    void drop_unused();
};

} // namespace MAC

#endif
