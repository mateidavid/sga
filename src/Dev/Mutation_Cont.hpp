//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __MUTATION_CONT_HPP
#define __MUTATION_CONT_HPP

#include "Mutation.hpp"


namespace MAC
{

namespace detail
{

struct Mutation_Set_Node_Traits
{
    typedef Mutation node;
    typedef Mutation_Fact fact_type;
    typedef typename fact_type::ptr_type node_ptr;
    typedef typename fact_type::const_ptr_type const_node_ptr;
    typedef uint8_t color;
    typedef Size_Type key_type;

    static node_ptr get_parent(const_node_ptr n) { return n->_parent; }
    static void set_parent(node_ptr n, node_ptr ptr) { n->_parent = ptr; }
    static node_ptr get_left(const_node_ptr n) { return n->_l_child; }
    static void set_left(node_ptr n, node_ptr ptr) { n->_l_child = ptr; }
    static node_ptr get_right(const_node_ptr n) { return n->_r_child; }
    static void set_right(node_ptr n, node_ptr ptr) { n->_r_child = ptr; }
    static color black() { return 0; }
    static color red() { return 1; }
    static color get_color(const_node_ptr n) { return n->_color; }
    static void set_color(node_ptr n, color c) { n->_color = c; }
}; // struct Mutation_Set_Node_Traits

struct Mutation_Set_Value_Traits
{
    typedef Mutation value_type;
    typedef Mutation_Set_Node_Traits node_traits;
    typedef typename node_traits::key_type key_type;
    typedef typename node_traits::node_ptr node_ptr;
    typedef typename node_traits::const_node_ptr const_node_ptr;
    typedef node_ptr pointer;
    typedef const_node_ptr const_pointer;
    typedef typename node_traits::fact_type::ref_type reference;
    typedef typename node_traits::fact_type::const_ref_type const_reference;
    typedef Mutation_Set_Value_Traits* value_traits_ptr;

    static const bi::link_mode_type link_mode = bi::safe_link;

    static node_ptr to_node_ptr (reference value) { return &value; }
    static const_node_ptr to_node_ptr (const_reference value) { return &value; }
    static pointer to_value_ptr(node_ptr n) { return n; }
    static const_pointer to_value_ptr(const_node_ptr n) { return n; }
}; // struct Mutation_Set_Value_Traits

} // namespace detail

class Mutation_Cont
    : private bi::multiset< Mutation,
                            bi::value_traits< detail::Mutation_Set_Value_Traits >,
                            bi::constant_time_size< false >,
                            bi::header_holder_type< bounded::Pointer_Holder< Mutation > > >
{
private:
    typedef bi::multiset< Mutation,
                          bi::value_traits< detail::Mutation_Set_Value_Traits >,
                          bi::constant_time_size< false >,
                          bi::header_holder_type< bounded::Pointer_Holder< Mutation > > > Base;
public:
    // allow move only
    DEFAULT_DEF_CTOR(Mutation_Cont);
    DELETE_COPY_CTOR(Mutation_Cont);
    DEFAULT_MOVE_CTOR(Mutation_Cont);
    DELETE_COPY_ASOP(Mutation_Cont);
    DEFAULT_MOVE_ASOP(Mutation_Cont);

    USING_INTRUSIVE_CONT(Base)
    using Base::check;
    using Base::lower_bound;
    friend class Graph;

    // check it is empty before deallocating
    ~Mutation_Cont() { ASSERT(empty()); }

    Base::size_type size() = delete;
    Base::size_type nonconst_size() const { return Base::size(); }

    /**
     * Create a Mutation container using mutations from a cigar string.
     * Pre: Cigar contains no 'M' operations (use disambiguate() first).
     * Post: Adjacent non-match operations are merged.
     * @param cigar Cigar object describing the match.
     * @param qr Query string.
     */
    static Mutation_Cont make_from_cigar(const Cigar& cigar, Contig_Entry_CBPtr ce_cbptr, const Seq_Proxy_Type& qr);

    /// Insert Mutation in container.
    iterator insert(Mutation_BPtr mut_bptr) { return Base::insert(*mut_bptr); }
    /// Erase Mutation from container.
    void erase(Mutation_CBPtr mut_cbptr) { Base::erase(iterator_to(*mut_cbptr)); }

    /**
     * Find an equivalent Mutation in container.
     * @param find_exact Bool; if true, look for that Mutation only; if false,
     * look for equivalent Mutations as well.
     * @return An exact (if find_exact is true) or equivalent Mutation in this
     * container; or NULL if none exists.
     *
    Mutation_CBPtr find(Mutation_CBPtr mut_cbptr, bool find_exact) const
    {
        const_iterator it;
        const_iterator it_end;
        for (std::tie(it, it_end) = Base::equal_range(*mut_cbptr); it != it_end; ++it)
        {
            if (&*it == mut_cbptr or (not find_exact and *mut_cbptr == *it))
            {
                return &*it;
            }
        }
        return nullptr;
    }
    */

    /**
     * Get equivalent Mutation from container;
     * if one exists, the given Mutation is deallocated;
     * if one doesn't exist, the given Mutation is added and returned.
     * Pre: Mutation is unlinked.
     * Pre: Mutation has no chunk pointers.
     *
    Mutation_CBPtr find_equiv_or_add(Mutation_BPtr mut_bptr)
    {
        ASSERT(mut_bptr->chunk_ptr_cont().empty());
        Mutation_CBPtr equiv_mut_cbptr = find(mut_bptr, false);
        if (equiv_mut_cbptr)
        {
            // an equivalent Mutation exists; we deallocate new one
            Mutation_Fact::del_elem(mut_bptr);
            return equiv_mut_cbptr;
        }
        else
        {
            // an equivalent Mutation does not exist; we add and return this one
            insert(mut_bptr);
            return mut_bptr;
        }
    }
    */

    /**
     * Search for a Mutation that completely spans a given contig position.
     * @param c_pos Contig position, 0-based.
     * @return Pointer to Mutation, or NULL if none exists.
     *
    Mutation_CBPtr find_span_pos(Size_Type c_pos) const;
    */

    /**
     * Move second half mutations from given container to this one.
     * Pre: No Mutation objects cross the cut.
     * @param other_cont Source container.
     * @param c_brk Position of the cut.
     * @param mut_left_cbptr Insertion at c_brk to remain on the left of the cut, if any.
     */
    void splice(Mutation_Cont& lhs_cont, const_iterator rhs_mut_cit)
    { splice(begin(), lhs_cont, rhs_mut_cit, lhs_cont.end()); }

    /**
     * Shift all Mutations in this container.
     * @param delta Signed integer value to add to all start points.
     */
    template < typename delta_type >
    void shift(delta_type delta)
    { for (auto mut_bptr : *this | referenced) mut_bptr->shift(delta); }

    /// Drop unused Mutation objects.
    //void drop_unused();

    /// Dispose of Mutation object
    static void dispose(Mutation_BPtr mut_bptr) { Mutation_Fact::del_elem(mut_bptr); }
    /// Clear the container and deallocate Mutation objects.
    void clear_and_dispose() { Base::clear_and_dispose(&dispose); }

    /** Reverse Mutations in this container (and their order). */
    void reverse_mutations()
    {
        Mutation_Cont new_mut_cont;
        Base::clear_and_dispose([&] (Mutation_BPtr mut_bptr)
        {
            mut_bptr->reverse();
            new_mut_cont.insert(mut_bptr);
        });
        *this = move(new_mut_cont);
    }

    /**
     * Move Mutations from given container into this one.
     * @param other_cont Container to clear.
     */
    void splice(iterator pos, Mutation_Cont& other_cont, const_iterator it_start, const_iterator it_end)
    {
        static_cast< Base& >(other_cont).erase_and_dispose(
            it_start, it_end,
            [&] (Mutation_BPtr mut_bptr)
            {
                Base::insert(pos, *mut_bptr);
            });
    }

    /**
     * Acquire Mutations from given container.
     * In contrast to splice(), this method looks for common Mutations. If such Mutations are found,
     * the existing copy is used, and the one in other_cont is deallocated. MCA-s are also adjusted.
     * @param other_cont Container to clear.
     *
    void merge(Mutation_Cont& other_cont)
    {
        if (empty())
        {
            *this = std::move(other_cont);
        }
        else
        {
            static_cast< Base& >(other_cont).clear_and_dispose([&] (Mutation_BPtr mut_bptr)
            {
                auto equiv_mut_bptr = find(mut_bptr, false).unconst();
                if (not equiv_mut_bptr)
                {
                    insert(mut_bptr);
                }
                else
                {
                    // an equivalent mutation exists in this container
                    // adjust mca-s to have mut_ptr point to existing mutation
                    for (auto mca_bptr : mut_bptr->chunk_ptr_cont() | referenced)
                    {
                        mca_bptr->mut_cbptr() = equiv_mut_bptr;
                    }
                    // merge the chunk_ptr containers
                    equiv_mut_bptr->chunk_ptr_cont().splice(mut_bptr->chunk_ptr_cont());
                    // deallocate mutation
                    Mutation_Fact::del_elem(mut_bptr);
                }
            });
        }
    }
    */
};

} // namespace MAC

#endif
