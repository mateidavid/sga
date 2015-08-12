#ifndef __ALLELE_CONT_HPP
#define __ALLELE_CONT_HPP

#include "Allele.hpp"

namespace MAC
{

namespace detail
{

struct Allele_List_Node_Traits
{
    typedef Allele node;
    typedef Allele_Fact fact_type;
    typedef fact_type::ptr_type node_ptr;
    typedef fact_type::const_ptr_type const_node_ptr;

    static node_ptr get_previous(const_node_ptr n) { return n->_previous; }
    static void set_previous(node_ptr n, node_ptr ptr) { n->_previous = ptr; }
    static node_ptr get_next(const_node_ptr n) { return n->_next; }
    static void set_next(node_ptr n, node_ptr ptr) { n->_next = ptr; }
}; // struct Allele_List_Node_Traits

struct Allele_List_Value_Traits
{
    typedef Allele value_type;
    typedef Allele_List_Node_Traits node_traits;
    typedef node_traits::node_ptr node_ptr;
    typedef node_traits::const_node_ptr const_node_ptr;
    typedef node_ptr pointer;
    typedef const_node_ptr const_pointer;
    typedef node_traits::fact_type::ref_type reference;
    typedef node_traits::fact_type::const_ref_type const_reference;
    typedef Allele_List_Value_Traits* value_traits_ptr;

    static const bi::link_mode_type link_mode = bi::safe_link;

    static node_ptr to_node_ptr (reference value) { return &value; }
    static const_node_ptr to_node_ptr (const_reference value) { return &value; }
    static pointer to_value_ptr(node_ptr n) { return n; }
    static const_pointer to_value_ptr(const_node_ptr n) { return n; }
}; // struct Allele_List_Value_Traits

} // namespace detail

class Allele_Cont
    : private bi::list< Allele,
                        bi::value_traits< detail::Allele_List_Value_Traits >,
                        bi::constant_time_size< false >,
                        bi::header_holder_type< bounded::Pointer_Holder< Allele > >
                      >
{
private:
    typedef bi::list< Allele,
                      bi::value_traits< detail::Allele_List_Value_Traits >,
                      bi::constant_time_size< false >,
                      bi::header_holder_type< bounded::Pointer_Holder< Allele > >
                      > Base;
public:
    // allow move only
    DEFAULT_DEF_CTOR(Allele_Cont);
    DELETE_COPY_CTOR(Allele_Cont);
    DEFAULT_MOVE_CTOR(Allele_Cont);
    DELETE_COPY_ASOP(Allele_Cont);
    DEFAULT_MOVE_ASOP(Allele_Cont);

    USING_INTRUSIVE_CONT(Base)
    using Base::check;

    // check it is empty when deallocating
    ~Allele_Cont() { ASSERT(empty()); }

    Base::size_type size() = delete;
    Base::size_type nonconst_size() const { return Base::size(); }

    /**
     * Find allele in container.
     */
    pair< const_iterator, unsigned > find(const Seq_Proxy_Type& seq) const
    {
        unsigned i = 0;
        for (auto it = begin(); it != end(); ++it)
        {
            if (it->seq() == seq)
            {
                return make_pair(it, i);
            }
        }
        return make_pair(end(), 0);
    }

    /**
     * Insert new sequence in container.
     */
    void insert(Seq_Type&& seq)
    {
        Allele_BPtr allele_bptr = Allele_Fact::new_elem(move(seq));
        Base::push_back(allele_bptr);
    }

}; // class Allele_Cont

}

#endif
