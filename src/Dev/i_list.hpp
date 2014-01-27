#ifndef __I_LIST_HPP
#define __I_LIST_HPP

#include <limits>
#include <boost/iterator/iterator_facade.hpp>

#include "common.hpp"


namespace detail
{
    /** Iterator base class, shared by const and non-const versions.
     * @param Traits Struct giving typedefs and static methods to access intrusive list fields.
     * @param dir_tag Bool; true iff forward traversal.
     */
    template <class Traits, bool dir_tag>
    class I_List_Iterator_Base : public boost::iterator_facade<
        I_List_Iterator_Base< Traits, dir_tag >,
        typename Traits::val_type,
        boost::bidirectional_traversal_tag
        >
    {
    public:
        typedef Traits traits;
        typedef typename traits::val_type val_type;
        typedef typename traits::val_ptr_type val_ptr_type;
        typedef typename traits::node_ptr_type node_ptr_type;
        typedef typename traits::const_node_ptr_type const_node_ptr_type;

        /** Default constructor. */
        I_List_Iterator_Base() : _node_ptr(0) {}
        /** Copy constructor. */
        template <bool other_dir_tag>
        I_List_Iterator_Base(const I_List_Iterator_Base< Traits, other_dir_tag >& rhs) : _node_ptr(rhs._node_ptr) {}
        /** Constructor from node pointer; explicit. */
        explicit I_List_Iterator_Base(const_node_ptr_type node_ptr) : _node_ptr(node_ptr) {}

        /** Retrieve stored pointer without dereferencing it. */
        val_ptr_type get() const { return traits::to_val_ptr(_node_ptr); }

    private:
        friend class boost::iterator_core_access;
        friend class I_List_Iterator_Base< Traits, not dir_tag >;

        template <bool other_dir_tag>
        bool equal(const I_List_Iterator_Base< Traits, other_dir_tag>& rhs) const { return _node_ptr == rhs._node_ptr; }

        val_type& dereference() const { return *get(); }

        void advance(bool dir) { _node_ptr = (dir? traits::get_next(_node_ptr) : traits::get_prev(_node_ptr)); }
        void increment() { advance(dir_tag); }
        void decrement() { advance(not dir_tag); }

        node_ptr_type _node_ptr;
    };

    /** Const iterator. */
    template <class Traits, bool dir_tag>
    class I_List_Const_Iterator : public I_List_Iterator_Base< Traits, dir_tag  >
    {
    public:
        typedef typename Traits::val_type val_type;
        typedef typename Traits::node_ptr_type node_ptr_type;
        typedef typename Traits::const_node_ptr_type const_node_ptr_type;

        /** Default constructor. */
        I_List_Const_Iterator() {}
        /** Copy constructor. */
        template <bool other_dir_tag>
        I_List_Const_Iterator(const I_List_Const_Iterator< Traits, other_dir_tag >& rhs) : I_List_Iterator_Base< Traits, dir_tag >(rhs) {}
        /** Constructor from node pointer; explicit. */
        explicit I_List_Const_Iterator(const_node_ptr_type node_ptr) : I_List_Iterator_Base< Traits, dir_tag >(node_ptr) {}

        const val_type& operator * () { return ((I_List_Iterator_Base< Traits, dir_tag >*)this)->operator *(); }
    };

    template <class Traits, bool dir_tag>
    class I_List_Iterator : public I_List_Iterator_Base< Traits, dir_tag  >
    {
    public:
        typedef typename Traits::val_type val_type;
        typedef typename Traits::node_ptr_type node_ptr_type;
        typedef typename Traits::const_node_ptr_type const_node_ptr_type;

        /** Default constructor. */
        I_List_Iterator() {}
        /** Copy constructor. */
        template <bool other_dir_tag>
        I_List_Iterator(const I_List_Iterator< Traits, other_dir_tag >& rhs) : I_List_Iterator_Base< Traits, dir_tag >(rhs) {}
        /** Constructor from node pointer; explicit. */
        explicit I_List_Iterator(const node_ptr_type node_ptr) : I_List_Iterator_Base< Traits, dir_tag >(node_ptr) {}

        /** Constructor from const iterator; explicit. */
        explicit I_List_Iterator(const I_List_Const_Iterator< Traits, dir_tag >& rhs) : I_List_Iterator_Base< Traits, dir_tag >(rhs) {}
        /** Conversion to const iterator. */
        operator const I_List_Const_Iterator< Traits, dir_tag >& () const { return *(I_List_Const_Iterator< Traits, dir_tag >*)this; }
        operator I_List_Const_Iterator< Traits, dir_tag >& () { return *(I_List_Const_Iterator< Traits, dir_tag >*)this; }

        val_type& operator * () { return ((I_List_Iterator_Base< Traits, dir_tag >*)this)->operator *(); }
    };

}
using detail::I_List_Iterator;
using detail::I_List_Const_Iterator;


/** Intrusive list. */
template <class Traits, class Size_Type = size_t>
class I_List
{
public:
    typedef Traits traits;
    typedef Size_Type size_type;

    typedef typename traits::val_type val_type;
    typedef typename traits::val_ptr_type val_ptr_type;
    typedef typename traits::const_val_ptr_type const_val_ptr_type;
    typedef typename traits::node_type node_type;
    typedef typename traits::node_ptr_type node_ptr_type;
    typedef typename traits::const_node_ptr_type const_node_ptr_type;

    typedef I_List_Iterator< traits, true > iterator;
    typedef I_List_Iterator< traits, false > reverse_iterator;
    typedef I_List_Const_Iterator< traits, true > const_iterator;
    typedef I_List_Const_Iterator< traits, false > const_reverse_iterator;

    /** Constructor.
     * @param head_node_ptr Pointer to head node, required.
     */
    I_List(node_ptr_type head_node_ptr) : _size(0), _head_node_ptr(head_node_ptr)
    {
        traits::set_next(head_node_ptr, head_node_ptr);
        traits::set_prev(head_node_ptr, head_node_ptr);
    }

    size_type size() const { return _size; }
    bool empty() const { return size() == 0; }

    const_iterator begin() const { return const_iterator(traits::get_next(_head_node_ptr)); }
    const_iterator end() const { return const_iterator(_head_node_ptr); }
    const_reverse_iterator rbegin() const { return const_reverse_iterator(traits::get_prev(_head_node_ptr)); }
    const_reverse_iterator rend() const { return const_reverse_iterator(_head_node_ptr); }
    iterator begin() { return iterator(static_cast<const I_List*>(this)->begin()); }
    iterator end() { return iterator(static_cast<const I_List*>(this)->end()); }
    reverse_iterator rbegin() { return reverse_iterator(static_cast<const I_List*>(this)->rbegin()); }
    reverse_iterator rend() { return reverse_iterator(static_cast<const I_List*>(this)->rend()); }

    const val_type& front() const { return *(traits::to_val_ptr(traits::get_next(_head_node_ptr))); }
    const val_type& back() const { return *(traits::to_val_ptr(traits::get_prev(_head_node_ptr))); }
    val_type& front() { return const_cast<val_type&>(static_cast<const I_List*>(this)->front()); }
    val_type& back() { return const_cast<val_type&>(static_cast<const I_List*>(this)->back()); }

    void push_back(val_ptr_type val_ptr)
    { insert_after(traits::get_prev(_head_node_ptr), val_ptr); }
    void push_front(val_ptr_type val_ptr)
    { insert_after(_head_node_ptr, val_ptr); }
    node_ptr_type insert(const_iterator it, val_ptr_type val_ptr)
    { insert_after(it.get(), val_ptr); }

    iterator erase(const_iterator it) { return iterator(erase(it.get())); }

    bool check() const
    {
        size_t n_elem = 0;
        node_ptr_type crt = _head_node_ptr;
        do
        {
            ASSERT(traits::get_prev(traits::get_next(crt)) == crt);
            crt = traits::get_next(crt);
            ++n_elem;
        } while (crt != _head_node_ptr);
        ASSERT(_size == n_elem - 1);
        return true;
    }

private:
    I_List& operator = (const I_List&) { return *this; }

    node_ptr_type insert_after(node_ptr_type node_ptr, val_ptr_type val_ptr)
    {
        ASSERT(_size < std::numeric_limits< Size_Type >::max());
        node_ptr_type next_node_ptr = traits::get_next(node_ptr);
        node_ptr_type new_node_ptr = traits::to_node_ptr(val_ptr);
        traits::set_prev(next_node_ptr, new_node_ptr);
        traits::set_next(new_node_ptr, next_node_ptr);
        traits::set_prev(new_node_ptr, node_ptr);
        traits::set_next(node_ptr, new_node_ptr);
        ++_size;
        return new_node_ptr;
    }

    node_ptr_type erase(node_ptr_type node_ptr)
    {
        ASSERT(_size > 0);
        traits::set_prev(traits::get_next(node_ptr), traits::get_prev(node_ptr));
        traits::set_next(traits::get_prev(node_ptr), traits::get_next(node_ptr));
        --_size;
        return traits::get_next(node_ptr);
    }

    size_type _size;
    node_ptr_type _head_node_ptr;
};


#endif
