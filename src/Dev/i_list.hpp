#ifndef __I_LIST_HPP
#define __I_LIST_HPP

#include <boost/iterator/iterator_facade.hpp>


template <class Traits>
class I_List_Iterator : public boost::iterator_facade<
    I_List_Iterator< Traits >,
    typename Traits::val_type,
    boost::bidirectional_traversal_tag
    >
{
public:
    typedef Traits traits;
    typedef typename traits::val_type val_type;
    typedef typename traits::val_ref_type val_ref_type;
    typedef typename traits::node_ptr_type node_ptr_type;

    I_List_Iterator() : _node_ptr(0) {}
    explicit I_List_Iterator(node_ptr_type node_ptr) : _node_ptr(node_ptr) {}

private:
    friend class boost::iterator_core_access;

    val_ref_type dereference() const
    { return *(traits::to_val_ptr(_node_ptr)); }

    bool equal(const I_List_Iterator& rhs) const
    { return _node_ptr == rhs._node_ptr; }

    void increment() { _node_ptr = traits::get_next(_node_ptr); }
    void decrement() { _node_ptr = traits::get_prev(_node_ptr); }

    node_ptr_type _node_ptr;
};


template <class Traits, class Size_Type = size_t>
class I_List
{
public:
    typedef Traits traits;
    typedef Size_Type size_type;

    typedef typename traits::val_type val_type;
    typedef typename traits::val_ptr_type val_ptr_type;
    typedef typename traits::val_ref_type val_ref_type;

    typedef typename traits::node_type node_type;
    typedef typename traits::node_ptr_type node_ptr_type;

    typedef I_List_Iterator< traits > iterator;

    I_List(node_ptr_type head_node_ptr) : _size(0), _head_node_ptr(head_node_ptr)
    {
        traits::set_next(head_node_ptr, head_node_ptr);
        traits::set_prev(head_node_ptr, head_node_ptr);
    }

    size_t size() const { return _size; }

    iterator begin() { return iterator(traits::get_next(_head_node_ptr)); }
    iterator end() { return iterator(_head_node_ptr); }

    val_ref_type front()
    { return *(traits::to_val_ptr(traits::get_next(_head_node_ptr))); }
    val_ref_type back()
    { return *(traits::to_val_ptr(traits::get_prev(_head_node_ptr))); }

    void push_back(val_ptr_type val_ptr)
    { insert_after(traits::get_prev(_head_node_ptr), val_ptr); }
    void push_front(val_ptr_type val_ptr)
    { insert_after(_head_node_ptr, val_ptr); }

private:
    I_List& operator = (const I_List&) { return *this; }

    node_ptr_type insert_after(node_ptr_type node_ptr, val_ptr_type val_ptr)
    {
        node_ptr_type next_node_ptr = traits::get_next(node_ptr);
        node_ptr_type new_node_ptr = traits::to_node_ptr(val_ptr);
        traits::set_prev(next_node_ptr, new_node_ptr);
        traits::set_next(new_node_ptr, next_node_ptr);
        traits::set_prev(new_node_ptr, node_ptr);
        traits::set_next(node_ptr, new_node_ptr);
        return new_node_ptr;
    }

    size_type _size;
    node_ptr_type _head_node_ptr;
};


#endif
