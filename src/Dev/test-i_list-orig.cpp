#include <iostream>
#include <boost/intrusive/list.hpp>
#include <boost/intrusive/pointer_traits.hpp>
#include <typeinfo>

#include "common.hpp"
#include "factory.hpp"


using namespace std;


struct A
{
    typedef factory::Factory< A > fact_type;
    typedef fact_type::ptr_type ptr_type;

    size_t _r_start;
    size_t _c_start;

    ptr_type _prev;
    ptr_type _next;

    A() { clog << "constructing A at: " << (void*) this << '\n'; }
};

ostream& operator <<(ostream& os, const A& rhs)
{
    os << "[_r_start=" << rhs._r_start << ",_c_start=" << rhs._c_start
       << ",_prev=" << rhs._prev << ",_next=" << rhs._next << "]";
    return os;
}

template <class T>
struct Node_Traits
{
    typedef factory::Holder< T > node;
    typedef factory::Factory< T > fact_type;
    typedef typename fact_type::ptr_type node_ptr;
    typedef typename fact_type::const_ptr_type const_node_ptr;

    static node_ptr get_next(const_node_ptr n) { return n->_next; }
    static void set_next(node_ptr n, node_ptr next) { n->_next = next; }
    static node_ptr get_previous(const_node_ptr n) { return n->_prev; }
    static void set_previous(node_ptr n, node_ptr prev) { n->_prev = prev; }
};

template <class T>
struct Value_Traits
{
    typedef Node_Traits< T > node_traits;
    typedef typename node_traits::node_ptr node_ptr;
    typedef typename node_traits::const_node_ptr const_node_ptr;
    typedef T value_type;
    typedef node_ptr pointer;
    typedef const_node_ptr const_pointer;
    typedef typename node_traits::fact_type::ref_type reference;
    typedef typename node_traits::fact_type::const_ref_type const_reference;

    static const boost::intrusive::link_mode_type link_mode = boost::intrusive::normal_link;

#ifndef NONSTATIC_VALUE_TRAITS
    static node_ptr to_node_ptr (reference value) {  return &value; }
    static const_node_ptr to_node_ptr (const_reference value) { return &value; }
    static pointer to_value_ptr(node_ptr n) { return n; }
    static const_pointer to_value_ptr(const_node_ptr n) { return n; }
#else
    node_ptr to_node_ptr (reference value) const { ASSERT(_val == 42); return &value; }
    const_node_ptr to_node_ptr (const_reference value) const { ASSERT(_val == 42); return &value; }
    pointer to_value_ptr(node_ptr n) const { ASSERT(_val == 42); return n; }
    const_pointer to_value_ptr(const_node_ptr n) const { ASSERT(_val == 42); return n; }
    Value_Traits(int val) : _val(val) {}
    int _val;
#endif
};

typedef factory::Factory< A > fact_type;
typedef fact_type::ptr_type ptr_type;
typedef fact_type::const_ptr_type const_ptr_type;
typedef fact_type::ref_type ref_type;
typedef fact_type::const_ref_type const_ref_type;
typedef boost::intrusive::pointer_traits< ptr_type >::rebind_pointer< void >::type void_ptr_type;

typedef boost::intrusive::list< A, boost::intrusive::value_traits< Value_Traits < A > > > ilist_type;


int main()
{
    clog << "--- type sizes:\n";
    clog << "sizeof(val_type)=" << sizeof(fact_type::val_type) << '\n';
    clog << "sizeof(ptr_type)=" << sizeof(ptr_type) << '\n';
    clog << "sizeof(wrapper_type)=" << sizeof(factory::detail::Factory_Wrapper< fact_type::val_type >) << '\n';

    clog << "--- type names:"
         << "\nA: " << typeid(A).name()
         << "\nconst A: " << typeid(const A).name()
         << "\nFactory< A >: " << typeid(factory::Factory< A >).name()
         << "\nBounded_Pointer< A >: " << typeid(factory::Bounded_Pointer< A >).name()
         << "\nBounded_Pointer< const A >: " << typeid(factory::Bounded_Pointer< const A >).name()
         << "\nBounded_Reference< A >: " << typeid(factory::Bounded_Reference< A >).name()
         << "\nBounded_Reference< const A >: " << typeid(factory::Bounded_Reference< const A >).name()
         << "\nHolder< A >: " << typeid(factory::Holder< A >).name()
         << "\nValue_Traits< A >: " << typeid(Value_Traits< A >).name()
         << "\nNode_Traits< A >: " << typeid(Node_Traits< A >).name()
         << "\npointer_traits<ptr_type>::rebind_pointer<void>::type: " << typeid(void_ptr_type).name()
         << "\nvoid: " << typeid(void).name()
         << "\nconst void: " << typeid(const void).name()
         << '\n';

    clog << "--- constructing factory\n";
    fact_type f(true);

    clog << "--- constructing a holder object\n";
    factory::Holder< A > h;
    A& a_ref = h;
    clog << "&h._r_start = " << &(((A&)h)._r_start) << "\n";
    clog << "&a_ref = " << &a_ref << "\n";
    clog << "&a_ref._r_start = " << &a_ref._r_start<< "\n";
    const void* vp = &h;
    clog << "vp = " << vp << "\n";
    ASSERT(vp == &a_ref);

    clog << "--- constructing ilist\n";
#ifndef NONSTATIC_VALUE_TRAITS
    clog << "using static value traits\n";
    ilist_type l;
#else
    clog << "using non-static value traits\n";
    ilist_type l(Value_Traits< A >(42));
#endif

    size_t n = 4;

    clog << "--- constructing vector of pointers\n";
    vector< ptr_type > ptr_a(n);

    clog << "--- allocating elements\n";
    for (size_t i = 0; i < n; ++i)
    {
        clog << "--- allocating element at index " << i << '\n';
        ptr_type& a = ptr_a[i];
        a = f.new_elem();
        a->_r_start = i;
        a->_c_start = 50 + i;
    }

    clog << "--- inserting elements in list\n";
    for (size_t i = 0; i < n; ++i)
    {
        clog << "--- inserting element at index " << i << '\n';
        if (i < n / 2)
        {
            l.push_back(*ptr_a[i]);
        }
        else
        {
            l.push_front(*ptr_a[i]);
        }
    }

    clog << "--- factory:\n" << f;

    clog << "--- list:\n";
    for (auto it = l.begin(); it != l.end(); ++it)
    {
        clog << "--- next element\n";
        clog << *it << '\n';
    }

    clog << "--- removing even index elements from list\n";
    for (size_t i = 0; i < n; ++i)
    {
        if (i % 2 > 0)
        {
            continue;
        }
        ilist_type::const_iterator it = l.iterator_to(*ptr_a[i]);
        l.erase(it);
    }

    clog << "--- factory:\n" << f;
    clog << "--- list:\n";
    for (auto it = l.begin(); it != l.end(); ++it)
    {
        clog << *it << '\n';
    }

    clog << "--- deallocating elements\n";
    for (size_t i = 0; i < n; ++i)
    {
        f.del_elem(ptr_a[i]);
    }

    clog << "--- exiting\n";

    return 0;
}
