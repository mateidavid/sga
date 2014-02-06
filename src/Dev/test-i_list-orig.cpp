#include <iostream>
#include <boost/intrusive/list.hpp>

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

    A() { cerr << "constructing A at: " << (void*) this << '\n'; }

    A* operator -> () { return this; }
    const A* operator -> () const { return this; }

    A& operator ++ () { ++_r_start; return *this; }
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
    typedef A value_type;
    typedef node_ptr pointer;
    typedef const_node_ptr const_pointer;
    typedef typename node_traits::fact_type::ref_type reference;
    typedef typename node_traits::fact_type::const_ref_type const_reference;

    static const boost::intrusive::link_mode_type link_mode = boost::intrusive::normal_link;

    static node_ptr to_node_ptr (reference value) {  return &value; }
    static const_node_ptr to_node_ptr (const_reference value) { return &value; }
    static pointer to_value_ptr(node_ptr n) { return n; }
    static const_pointer to_value_ptr(const_node_ptr n) { return n; }
};

typedef factory::Factory< A > fact_type;
typedef fact_type::ptr_type ptr_type;
typedef fact_type::const_ptr_type const_ptr_type;
typedef fact_type::ref_type ref_type;
typedef fact_type::const_ref_type const_ref_type;

typedef boost::intrusive::list< A, boost::intrusive::value_traits< Value_Traits < A > > > ilist_type;

int main()
{
    /*
    Factory_Wrapper<int> w;
    int& i = w._key;
    */

    /*
    typedef Factory< size_t > factory_type;
    factory_type f;
    factory_type::ptr_type::_factory_ptr = &f;
    factory_type::ptr_type p1, p2, p3;

    clog << "sizeof(val_type)=" << sizeof(factory_type::val_type) << '\n';
    clog << "sizeof(ptr_type)=" << sizeof(factory_type::ptr_type) << '\n';
    clog << "sizeof(wrapper_type)=" << sizeof(factory_type::wrapper_type) << '\n';

    p1 = f.new_elem();
    *p1 = 17;
    p2 = f.new_elem();
    *p2 = 9;
    p3 = f.new_elem();
    *p3 = 24;
    f.print();

    f.del_elem(p2);
    f.print();

    p2 = f.new_elem();
    *p2 = 15;
    p3 = f.new_elem();
    *p3 = 42;
    f.print();
    */

    fact_type f(true);
    ilist_type l;
    
    ptr_type a;

    a = f.new_elem();
    a->_r_start = 5;
    a->_c_start = 17;
    l.push_back(*a);

    a = f.new_elem();
    a->_r_start = 15;
    a->_c_start = 23;
    l.push_back(*a);

    a = f.new_elem();
    a->_r_start = 8;
    a->_c_start = 1;
    l.push_front(*a);

    cout << f;

    cout << "l: beg->end\n";
    for (auto it = l.begin(); it != l.end(); ++it)
    {
        clog << *it << '\n';
    }
    cout << "l: rbeg->rend\n";
    for (auto it = l.rbegin(); it != l.rend(); ++it)
    {
        clog << *it << '\n';
    }

    return 0;
}
