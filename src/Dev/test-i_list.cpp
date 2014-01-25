#include <iostream>
#include <deque>
#include <cassert>
#include "factory.hpp"
#include "i_list.hpp"

using namespace std;


struct A
{
    typedef factory::Identifier< A > idn_type;
    typedef factory::Const_Identifier< A > const_idn_type;

    size_t _r_start;
    size_t _c_start;

    idn_type _prev;
    idn_type _next;

    A() { cerr << "constructing A at: " << (void*) this << '\n'; }

    A* operator -> () { return this; }
    const A* operator -> () const { return this; }
};

ostream& operator <<(ostream& os, const A& rhs)
{
    os << "[_r_start=" << rhs._r_start << ",_c_start=" << rhs._c_start
       << ",_prev=" << rhs._prev << ",_next=" << rhs._next << "]";
    return os;
}

struct traits
{
    typedef A val_type;
    typedef A::idn_type val_ptr_type;
    typedef A& val_ref_type;

    typedef A node_type;
    typedef A::idn_type node_ptr_type;
    typedef A::const_idn_type const_node_ptr_type;

    static node_ptr_type get_next(const node_ptr_type n) { return n->_next; }
    static void set_next(node_ptr_type n, node_ptr_type next) { n->_next = next; }
    static node_ptr_type get_prev(const node_ptr_type n) { return n->_prev; }
    static void set_prev(node_ptr_type n, node_ptr_type prev) { n->_prev = prev; }

    static node_ptr_type to_node_ptr (val_ptr_type val_ptr) { return node_ptr_type(val_ptr); }
    static val_ptr_type to_val_ptr(node_ptr_type n) { return val_ptr_type(n); }
};


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

    typedef factory::Factory< A > factory_type;
    typedef factory_type::idn_type idn_type;
    factory_type f;
    f.set_active();
    
    idn_type a;
    a = f.add_elem();
    I_List<traits> l(a);

    a = f.add_elem();
    a->_r_start = 5;
    a->_c_start = 17;
    l.push_back(a);

    a = f.add_elem();
    a->_r_start = 15;
    a->_c_start = 23;
    l.push_back(a);

    a = f.add_elem();
    a->_r_start = 8;
    a->_c_start = 1;
    l.push_front(a);

    cout << f;

    cout << "l:\n";
    for (auto it = l.begin(); it != l.end(); ++it)
    {
        clog << *it << '\n';
    }

    return 0;
}
