#include <iostream>
#include <deque>
#include <cassert>
#include "factory.hpp"

using namespace std;


void fcn_const_ref(factory::Bounded_Reference< size_t, true > r)
{
    cout << "fcn: r=" << r << '\n';
}

int main()
{
    /*
    Factory_Wrapper<int> w;
    int& i = w._key;
    */

    typedef factory::Factory< size_t > factory_type;
    factory_type f;
    f.set_active();
    factory_type::ptr_type p1, p2, p3;

    clog << "sizeof(val_type)=" << sizeof(factory_type::val_type) << '\n';
    clog << "sizeof(ptr_type)=" << sizeof(factory_type::ptr_type) << '\n';
    clog << "sizeof(wrapper_type)=" << sizeof(factory::detail::Factory_Wrapper<factory_type::val_type>) << '\n';

    p1 = f.new_elem();
    (*p1) = 17;
    p2 = f.new_elem();
    *p2 = 9;
    p3 = f.new_elem();
    *p3 = 24;
    cout << f;

    factory_type::const_ptr_type p4(p3);
    factory_type::const_ptr_type p5(p4);
    f.del_elem((factory_type::const_ptr_type)p2);
    cout << f;

    p2 = f.new_elem();
    *p2 = 15;
    p3 = f.new_elem();
    *p3 = 42;
    cout << f;

    //factory_type::const_idn_type p4(p3);
    cout << "p4=" << p4 << '\n';

    //factory_type::idn_type p5(p4);
    cout << "p5=" << p5 << '\n';

    p4 = p3;
    //p2 = p4;

    cout << "p1=" << p1 << '\n';
    cout << "p2=" << p2 << '\n';
    p1++;
    ++p2;
    cout << "p1=" << p1 << '\n';
    cout << "p2=" << p2 << '\n';

    cout << "-----\n";

    factory::Bounded_Reference< size_t, false > r1(*p1);
    cout << "p1=" << p1 << '\n';
    cout << "r1(*p1)=" << r1 << '\n';

    size_t& q1 = r1;
    size_t q2 = r1;
    cout << "r1=" << r1 << '\n';
    cout << "q1=" << q1 << '\n';
    cout << "q2=" << q2 << '\n';

    factory::Bounded_Reference< size_t, false > r2 = *p2;
    factory::Bounded_Reference< size_t, false > r3 = r2;
    cout << "r1=" << r1 << '\n';
    cout << "r2=" << r2 << '\n';
    cout << "r3=" << r3 << '\n';
    cout << "q1=" << q1 << '\n';

    r2 = r1;
    cout << "r1=" << r1 << '\n';
    cout << "r2=" << r2 << '\n';

    size_t& q3 = (size_t&)r3;
    q3++;
    cout << "r3=" << r3 << '\n';

    cout << f;
    fcn_const_ref(*p1);
    fcn_const_ref(r2);

    factory::Holder< size_t > h1;
    cout << "sizeof(holder)=" << sizeof(h1) << "\n";
    h1 = 59;
    factory::Bounded_Reference< size_t, false > r4 = h1;
    cout << "r4=" << r4 << "\n";
    cout << f;

    /*
    typedef Factory< A > factory_type;
    typedef factory_type::ptr_type ptr_type;
    factory_type f;
    ptr_type::_factory_ptr = &f;
    
    ptr_type a;
    a = f.new_elem();
    I_List<traits> l(a);

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

    f.print();

    clog << "l:\n";
    for (auto it = l.begin(); it != l.end(); ++it)
    {
        clog << *it << '\n';
    }
    */

    return 0;
}
