#include <iostream>
#include <deque>
#include <cassert>
#include "factory.hpp"

using namespace std;

typedef factory::Factory< size_t > fact_type;
typedef fact_type::ptr_type ptr_type;
typedef fact_type::const_ptr_type const_ptr_type;
typedef fact_type::ref_type ref_type;
typedef fact_type::const_ref_type const_ref_type;

void fcn_const_ref(const_ref_type r)
{
    cout << "fcn: r=" << r << '\n';
}

int main()
{
    fact_type f;
    f.set_active();
    ptr_type p1, p2, p3;

    clog << "sizeof(val_type)=" << sizeof(fact_type::val_type) << '\n';
    clog << "sizeof(ptr_type)=" << sizeof(ptr_type) << '\n';
    clog << "sizeof(wrapper_type)=" << sizeof(factory::detail::Factory_Wrapper<fact_type::val_type>) << '\n';

    p1 = f.new_elem();
    *p1 = 17;
    p2 = f.new_elem();
    *p2 = 9;
    p3 = f.new_elem();
    *p3 = 24;
    cout << "factory after new 17,9,24:\n" << f;

    const_ptr_type p4(p3);
    cout << "(p4 == p3): " << (p4 == p3) << '\n';
    cout << "(p4 != p3): " << (p4 != p3) << '\n';
    cout << "(p4): " << (bool)p4 << '\n';
    const_ptr_type p5(p4);
    f.del_elem(p2);
    p2 = f.new_elem();
    *p2 = 15;
    p3 = f.new_elem();
    *p3 = 42;
    cout << "factory after del 9, new 15,42:\n" << f;

    //fact_type::const_idn_type p4(p3);
    cout << "p4=" << p4 << '\n';
    //fact_type::idn_type p5(p4);
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

    ref_type r1(*p1);
    cout << "p1=" << p1 << '\n';
    cout << "r1(*p1)=" << r1 << '\n';

    size_t& q1 = r1;
    size_t q2 = r1;
    cout << "r1=" << r1 << '\n';
    cout << "q1=" << q1 << '\n';
    cout << "q2=" << q2 << '\n';

    ref_type r2 = *p2;
    ref_type r3 = r2;
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

    const_ref_type cr1(r1);
    //cr1 = 18;
    cout << "cr1=" << cr1 << '\n';
    r2 = cr1;
    cout << "r2=" << r2 << '\n';
    ref_type r4(cr1);
    r4 = 33;
    cout << "r4=" << r4 << '\n';

    //*p4 = 28;
    cout << "p4=" << p4 << '\n';

    factory::Bounded_Pointer< const void > vp1(p1);
    ptr_type p6(vp1);

    const ptr_type cp7(p1);
    ptr_type p8(cp7);

    /*
    factory::Holder< size_t > h1;
    cout << "sizeof(holder)=" << sizeof(h1) << "\n";
    h1 = 59;
    factory::Bounded_Reference< size_t > r4 = h1;
    cout << "r4=" << r4 << "\n";
    cout << f;
    */

    /*
    typedef Factory< A > fact_type;
    typedef fact_type::ptr_type ptr_type;
    fact_type f;
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
