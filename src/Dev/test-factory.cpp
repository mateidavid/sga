#include <iostream>
#include <deque>
#include <cassert>
#include "factory.hpp"

using namespace std;


int main()
{
    /*
    Factory_Wrapper<int> w;
    int& i = w._key;
    */

    typedef factory::Factory< size_t > factory_type;
    factory_type f;
    f.set_active();
    factory_type::idn_type p1, p2, p3;

    clog << "sizeof(val_type)=" << sizeof(factory_type::val_type) << '\n';
    clog << "sizeof(idn_type)=" << sizeof(factory_type::idn_type) << '\n';
    clog << "sizeof(wrapper_type)=" << sizeof(factory::detail::Factory_Wrapper<factory_type::val_type>) << '\n';

    p1 = f.new_elem();
    *p1 = 17;
    p2 = f.new_elem();
    *p2 = 9;
    p3 = f.new_elem();
    *p3 = 24;
    cout << f;

    f.del_elem(p2);
    cout << f;

    p2 = f.new_elem();
    *p2 = 15;
    p3 = f.new_elem();
    *p3 = 42;
    cout << f;

    factory_type::const_idn_type p4(p3);
    cout << "p4=" << p4 << '\n';

    //factory_type::idn_type p5(p4);
    //cout << "p5=" << p5 << '\n';

    //p4 = p3;
    //p2 = p4;

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
