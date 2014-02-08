#include <iostream>
#include <boost/intrusive/set.hpp>

#include "common.hpp"
#include "factory.hpp"


using namespace std;


struct B
{
    typedef factory::Factory< B > fact_type;
    typedef fact_type::ptr_type ptr_type;

    size_t _r_start;
    size_t _c_start;

    ptr_type _parent;
    ptr_type _l_child;
    ptr_type _r_child;
    int _col;

    B() { cerr << "constructing B at: " << (void*) this << '\n'; }
    B(const B& rhs)
    : _r_start(rhs._r_start), _c_start(rhs._c_start),
    _parent(rhs._parent), _l_child(rhs._l_child), _r_child(rhs._r_child), _col(rhs._col)
    {}

    B* operator -> () { return this; }
    const B* operator -> () const { return this; }

    B& operator ++ () { ++_r_start; return *this; }
};

bool operator < (const B& lhs, const B& rhs) { return lhs._r_start < rhs._r_start; }

ostream& operator <<(ostream& os, const B& rhs)
{
    os << "[_r_start=" << rhs._r_start << ",_c_start=" << rhs._c_start
       << ",_parent=" << rhs._parent
       << ",_l_child=" << rhs._l_child
       << ",_r_child=" << rhs._r_child
       << ",_col=" << rhs._col
       << "]";
    return os;
}

template <class T>
struct Node_Traits
{
    typedef factory::Holder< T > node;
    typedef factory::Factory< T > fact_type;
    typedef typename fact_type::ptr_type node_ptr;
    typedef typename fact_type::const_ptr_type const_node_ptr;
    typedef int color;

    static node_ptr get_parent(const_node_ptr n) { return n->_parent; }
    static void set_parent(node_ptr n, node_ptr ptr) { n->_parent = ptr; }
    static node_ptr get_left(const_node_ptr n) { return n->_l_child; }
    static void set_left(node_ptr n, node_ptr ptr) { n->_l_child = ptr; }
    static node_ptr get_right(const_node_ptr n) { return n->_r_child; }
    static void set_right(node_ptr n, node_ptr ptr) { n->_r_child = ptr; }
    static color get_color(const_node_ptr n) { return n->_col; }
    static void set_color(node_ptr n, color c) { n->_col = c ; }
    static color black() { return 0; }
    static color red() { return 1; }
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

    static node_ptr to_node_ptr (reference value) {  return &value; }
    static const_node_ptr to_node_ptr (const_reference value) { return &value; }
    static pointer to_value_ptr(node_ptr n) { return n; }
    static const_pointer to_value_ptr(const_node_ptr n) { return n; }

    Value_Traits* operator -> () { return this; }
    const Value_Traits* operator -> () const { return this; }
};

typedef factory::Factory< B > fact_type;
typedef fact_type::ptr_type ptr_type;
typedef fact_type::const_ptr_type const_ptr_type;
typedef fact_type::ref_type ref_type;
typedef fact_type::const_ref_type const_ref_type;

typedef boost::intrusive::multiset< B, boost::intrusive::value_traits< Value_Traits < B > > > itree_type;

int main()
{
    fact_type f(true);
    itree_type t;

    factory::detail::Bounded_Pointer<const Value_Traits< B >, unsigned int> wtf1;
    const factory::detail::Bounded_Pointer<const void, unsigned int>& wtf2(wtf1);

    ptr_type a;

    a = f.new_elem();
    a->_r_start = 5;
    a->_c_start = 17;
    t.insert(*a);

    a = f.new_elem();
    a->_r_start = 15;
    a->_c_start = 23;
    t.insert(*a);

    a = f.new_elem();
    a->_r_start = 8;
    a->_c_start = 1;
    t.insert(*a);

    cout << "factory:\n" << f;

    cout << "tree:\n";
    for (auto it = t.begin(); it != t.end(); ++it)
    {
        clog << *it << '\n';
    }

    return 0;
}
