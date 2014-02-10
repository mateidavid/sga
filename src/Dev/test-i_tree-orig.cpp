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

    B() { clog << "constructing B at: " << (void*) this << '\n'; }
    B(const B& rhs)
    : _r_start(rhs._r_start), _c_start(rhs._c_start),
    _parent(rhs._parent), _l_child(rhs._l_child), _r_child(rhs._r_child), _col(rhs._col)
    {}
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

typedef factory::Factory< B > fact_type;
typedef fact_type::ptr_type ptr_type;
typedef fact_type::const_ptr_type const_ptr_type;
typedef fact_type::ref_type ref_type;
typedef fact_type::const_ref_type const_ref_type;

typedef boost::intrusive::multiset< B, boost::intrusive::value_traits< Value_Traits < B > > > itree_type;

int main()
{
    clog << "--- type sizes:\n";
    clog << "sizeof(val_type)=" << sizeof(fact_type::val_type) << '\n';
    clog << "sizeof(ptr_type)=" << sizeof(ptr_type) << '\n';
    clog << "sizeof(wrapper_type)=" << sizeof(factory::detail::Factory_Wrapper< fact_type::val_type >) << '\n';

    clog << "--- type names:"
         << "\nB: " << typeid(B).name()
         << "\nconst B: " << typeid(const B).name()
         << "\nFactory< B >: " << typeid(factory::Factory< B >).name()
         << "\nBounded_Pointer< B >: " << typeid(factory::Bounded_Pointer< B >).name()
         << "\nBounded_Pointer< const B >: " << typeid(factory::Bounded_Pointer< const B >).name()
         << "\nBounded_Reference< B >: " << typeid(factory::Bounded_Reference< B >).name()
         << "\nBounded_Reference< const B >: " << typeid(factory::Bounded_Reference< const B >).name()
         << "\nHolder< B >: " << typeid(factory::Holder< B >).name()
         << "\nValue_Traits< B >: " << typeid(Value_Traits< B >).name()
         << "\nNode_Traits< B >: " << typeid(Node_Traits< B >).name()
         << "\nvoid: " << typeid(void).name()
         << "\nconst void: " << typeid(const void).name()
         << '\n';

    clog << "--- constructing factory\n";
    fact_type f(true);

    clog << "--- constructing itree\n";
#ifndef NONSTATIC_VALUE_TRAITS
    clog << "using static value traits\n";
    itree_type l;
#else
    clog << "using non-static value traits\n";
    itree_type l(std::less< B >(), Value_Traits< B >(42));
#endif

    size_t n = 10;

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

    clog << "--- inserting elements in tree\n";
    for (size_t i = 0; i < n; ++i)
    {
        clog << "--- inserting element at index " << i << '\n';
        l.insert(*ptr_a[i]);
    }

    clog << "--- factory:\n" << f;
    clog << "--- tree:\n";
    auto it = l.begin();
    auto it_end = l.end();
    while (it != it_end)
    {
        clog << "--- next element\n";
        clog << *it << '\n';
        ++it;
    }

    clog << "--- removing even index elements from tree\n";
    for (size_t i = 0; i < n; ++i)
    {
        if (i % 2 > 0)
        {
            continue;
        }
        itree_type::const_iterator it = l.iterator_to(*ptr_a[i]);
        l.erase(it);
    }

    clog << "--- factory:\n" << f;
    clog << "--- tree:\n";
    for (it = l.begin(); it != l.end(); ++it)
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
