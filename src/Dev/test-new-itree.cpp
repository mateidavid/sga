#include <iostream>

#include "common.hpp"
#include "factory.hpp"
#include "itree.hpp"


using namespace std;


struct B
{
    typedef factory::Factory< B > fact_type;
    typedef fact_type::ptr_type ptr_type;

    size_t _start;
    size_t _end;

    ptr_type _parent;
    ptr_type _l_child;
    ptr_type _r_child;
    int _col;
    size_t _max_end;
};

ostream& operator <<(ostream& os, const B& rhs)
{
    os << "[_start=" << rhs._start << ",_end=" << rhs._end
       << ",_parent=" << rhs._parent
       << ",_l_child=" << rhs._l_child
       << ",_r_child=" << rhs._r_child
       << ",_col=" << rhs._col
       << ",_max_end=" << rhs._max_end
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
    typedef size_t key_type;

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
    static key_type get_max_end(const_node_ptr n) { return n->_max_end; }
    static void set_max_end(node_ptr n, key_type k) { n->_max_end = k ; }
};

template <class T>
struct Value_Traits
{
    typedef Node_Traits< T > node_traits;
    typedef typename node_traits::node_ptr node_ptr;
    typedef typename node_traits::const_node_ptr const_node_ptr;
    typedef T value_type;
    typedef typename node_traits::key_type key_type;
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
    static key_type get_start(const_pointer n) { return n->_start; }
    static key_type get_start(const value_type* n) { return n->_start; }
    static key_type get_end(const_pointer n) { return n->_end; }
    static key_type get_end(const value_type* n) { return n->_end; }
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

typedef itree< Value_Traits< B > > iitree_type;
typedef boost::intrusive::bstree_algorithms< detail::ITree_Node_Traits< Value_Traits< B > > > bstree_algorithms;

ptr_type get_root(iitree_type& t)
{
    return bstree_algorithms::get_header(&*t.begin())->_parent;
}

void print_sub_tree(ptr_type r, size_t depth)
{
    if (!r)
    {
        return;
    }
    print_sub_tree(r->_l_child, depth + 1);
    for (size_t i = 0; i < depth; ++i)
    {
        clog << "  ";
    }
    clog << '[' << r->_start << ',' << r->_end << "] " << r->_max_end << '\n';
    print_sub_tree(r->_r_child, depth + 1);
}

void print_tree(iitree_type& t)
{
    print_sub_tree(get_root(t), 0);
}

bool check_max_ends(ptr_type node_ptr, size_t& max_end)
{
    if (!node_ptr)
    {
        max_end = 0;
        return true;
    }
    size_t max_end_left;
    size_t max_end_right;
    if (not check_max_ends(node_ptr->_l_child, max_end_left) or not check_max_ends(node_ptr->_r_child, max_end_right))
    {
        return false;
    }
    if (node_ptr->_max_end != max(node_ptr->_end, max(max_end_left, max_end_right)))
    {
        clog << "_max_end error: " << *node_ptr << '\n';
        return false;
    }
    max_end = node_ptr->_max_end;
    return true;
}

void check_max_ends(iitree_type& t, fact_type& f)
{
    ptr_type root_node = get_root(t);
    size_t max_end;
    if (not check_max_ends(root_node, max_end))
    {
        clog << "factory:\n" << f;
        clog << "tree:\n";
        for (auto const& e :t)
        {
            clog << e << '\n';
        }
        print_tree(t);
        exit(1);
    }
}


int main()
{
    clog << "----- type sizes:\n";
    clog << "sizeof(val_type)=" << sizeof(fact_type::val_type) << '\n';
    clog << "sizeof(ptr_type)=" << sizeof(ptr_type) << '\n';
    clog << "sizeof(wrapper_type)=" << sizeof(factory::detail::Factory_Wrapper< fact_type::val_type >) << '\n';

    clog << "----- type names:"
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
         << "\nconst_ptr_type: " << typeid(const_ptr_type).name()
         << "\nuncast_types< const_ptr_type >::element_type: "
         << typeid(boost::intrusive::detail::uncast_types<const_ptr_type >::element_type).name()
         << "\nuncast_types< const_ptr_type >::non_const_type: "
         << typeid(boost::intrusive::detail::uncast_types<const_ptr_type >::non_const_type).name()
         << "\nuncast_types< const_ptr_type >::non_const_pointer: "
         << typeid(boost::intrusive::detail::uncast_types<const_ptr_type >::non_const_pointer).name()
         << "\nconst_ptr_type::rebind<ptr_type>::type: " << typeid(const_ptr_type::rebind<ptr_type>::type).name()
         << "\nptr_type::rebind<const_ptr_type>::type: " << typeid(ptr_type::rebind<const_ptr_type>::type).name()
         << "\npointer_traits<const_ptr_type>::rebind_pointer<ptr_type>::type: "
         << typeid(boost::intrusive::pointer_traits<const_ptr_type>::rebind_pointer<ptr_type>::type).name()
         << "\npointer_traits<ptr_type>::rebind_pointer<const_ptr_type>::type: "
         << typeid(boost::intrusive::pointer_traits<ptr_type>::rebind_pointer<const_ptr_type>::type).name()
         << '\n';

    clog << "----- constructing factory\n";
    fact_type f(true);

    clog << "----- constructing iitree\n";
#ifndef NONSTATIC_VALUE_TRAITS
    clog << "using static value traits\n";
    iitree_type t;
#else
    clog << "using non-static value traits\n";
    iitree_type l(std::less< B >(), Value_Traits< B >(42));
#endif

    size_t n = 20;

    clog << "----- constructing vector of pointers\n";
    vector< ptr_type > ptr_a(n);

    clog << "----- allocating elements\n";
    for (size_t i = 0; i < n; ++i)
    {
        clog << "--- allocating element at index " << i << '\n';
        ptr_type& a = ptr_a[i];
        a = f.new_elem();
        a->_start = i;
        a->_end = i + 1 + (i % 3);
    }

    clog << "----- inserting elements in tree\n";
    for (size_t i = 0; i < n; ++i)
    {
        clog << "--- inserting element: " << *ptr_a[i] << '\n';
        B& b_ref = *ptr_a[i];
        (void)b_ref;
        t.insert(*ptr_a[i]);
        check_max_ends(t, f);
    }

    clog << "----- factory:\n" << f;
    clog << "----- tree:\n";
    auto it = t.begin();
    auto it_end = t.end();
    while (it != it_end)
    {
        clog << "--- next element\n";
        clog << *it << '\n';
        ++it;
    }

    clog << "----- removing even index elements from tree\n";
    for (size_t i = 0; i < n; ++i)
    {
        if (i % 2 > 0)
        {
            continue;
        }
        auto it = t.iterator_to(*ptr_a[i]);
        t.erase(it);
        check_max_ends(t, f);
    }

    clog << "----- factory:\n" << f;
    clog << "----- tree:\n";
    for (it = t.begin(); it != t.end(); ++it)
    {
        clog << *it << '\n';
    }
    print_tree(t);

    clog << "----- deallocating elements\n";
    for (size_t i = 0; i < n; ++i)
    {
        f.del_elem(ptr_a[i]);
    }

    clog << "----- exiting\n";


    return 0;
}
