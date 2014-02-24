#include <iostream>

#include <boost/intrusive/list.hpp>
#include "common.hpp"
#include "factory.hpp"
#include "itree.hpp"


using namespace std;


struct B
{
    typedef Factory< B > fact_type;
    typedef fact_type::ptr_type ptr_type;

    size_t _start;
    size_t _end;

    ptr_type _parent;
    ptr_type _l_child;
    ptr_type _r_child;
    int _col;
    size_t _max_end;

    ptr_type _list_next;
    ptr_type _list_prev;
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

bool intersect(const B& lhs, const B& rhs)
{
    return (lhs._start <= rhs._start and rhs._start <= lhs._end)
           or (rhs._start <= lhs._start and lhs._start <= rhs._end);
}

template <class T>
struct Node_Traits
{
    typedef Holder< T > node;
    typedef Factory< T > fact_type;
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

    static node_ptr to_node_ptr (reference value) { return &value; }
    static const_node_ptr to_node_ptr (const_reference value) { return &value; }
    static pointer to_value_ptr(node_ptr n) { return n; }
    static const_pointer to_value_ptr(const_node_ptr n) { return n; }
    static key_type get_start(const_pointer n) { return n->_start; }
    static key_type get_start(const value_type* n) { return n->_start; }
    static key_type get_end(const_pointer n) { return n->_end; }
    static key_type get_end(const value_type* n) { return n->_end; }
};


template <class T>
struct List_Node_Traits
{
    typedef Holder< T > node;
    typedef Factory< T > fact_type;
    typedef typename fact_type::ptr_type node_ptr;
    typedef typename fact_type::const_ptr_type const_node_ptr;

    static node_ptr get_next(const_node_ptr n) { return n->_list_next; }
    static void set_next(node_ptr n, node_ptr next) { n->_list_next = next; }
    static node_ptr get_previous(const_node_ptr n) { return n->_list_prev; }
    static void set_previous(node_ptr n, node_ptr prev) { n->_list_prev = prev; }
};

template <class T>
struct List_Value_Traits
{
    typedef List_Node_Traits< T > node_traits;
    typedef typename node_traits::node_ptr node_ptr;
    typedef typename node_traits::const_node_ptr const_node_ptr;
    typedef T value_type;
    typedef node_ptr pointer;
    typedef const_node_ptr const_pointer;
    typedef typename node_traits::fact_type::ref_type reference;
    typedef typename node_traits::fact_type::const_ref_type const_reference;

    static const boost::intrusive::link_mode_type link_mode = boost::intrusive::normal_link;

    static node_ptr to_node_ptr (reference value) { return &value; }
    static const_node_ptr to_node_ptr (const_reference value) { return &value; }
    static pointer to_value_ptr(node_ptr n) { return n; }
    static const_pointer to_value_ptr(const_node_ptr n) { return n; }
};


typedef Factory< B > fact_type;
typedef fact_type::ptr_type ptr_type;
typedef fact_type::const_ptr_type const_ptr_type;
typedef fact_type::ref_type ref_type;
typedef fact_type::const_ref_type const_ref_type;

typedef boost::intrusive::itree< Value_Traits< B > > iitree_type;
typedef iitree_type::itree_algo itree_algo;
typedef boost::intrusive::list< B, boost::intrusive::value_traits< List_Value_Traits< B > > > ilist_type;

const_ptr_type get_root(iitree_type& t)
{
    return itree_algo::get_header(&*t.begin())->_parent;
}

void print_sub_tree(const_ptr_type r, size_t depth)
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

bool check_max_ends(const_ptr_type node_ptr, size_t& max_end)
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
    const_ptr_type root_node = get_root(t);
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
    clog << "----- constructing factory\n";
    fact_type f(true);

    clog << "----- constructing iitree & ilist\n";
    iitree_type t;
    ilist_type l;

    iitree_type::intersection_iterator it = t.interval_intersect_end();
    boost::iterator_range<iitree_type::intersection_iterator> r1(it, it);
    (void)r1;

    typedef boost::iterator_range<iitree_type::intersection_iterator> it_range_type;
    it_range_type r2 = boost::make_iterator_range<iitree_type::intersection_iterator>(it, it);
    (void)r2;

    typedef boost::range_iterator< it_range_type >::type type_1;
    type_1 obj_1;
    typedef boost::range_difference< iitree_type::intersection_iterator >::type type_2;

    /*
    const size_t max_load = 100;
    const size_t range_max = 20;
    const size_t n_ops = 1000000;
    srand48(88473124);

    clog << "----- main loop\n";
    for (size_t i = 0; i < n_ops; ++i)
    {
        int op = int(drand48()*3);
        if (op == 0)
        {
            // insert new element
            if (l.size() >= max_load)
            {
                continue;
            }
            ptr_type a = f.new_elem();
            size_t e1 = size_t(drand48() * range_max);
            size_t e2 = size_t(drand48() * range_max);
            a->_start = min(e1, e2);
            a->_end = max(e1, e2);
            clog << "adding: " << *a << '\n';
            l.push_back(*a);
            t.insert(*a);
        }
        else if (op == 1)
        {
            // delete existing element
            if (l.size() == 0)
            {
                continue;
            }
            size_t idx = size_t(drand48() * l.size());
            auto it = l.begin();
            while (idx > 0)
            {
                ++it;
                --idx;
            }
            ptr_type a = &*it;
            clog << "deleting: " << *a << '\n';
            l.erase(it);
            t.erase(t.iterator_to(*a));
            f.del_elem(a);
        }
        else if (op == 2)
        {
            // compute intersection with some interval
            size_t e1 = size_t(drand48() * range_max);
            size_t e2 = size_t(drand48() * range_max);
            if (e1 > e2)
            {
                swap(e1, e2);
            }
            ptr_type a = f.new_elem();
            a->_start = e1;
            a->_end = e2;
            clog << "checking: " << *a << '\n';
            // first count intersections using list
            size_t res_list = 0;
            for (const auto& v : l)
            {
                if (intersect(v, *a))
                {
                    ++res_list;
                }
            }
            // next, count intersections using interval tree
            auto tmp_v = t.interval_intersect(e1, e2);
            size_t res_itree = tmp_v.size();
            // count also using intersection_iterator
            size_t res_iterator = 0;
            auto it_crt = t.interval_intersect_begin(e1, e2);
            auto it_end = t.interval_intersect_end();
            while (it_crt != it_end)
            {
                if (&*it_crt != tmp_v[res_iterator])
                {
                    break;
                }
                clog << "[" << it_crt->_start << "," << it_crt->_end << "]\n";
                ++it_crt;
                ++res_iterator;
            }
            // count with iterator range
            size_t res_iterator_range = 0;
            for (const auto& r : t.interval_intersect_range(e1, e2))
            {
                (void)r;
                ++res_iterator_range;
            }
            if (res_list != res_itree or res_list != res_iterator or res_list != res_iterator_range)
            {
                clog << "wrong intersection with " << *a << '\n';
                clog << "factory:\n" << f;
                clog << "list:\n";
                for (const auto& v : l)
                {
                    clog << v << '\n';
                }
                clog << "itree:\n";
                for (const auto& v : t)
                {
                    clog << v << '\n';
                }
                print_tree(t);
                exit(EXIT_FAILURE);
            }
            clog << "intersection ok, size = " << res_list << " / " << l.size() << '\n';
            f.del_elem(a);
        }
    }
    clog << "----- clearing list\n";
    // clear list
    while (l.size() > 0)
    {
        auto it = l.begin();
        ptr_type a = &*it;
        l.erase(it);
        t.erase(t.iterator_to(*a));
        f.del_elem(a);
    }
    if (f.unused() != f.size() - 3)
    {
        clog << "factory not empty on exit:\n" << f;
        exit(EXIT_FAILURE);
    }
    clog << "----- success\n";
    */
}
