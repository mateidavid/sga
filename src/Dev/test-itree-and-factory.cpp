#include "shortcuts.hpp"
#include <iostream>
#include <time.h>
#include <boost/intrusive/list.hpp>
#include <boost/intrusive/itree.hpp>
#include <boost/tti/tti.hpp>
#include <tclap/CmdLine.h>
#include "factory.hpp"
#include "ref_range.hpp"
#include "type_str.hpp"

using namespace std;
namespace bi = boost::intrusive;

namespace global
{
    string prog_desc =
        "Test intrusive itree and list with bounded pointers. "
        "The program tests the data structures by performing a series of random operations, such as:\n"
        "- insert new element\n"
        "- delete existing element\n"
        "- compute intersection with some interval\n"
        "- check max_end fields\n"
        "- clone tree\n";
    TCLAP::CmdLine cmd_parser(prog_desc);

    TCLAP::ValueArg< long > seed("", "seed", "Random seed (-1: use time).", false, -1, "int", cmd_parser);
    TCLAP::ValueArg< unsigned > max_load("", "max-load", "Maximum load of the interval tree.", false, 100, "int", cmd_parser);
    TCLAP::ValueArg< unsigned > range_max("", "range-max", "Maximum endpoint.", false, 20, "int", cmd_parser);
    TCLAP::ValueArg< unsigned > n_ops("", "n-ops", "Number of operations.", false, 1000, "int", cmd_parser);
    TCLAP::SwitchArg print_tree("", "print-tree", "Print tree after each operation.", cmd_parser, false);
}

namespace detail
{
    //check Boost Intrusive contains hooks to maintain extra data
    typedef bi::detail::extra_data_manager< void > extra_data_manager_check;
}

struct Value
{
//private:
//    friend class Factory< Value >;
    Value() = default;
    Value(const Value& rhs)
        : _start(rhs._start), _end(rhs._end),
          _parent(), _l_child(), _r_child(), _list_prev(), _list_next() {}
    Value(Value&&) = delete;
    ~Value() { ASSERT(is_unlinked()); }
//public:
    typedef bounded::Factory< Value > fact_type;
    typedef fact_type::ptr_type ptr_type;

    size_t _start;
    size_t _end;

    ptr_type _parent;
    ptr_type _l_child;
    ptr_type _r_child;
    int _col;
    size_t _max_end;

    ptr_type _list_prev;
    ptr_type _list_next;
    bool is_unlinked() const { return not(_parent or _l_child or _r_child or _list_prev or _list_next); }
};

typedef bounded::Factory< Value > fact_type;
typedef fact_type::ptr_type ptr_type;
typedef fact_type::const_ptr_type const_ptr_type;
typedef fact_type::ref_type ref_type;
typedef fact_type::const_ref_type const_ref_type;

ostream& operator <<(ostream& os, const Value& rhs)
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

bool intersect(const Value& lhs, const Value& rhs)
{
    return (lhs._start <= rhs._start and rhs._start <= lhs._end)
           or (rhs._start <= lhs._start and lhs._start <= rhs._end);
}

template <class T>
struct ITree_Node_Traits
{
    typedef T node;
    typedef bounded::Factory< T > fact_type;
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
struct ITree_Value_Traits
{
    typedef T value_type;
    typedef ITree_Node_Traits< T > node_traits;
    typedef typename node_traits::key_type key_type;
    typedef typename node_traits::node_ptr node_ptr;
    typedef typename node_traits::const_node_ptr const_node_ptr;
    typedef node_ptr pointer;
    typedef const_node_ptr const_pointer;
    typedef typename node_traits::fact_type::ref_type reference;
    typedef typename node_traits::fact_type::const_ref_type const_reference;
    typedef ITree_Value_Traits* value_traits_ptr;

    static const bi::link_mode_type link_mode = bi::safe_link;

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
    typedef T node;
    typedef bounded::Factory< T > fact_type;
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
    typedef T value_type;
    typedef List_Node_Traits< T > node_traits;
    typedef typename node_traits::node_ptr node_ptr;
    typedef typename node_traits::const_node_ptr const_node_ptr;
    typedef node_ptr pointer;
    typedef const_node_ptr const_pointer;
    typedef typename node_traits::fact_type::ref_type reference;
    typedef typename node_traits::fact_type::const_ref_type const_reference;
    typedef List_Value_Traits* value_traits_ptr;

    static const bi::link_mode_type link_mode = bi::safe_link;

    static node_ptr to_node_ptr (reference value) { return &value; }
    static const_node_ptr to_node_ptr (const_reference value) { return &value; }
    static pointer to_value_ptr(node_ptr n) { return n; }
    static const_pointer to_value_ptr(const_node_ptr n) { return n; }
};

typedef bi::itree< Value,
                   bi::value_traits< ITree_Value_Traits< Value > >,
                   bi::header_holder_type< bounded::Pointer_Holder< Value > >
                 > itree_type;
typedef itree_type::itree_algo itree_algo;
typedef bi::list< Value,
                  bi::value_traits< List_Value_Traits< Value > >,
                  bi::header_holder_type< bounded::Pointer_Holder< Value > >
                  > list_type;

static_assert(
    bi::detail::extra_data_manager<
        bi::detail::itree_node_traits < ITree_Value_Traits< Value > >
    >::enabled,
    "Extra data manager is not enabled");

// Check iterator traits
// ITree::iterator
static_assert(std::is_same<
              std::iterator_traits< itree_type::iterator >::value_type,
              Value
              >::value, "iterator_traits<itree_type::iterator>::value_type != Value");
static_assert(std::is_same<
              std::iterator_traits< itree_type::iterator >::reference,
              bounded::Reference< Value >
              >::value, "iterator_traits<itree_type::iterator>::reference != Reference<Value>");
static_assert(std::is_same<
              std::iterator_traits< itree_type::iterator >::pointer,
              bounded::Pointer< Value >
              >::value, "iterator_traits<itree_type::iterator>::pointer != Pointer<Value>");
// ITree::const_iterator
// NOTE: value_type does not inherit const-ness, only pointer&reference types
static_assert(std::is_same<
              std::iterator_traits< itree_type::const_iterator >::value_type,
              Value
              >::value, "iterator_traits<itree_type::const_iterator>::value_type != Value");
static_assert(std::is_same<
              std::iterator_traits< itree_type::const_iterator >::reference,
              bounded::Reference< const Value >
              >::value, "iterator_traits<itree_type::const_iterator>::reference != Reference<const Value>");
static_assert(std::is_same<
              std::iterator_traits< itree_type::const_iterator >::pointer,
              bounded::Pointer< const Value >
              >::value, "iterator_traits<itree_type::const_iterator>::pointer != Pointer<const Value>");
// ITree::intersection_iterator
static_assert(std::is_same<
              std::iterator_traits< itree_type::intersection_iterator >::value_type,
              Value
              >::value, "iterator_traits<itree_type::intersection_iterator>::value_type != Value");
static_assert(std::is_same<
              std::iterator_traits< itree_type::intersection_iterator >::reference,
              bounded::Reference< Value >
              >::value, "iterator_traits<itree_type::intersection_iterator>::reference != Reference<Value>");
static_assert(std::is_same<
              std::iterator_traits< itree_type::intersection_iterator >::pointer,
              bounded::Pointer< Value >
              >::value, "iterator_traits<itree_type::intersection_iterator>::pointer != Pointer<Value>");
// ITree::intersection_const_iterator
static_assert(std::is_same<
              std::iterator_traits< itree_type::intersection_const_iterator >::value_type,
              Value
              >::value, "iterator_traits<itree_type::intersection_const_iterator>::value_type != Value");
static_assert(std::is_same<
              std::iterator_traits< itree_type::intersection_const_iterator >::reference,
              bounded::Reference< const Value >
              >::value, "iterator_traits<itree_type::intersection_const_iterator>::reference != Reference<const Value>");
static_assert(std::is_same<
              std::iterator_traits< itree_type::intersection_const_iterator >::pointer,
              bounded::Pointer< const Value >
              >::value, "iterator_traits<itree_type::intersection_const_iterator>::pointer != Pointer<const Value>");

const_ptr_type get_root(itree_type& t)
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

void print_tree(itree_type& t)
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

void check_max_ends(itree_type& t, fact_type&)
{
    const_ptr_type root_node = get_root(t);
    size_t max_end;
    if (not check_max_ends(root_node, max_end))
    {
        //clog << "factory:\n" << f;
        clog << "tree:\n";
        for (auto const& e :t)
        {
            clog << e << '\n';
        }
        print_tree(t);
        exit(1);
    }
}

void real_main()
{
    clog << "----- program options:"
         << "\nmax_load=" << global::max_load
         << "\nrange_max=" << global::range_max
         << "\nn_ops=" << global::n_ops
         << "\nseed=" << global::seed << '\n';

    clog << "--- type names:"
         << "\nValue: " << type_str< Value >()
         << "\nconst Value: " << type_str< const Value >()
         << "\nFactory< Value >: " << type_str< bounded::Factory< Value > >()
         << "\nPointer< Value >: " << type_str< bounded::Pointer< Value > >()
         << "\nPointer< const Value >: " << type_str< bounded::Pointer< const Value > >()
         << "\nReference< Value >: " << type_str< bounded::Reference< Value > >()
         << "\nReference< const Value >: " << type_str< bounded::Reference< const Value > >()
         << "\niterator_traits<intersection_iterator>::pointer: " << type_str< std::iterator_traits< itree_type::intersection_iterator >::pointer >()
         << '\n';

    clog << "--- type sizes:\n"
         << "sizeof(val_type)=" << sizeof(fact_type::val_type) << '\n'
         << "sizeof(ptr_type)=" << sizeof(ptr_type) << '\n'
         << "sizeof(ref_type)=" << sizeof(fact_type::ref_type) << '\n';

    clog << "----- constructing factory\n";
    fact_type f(true);

    {
        clog << "----- constructing itree & list\n";
        itree_type t;
        list_type l;

        clog << "----- initializing random number generator\n";
        srand48(global::seed);

        clog << "----- main loop\n";
        for (size_t i = 0; i < global::n_ops; ++i)
        {
            int op = int(drand48()*5);
            if (op == 0)
            {
                // insert new element
                if (l.size() >= global::max_load)
                {
                    continue;
                }
                ptr_type a = f.new_elem();
                size_t e1 = size_t(drand48() * global::range_max);
                size_t e2 = size_t(drand48() * global::range_max);
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
                size_t e1 = size_t(drand48() * global::range_max);
                size_t e2 = size_t(drand48() * global::range_max);
                if (e1 > e2)
                {
                    swap(e1, e2);
                }
                ptr_type a = f.new_elem();
                a->_start = e1;
                a->_end = e2;
                clog << "checking intersection with: " << *a << '\n';
                // first count intersections using list
                size_t res_list = 0;
                for (const auto& v : l)
                {
                    if (intersect(v, *a))
                    {
                        ++res_list;
                    }
                }
                // count with iterator range
                size_t res_iterator_range = 0;
                for (const auto& r : t.iintersect(e1, e2))
                {
                    (void)r;
                    ++res_iterator_range;
                }
                if (res_iterator_range != res_list)
                {
                    clog << "wrong intersection with " << *a << '\n';
                    //clog << "factory:\n" << f;
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
                    abort();
                }
                clog << "intersection ok, size = " << res_list << " / " << l.size() << '\n';
                f.del_elem(a);
            }
            else if (op == 3)
            {
                // check max_end fields in the tree
                clog << "checking max_end fields\n";
                check_max_ends(t, f);
            }
            else if (op == 4)
            {
                // clone tree and check
                clog << "cloning tree of size: " << t.size() << '\n';
                itree_type t2;
                t2.clone_from(t, fact_type::cloner_type(), fact_type::disposer_type()); //, fact_type::cloner_type(), fact_type::disposer_type());
                clog << "checking max_end fields in clone of size: " << t2.size() << '\n';
                check_max_ends(t2, f);
                clog << "destroying clone\n";
                ptr_type tmp;
                t2.clear_and_dispose(fact_type::disposer_type());
            }
            if (global::print_tree)
            {
                print_tree(t);
            }
        }
        clog << "----- clearing list\n";
        t.clear();
        l.clear_and_dispose(fact_type::disposer_type());
    }
    if (f.unused() != f.size())
    {
        clog << "factory not empty on exit!\n"; // << f;
        abort();
    }
    clog << "----- success\n";
}

int main(int argc, char* argv[])
{
    global::cmd_parser.parse(argc, argv);
    global_assert::prog_name() = global::cmd_parser.getProgramName();
    if (global::seed < 0)
    {
        global::seed.get() = time(nullptr);
    }
    srand48(global::seed);
    real_main();
}
