#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <typeinfo>

using namespace std;
using namespace boost;

struct Vertex
{
    string my_id;
};

struct Edge
{
    int my_val;
};

typedef adjacency_list<listS, listS, undirectedS, Vertex, Edge> g_type;
g_type g;

typedef property_map<g_type, int Edge::*>::type pm_type;


void print_graph(g_type g)
{
    graph_traits<g_type>::vertex_iterator v_it, v_it_end;
    graph_traits<g_type>::out_edge_iterator e_it, e_it_end;

    for (tie(v_it, v_it_end) = vertices(g); v_it != v_it_end; ++v_it)
    {
        cout << "vertex descriptor=[" << *v_it << "] my_id=[" << g[*v_it].my_id << "]\n";
        for (tie(e_it, e_it_end) = out_edges(*v_it, g); e_it != e_it_end; ++e_it)
        {
            cout << "edge descriptor=[" << *e_it << "] my_val=[" << g[*e_it].my_val << "]"
                 << " source=[" << source(*e_it, g) << "]"
                 << " target=[" << target(*e_it, g) << "]\n";
        }
    }
}


int main()
{
    graph_traits<g_type>::vertex_descriptor u, v;
    graph_traits<g_type>::edge_descriptor e, f;
    graph_traits<g_type>::vertex_iterator v_it, v_it_end;
    graph_traits<g_type>::out_edge_iterator e_it, e_it_end;
    bool ok;

    u = add_vertex(g);
    g[u].my_id = "One";

    v = add_vertex(g);
    g[v].my_id = "Two";

    Edge e_p;
    e_p.my_val = 5;
    tie(e, ok) = add_edge(u, v, e_p, g);
    if (not ok)
        cout << "edge already present; returned descriptor=[" << e << "]\n";
    f = e;

    pm_type pm = get(&Edge::my_val, g);
    //pm[f] = 17;
    //pm = 0;
    u = .3;

    cout << "print 1:\n";
    print_graph(g);

    cout << "typeof(pm): " << typeid(pm).name() << "\n";
    cout << "pm[f]: " << pm[f] << "\n";

    v = add_vertex(g);
    g[v].my_id = "Three";

    e_p.my_val = 7;
    tie(e, ok) = add_edge(u, v, e_p, g);
    if (not ok)
        cout << "edge already present; returned descriptor=[" << e << "]\n";
    g[f].my_val = 9;

    cout << "print 2:\n";
    print_graph(g);
    cout << "pm[f]: " << pm[f] << "\n";

    return EXIT_SUCCESS;
}
