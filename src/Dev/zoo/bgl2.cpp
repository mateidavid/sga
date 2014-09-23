#include <iostream>
#include <boost/graph/adjacency_list.hpp>

using namespace std;
using namespace boost;

struct Vertex
{
    Vertex() {}
    Vertex(const string& name_) : name(name_) {}
    string name;
};

namespace boost
{
    namespace graph
    {
        template<>
        struct internal_vertex_name<Vertex>
        {
            typedef multi_index::member<Vertex, string, &Vertex::name> type;
        };
        
        template<>
        struct internal_vertex_constructor<Vertex>
        {
            typedef vertex_from_name<Vertex> type;
        };
    }
    
    struct unique_listS {};
    template <class T>
    struct container_gen<unique_listS, T> {
        typedef list<T> type;
    };
    template <>
    struct parallel_edge_traits<unique_listS> { 
        typedef disallow_parallel_edge_tag type;
    };
}


int main()
{
    typedef adjacency_list<unique_listS, listS, undirectedS, Vertex> g_type;
    
    g_type g;
    graph_traits<g_type>::vertex_descriptor u, v;
    graph_traits<g_type>::edge_descriptor e, f;
    graph_traits<g_type>::vertex_iterator v_it, v_it_end;
    graph_traits<g_type>::out_edge_iterator e_it, e_it_end;
    bool ok;

    u = add_vertex("One", g);
    v = add_vertex("Two", g);

    tie(e, ok) = add_edge("One", "Two", g);
    if (not ok)
        cout << "edge One->Two already present\n";
    cout << "returned descriptor=[" << e << "]\n";

    tie(e, ok) = add_edge(v, u, g);
    if (not ok)
        cout << "edge Two->One already present\n";
    cout << "returned descriptor=[" << e << "]\n";

    tie(e, ok) = add_edge("Two", "Three", g);
    if (not ok)
        cout << "edge Two->Three already present\n";
    cout << "returned descriptor=[" << e << "]\n";
    
    for (tie(v_it, v_it_end) = vertices(g); v_it != v_it_end; ++v_it)
    {
        cout << "vertex descriptor=[" << *v_it << "] name=[" << g[*v_it].name << "] -->";
        for (tie(e_it, e_it_end) = out_edges(*v_it, g); e_it != e_it_end; ++e_it)
        {
            cout << " " << *e_it;
        }
        cout << "\n";
    }

    return EXIT_SUCCESS;
}
