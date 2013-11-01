#include <iostream>
#include <vector>
#include <list>

#include "Graph.hpp"
#include "indent.hpp"

using namespace std;
using namespace MAC;


int main()
{
    Graph g;

    g.add_read(new string("001"), new Seq_Type("ACGT"));
    cerr << "g after add_read:" << indent::inc << indent::nl << g << indent::dec << indent::nl;

    return EXIT_SUCCESS;
}
