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

    g.add_read(new string("001"), new Seq_Type("ACGTGGGCATTTGACATAGCTGGATAAACCCTCAGATCGGAGTTACTTAG"));
    cerr << "g after adding read 001:" << indent::inc << indent::nl << g << indent::dec << indent::nl;
    g.add_read(new string("002"), new Seq_Type("AAACCCTCAGATCGGAGTTACTTAGGTCAGAGTTTAGCTTGATCCTTTAG"));
    cerr << "g after adding read 002:" << indent::inc << indent::nl << g << indent::dec << indent::nl;
    g.add_overlap(string("001"), string("002"), 25, 25, 0, 25, 0, string("25M"));
    cerr << "g after adding overlap 001-002:" << indent::inc << indent::nl << g << indent::dec << indent::nl;

    return EXIT_SUCCESS;
}
