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

    g.add_read(new string("001"), new Seq_Type("ACGTGGGCATTTGACATAGCTGGATAAACCCCTCAGATCGGAGTTACTTA"));
    cerr << indent::nl << "g after adding read 001:" << indent::inc << indent::nl << g << indent::dec << indent::nl;
    g.check();

    g.add_read(new string("002"), new Seq_Type("AAACCCTCAGATCGCAGTTAACTTAGGTCAGAGTTTAGCTTGATCCTTTA"));
    cerr << indent::nl << "g after adding read 002:" << indent::inc << indent::nl << g << indent::dec << indent::nl;
    g.check();

    g.add_overlap(string("001"), string("002"), 25, 25, 0, 25, 0, string("6M1D14M1I4M"));
    cerr << indent::nl << "g after adding overlap 001-002:" << indent::inc << indent::nl << g << indent::dec << indent::nl;
    g.check();

    return EXIT_SUCCESS;
}
