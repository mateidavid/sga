#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <list>

#include "Graph.hpp"
#include "indent.hpp"
#include "ixstream.hpp"

using namespace std;
using namespace MAC;

auto inc_tab = indent::inc;
auto dec_tab = indent::dec;
using indent::nl;
using indent::tab;


int main(int argc, char* argv[])
{
    Graph g;
    assert(argc >= 2);
    ixstream ixs(argv[1]);
    if (not ixs)
    {
        cerr << "error opening file: " << argv[1] << '\n';
        return EXIT_FAILURE;
    }

    string line;
    size_t line_count = 0;
    while (getline(ixs, line))
    {
        cout << line << '\n';
        istringstream iss(line + "\n");
        string rec_type;
        iss >> rec_type;
        if (rec_type == "HT")
        {
            // ignore header line for now
        }
        else if (rec_type == "VT")
        {
            string* r_id_ptr = new string();
            string* r_seq_ptr = new string();
            iss >> *r_id_ptr >> *r_seq_ptr;
            assert(not iss.eof());
            g.add_read(r_id_ptr, r_seq_ptr);
        }
        else if (rec_type == "ED")
        {
            string r1_id;
            string r2_id;
            size_t r1_start, r1_end, r1_len;
            size_t r2_start, r2_end, r2_len;
            int rc, tmp;
            string sam_cigar;
            string sam_pi;
            iss >> r1_id >> r2_id >> r1_start >> r1_end >> r1_len >> r2_start >> r2_end >> r2_len >> rc >> tmp >> sam_cigar >> sam_pi;
            assert(not iss.eof());
            ++r1_end;
            ++r2_end;
            g.add_overlap(r1_id, r2_id, r1_start, r1_end - r1_start, r2_start, r2_end - r2_start, rc, sam_cigar.substr(5));
        }
        //cerr << g;
        if ((++line_count % 10000) == 0)
            cerr << line_count << '\n';
    }
    cout << g;
    cerr << "success\n";

    return EXIT_SUCCESS;
}
