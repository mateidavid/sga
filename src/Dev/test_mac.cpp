#include <iostream>
#include <vector>
#include <list>

#include "Graph.hpp"
#include "indent.hpp"

using namespace std;
using namespace MAC;

auto inc_tab = indent::inc;
auto dec_tab = indent::dec;
using indent::nl;
using indent::tab;


int main()
{
    Graph g;

    string line;
    while (getline(cin, line))
    {
        cerr << line << '\n';
        istringstream is(line + "\n");
        string rec_type;
        is >> rec_type;
        if (rec_type == "HT")
        {
            // ignore header line for now
        }
        else if (rec_type == "VT")
        {
            string* r_id_ptr = new string();
            string* r_seq_ptr = new string();
            is >> *r_id_ptr >> *r_seq_ptr;
            assert(not is.eof());
            g.add_read(r_id_ptr, r_seq_ptr);
            g.check();
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
            is >> r1_id >> r2_id >> r1_start >> r1_end >> r1_len >> r2_start >> r2_end >> r2_len >> rc >> tmp >> sam_cigar >> sam_pi;
            assert(not is.eof());
            ++r1_end;
            ++r2_end;
            g.add_overlap(r1_id, r2_id, r1_start, r1_end - r1_start, r2_start, r2_end - r2_start, rc, sam_cigar.substr(5));
            g.check();
        }
        cerr << g;
    }

    return EXIT_SUCCESS;
}
