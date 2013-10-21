#include <iostream>
#include <vector>
#include <list>

#include "MAC.hpp"

using namespace std;
using namespace MAC;


int main()
{
    vector<Read_Entry> re_list;
    vector<Contig_Entry> ce_list;

    re_list.push_back(Read_Entry());
    ce_list.push_back(Contig_Entry());
    add_read(string("001"), Seq_Type("ACGT"), &re_list[0], &ce_list[0]);

    cout << ce_list[0];
    cout << re_list[0];

    return EXIT_SUCCESS;
}
