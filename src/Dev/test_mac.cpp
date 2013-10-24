#include <iostream>
#include <vector>
#include <list>

#include "MAC.hpp"

using namespace std;
using namespace MAC;


int main()
{
    vector<Read_Entry> re_cont;
    Contig_Entry_Cont ce_cont;
    Contig_Entry* tmp;

    re_cont.push_back(Read_Entry());
    ce_cont.push_back(Contig_Entry());
    tmp = const_cast<Contig_Entry*>(&(*ce_cont.begin()));
    add_read(string("001"), Seq_Type("ACGT"), &re_cont[0], tmp);

    cout << *(ce_cont.begin());
    cout << re_cont[0];

    return EXIT_SUCCESS;
}
