#include <iostream>
#include <vector>
#include <list>

#include "MAC.hpp"

using namespace std;
using namespace MAC;


int main()
{
    Read_Entry_Cont re_cont;
    Contig_Entry_Cont ce_cont;

    add_read(re_cont, ce_cont, new string("001"), new Seq_Type("ACGT"));

    cout << *(ce_cont.begin());
    cout << *(re_cont.begin());

    return EXIT_SUCCESS;
}
