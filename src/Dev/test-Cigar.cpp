#include <iostream>

#include "Cigar.hpp"

using namespace std;
using namespace cigar;


int main()
{
    typedef Cigar< unsigned > cigar_type;
    string s;
    cin >> s;
    bool comp;
    cin >> comp;
    cigar_type c(s, comp);

    cout << c << "\n";

    cigar_type c2 = c.complement();
    cout << "complement:" << c2 << "\n";

    cigar_type c3 = c2.complement();
    cout << "2x complement:" << c3 << "\n";

    cout << "c == c3? : " << (c == c3) << "\n";

    cout << "substring(c,0,0):" << c.subcigar(0, 0) << "\n";
    cout << "substring(c,0,len):" << c.subcigar(0, c.n_ops()) << "\n";
    cout << "substring(c,1,2):" << c.subcigar(1, 2) << "\n";

    c.cut_op(1,1);
    cout << "cut_op(1,1):" << c << "\n";

    return 0;
}
