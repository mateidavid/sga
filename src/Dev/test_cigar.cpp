#include <iostream>

#include "Cigar.hpp"

using namespace std;
using namespace MAC;


int main()
{
    string s;
    cin >> s;
    bool comp;
    cin >> comp;
    Cigar c(s, comp);

    cout << c << "\n";

    Cigar c2 = c.complement();
    cout << "complement:" << c2 << "\n";

    Cigar c3 = c2.complement();
    cout << "2x complement:" << c3 << "\n";

    cout << "c == c3? : " << (c == c3) << "\n";

    cout << "substring(c,0,0):" << c.substring(0, 0) << "\n";
    cout << "substring(c,0,len):" << c.substring(0, c.get_n_ops()) << "\n";
    cout << "substring(c,1,2):" << c.substring(1, 3) << "\n";

    c.cut_op(1,1);
    cout << "cut_op(1,1):" << c << "\n";

    return 0;
}
