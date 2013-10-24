#include <iostream>
#include <vector>
#include "print_seq.hpp"

using namespace std;


int main()
{
    vector<int> v(5);
    vector<int*> u(2);
    u[0] = new int(5);
    u[1] = new int(7);

    cout << "v: "; print_seq(cout, v); cout << "\n";
    cout << "u: "; print_seq(cout, u); cout << "\n";
    cout << "u*: (\n  "; print_ptr_seq(cout, u, "\n  "); cout << "\n)\n";

    return 0;
}
