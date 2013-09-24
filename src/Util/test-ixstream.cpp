#include <iostream>
#include "ixstream.hpp"
using namespace std;

int main(int argc, char * argv[])
{
    ixstream s1;
    ixstream s2(cin);
    ixstream s3(argv[1]);

    int i;
    while (s1 >> i)
        cout << "s1: " << i << "\n";
    while (s2 >> i)
        cout << "s2: " << i << "\n";
    while (s3 >> i)
        cout << "s3: " << i << "\n";

    return 0;
}
