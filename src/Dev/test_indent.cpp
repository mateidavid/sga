#include <iostream>
#include "indent.hpp"

using namespace std;


int main()
{
    cout << "level 0";
    cout << indent::inc << indent::endl << "level 1";
    cout << indent::endl << "still level 1";
    cout << indent::inc << indent::nl << "level 2";
    cout << indent::dec << indent::nl << "level 1 again";
    cout << '\n' << "not indented";
    cout << indent::endl << "indented again";
    cout << indent::dec << indent::nl << "done indenting\n";
    return 0;
}
