#include <iostream>
#include <cstdlib>
#include "bithacks.hpp"

using namespace std;


// Count set bits on stdin.
//
// Report: bytes count, set bit count, and pct of set bits.


int main()
{
    size_t byte_cnt = 0;
    size_t setbit_cnt = 0;

    uint64_t tmp;
    do
    {
        tmp = 0;
        cin.read((char*)&tmp, 8);
        byte_cnt += cin.gcount();
        setbit_cnt += bitcount(tmp);
    } while (cin.good());

    cout << byte_cnt << "\t" << setbit_cnt
         << "\t" << double(setbit_cnt)/(8*byte_cnt) << "\n";

    return EXIT_SUCCESS;
}
