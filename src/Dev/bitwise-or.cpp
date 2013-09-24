#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdint.h>

using namespace std;


// Compute bitwise OR of two files.


int main(int argc, char * argv[])
{
    if (argc != 3)
    {
        cerr << "use: " << argv[0] << " <file_1> <file_2>\n";
        exit(EXIT_FAILURE);
    }

    ifstream is1(argv[1], ios_base::binary);
    ifstream is2(argv[2], ios_base::binary);

    uint64_t tmp1;
    uint64_t tmp2;
    do
    {
        tmp1 = 0;
        is1.read((char*)&tmp1, 8);
        tmp2 = 0;
        is2.read((char*)&tmp2, 8);
        if (is1.gcount() != is2.gcount())
        {
            cerr << "inputs were of different length\n";
            exit(EXIT_FAILURE);
        }
        tmp1 |= tmp2;
        cout.write((char*)&tmp1, is1.gcount());
    } while (is1.good() and is2.good());

    return EXIT_SUCCESS;
}
