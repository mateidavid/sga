#include "Cigar_Op.hpp"

using cigar::Cigar_Op;
using std::string;

int main()
{
    Cigar_Op< size_t > c1;
    Cigar_Op< uint8_t > c2;
    string cigar_string("1M1I1=10S200H");
    std::istringstream is(cigar_string.c_str());
    while (is >> c1)
    {
        std::cout << "c1='" << c1 << "'\n";
    }
    is.clear();
    is.str(cigar_string.c_str());
    while (is >> c2)
    {
        std::cout << "c2='" << c2 << "'\n";
    }
}
