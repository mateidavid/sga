#include <iostream>
#include <vector>
#include "ref_range.hpp"


int main()
{
    const std::vector< int > v = { 1, 2, 3 };
    for (auto e : v | referenced)
    {
        std::cout << *e << "\n";
    }
}
