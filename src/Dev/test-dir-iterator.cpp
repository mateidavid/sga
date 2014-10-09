#include <iostream>
#include <vector>
#include "dir_iterator.hpp"


int main()
{
    std::vector<int> v = { 1, 2, 3 };
    typedef dir_iterator< std::vector<int>::iterator > dir_it_t;

    for (dir_it_t it(v.begin()); it != dir_it_t(v.end()); ++it)
        std::cout << *it << "\n";
    for (dir_it_t it(v.end(), true); it != dir_it_t(v.begin(), true); ++it)
        std::cout << *it << "\n";
}
