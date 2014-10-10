#include <iostream>
#include <vector>
#include "ref_iterator.hpp"


int main()
{
    typedef std::vector< int > cont_t;
    typedef cont_t::iterator it_t;
    typedef cont_t::const_iterator cit_t;
    typedef cont_t::reverse_iterator rit_t;
    typedef cont_t::const_reverse_iterator crit_t;
    typedef ref_iterator< it_t > ref_it_t;
    typedef ref_iterator< cit_t > ref_cit_t;
    typedef ref_iterator< rit_t > ref_rit_t;
    typedef ref_iterator< crit_t > ref_crit_t;

    std::vector< int > v = { 1, 2, 3 };

    for (ref_it_t it = ref_it_t(v.begin()); it != ref_it_t(v.end()); ++it)
    {
        std::cout << *it << "\n";
    }

    for (auto it = ref_rit_t(v.rbegin()); it != ref_rit_t(v.rend()); ++it)
    {
        std::cout << *it << "\n";
    }
}
