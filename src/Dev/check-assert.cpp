#include <iostream>
#include <cassert>
#include "global_assert.hpp"

int x;

bool f() { x += 1; return true; }

int main()
{
    x = 0;
    assert(f());
    std::cout << "assert: " << (x > 0? "enabled" : "disabled") << "\n";

    x = 0;
    ASSERT(f());
    std::cout << "ASSERT: " << (x > 0? "enabled" : "disabled") << "\n";

    x = 0;
    BOOST_ASSERT_MSG(f(), "");
    std::cout << "BOOST_ASSERT_MSG: " << (x > 0? "enabled" : "disabled") << "\n";
}
