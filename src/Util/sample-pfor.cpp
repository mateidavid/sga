#include <iostream>
#include <vector>

#include "pfor.hpp"

using namespace std;

struct Thread_Output_Storage
{
    Thread_Output_Storage() {}
    ~Thread_Output_Storage()
    {
        cout << std::hash<std::thread::id>()(std::this_thread::get_id()) << ": " << os.str() << flush;
    }
    ostringstream os;
}; // struct Thread_Output_Storage


int main()
{
    vector< int > v = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    auto it = v.begin();
    auto it_end = v.end();
    pfor< int, Thread_Output_Storage >(
        nullptr,
        [&] (int& i)
        {
            if (it == it_end)
            {
                return false;
            }
            i = *(it++);
            return true;
        },
        [&] (int& i, Thread_Output_Storage& tos)
        {
            tos.os << i << endl;
        },
        nullptr,
        4,
        1);
}
