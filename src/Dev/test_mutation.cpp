#include "MAC.hpp"

#include <vector>

using namespace std;


int main()
{
    typename MAC::Mutation m1;
    MAC::Mutation m2(3, 1);
    MAC::Mutation m3(5, 1, string("C"));

    vector<MAC::Mutation> v(3);
    v[0] = MAC::Mutation::ins(5, 'A');
    v[1] = MAC::Mutation::snp(7, 'T');
    v[2] = MAC::Mutation::del(20);

    cout << "m1:" << m1 << "\n";
    cout << "m2:" << m2 << "\n";
    cout << "m3:" << m3 << "\n";
    for (size_t i = 0; i < v.size(); ++i)
    {
        cout << "v[" << i << "]:" << v[i] << "\n";
    }

    return EXIT_SUCCESS;
}
