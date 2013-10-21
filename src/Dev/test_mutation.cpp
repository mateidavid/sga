#include "Mutation.hpp"

#include <vector>

using namespace std;


int main()
{
    Mutation<string,size_t> m1;
    Mutation<string> m2(3, 1);
    Mutation<> m3(5, 1, string("C"));

    vector<Mutation<string, unsigned char>> v(3);
    v[0] = Mutation<string, unsigned char>::ins(5, 'A');
    v[1] = Mutation<string, unsigned char>::snp(7, 'T');
    v[2] = Mutation<string, unsigned char>::del(20);

    cout << "m1:" << m1 << "\n";
    cout << "m2:" << m2 << "\n";
    cout << "m3:" << m3 << "\n";
    for (size_t i = 0; i < v.size(); ++i)
    {
        cout << "v[" << i << "]:" << v[i] << "\n";
    }

    return EXIT_SUCCESS;
}
