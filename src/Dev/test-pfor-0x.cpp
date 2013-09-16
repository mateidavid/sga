#include "pfor.hpp"

using namespace std;


int sum;

int main()
{
    pfor<int, TLS_pfor>(NULL,
                        [&cin] (TLS_pfor& tls, int& e) { return bool(cin >> e); },
                        [] (TLS_pfor& tls, int& e) { *tls.out_str_p << e << "\n"; },
                        NULL,
                        4,
                        10);
    cout << sum << "\n";

    return 0;
}
