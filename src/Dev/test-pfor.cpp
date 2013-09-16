#include "pfor.hpp"

using namespace std;


int sum;

void process_item(TLS_pfor& tls, int& e)
{
#pragma omp atomic
    sum += e;
    if (e % 5 == 0 and e % 10 != 0)
        *tls.out_str_p << e << "\n";
}

bool get_item(TLS_pfor& tls, int& e)
{
    return cin >> e;
}

int main()
{
    pfor<int, TLS_pfor>(NULL,
                        //[&cin] (TLS_pfor<int>& tls, int& e) {
                        //    return bool(cin >> e);
                        //},
                        &get_item,
                        &process_item,
                        NULL,
                        4,
                        100);
    cout << sum << "\n";

    return 0;
}
