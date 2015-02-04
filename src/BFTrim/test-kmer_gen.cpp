#include <iostream>
#include <string>
#include "kmer_gen.hpp"

using namespace std;

int main()
{
    string s = "abracadara";
    size_t k = 3;
    auto kg = kmer_gen(s.size(), k);
    for (auto it = kg.begin(); it != kg.end(); ++it)
    {
        cout << *it << ": " << s.substr(*it, k) << endl;
    }
    kg = kmer_gen(s.size(), k, [&] (size_t i) { return s[i] != 'r'; });
    for (auto it = kg.begin(); it != kg.end(); ++it)
    {
        cout << *it << ": " << s.substr(*it, k) << endl;
    }
}
