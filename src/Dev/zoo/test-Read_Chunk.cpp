#include <iostream>
#include <string>
#include "Read_Chunk.hpp"
#include "print_seq.hpp"

using namespace std;
using namespace MAC;


int main()
{
    cout << "cigar: ";
    string cigar_string;
    cin >> cigar_string;

    cout << "cigar reversed: ";
    int cigar_reversed;
    cin >> cigar_reversed;

    cout << "rf_start: ";
    Size_Type rf_start;
    cin >> rf_start;

    cout << "qr_start: ";
    Size_Type qr_start;
    cin >> qr_start;

    cout << "query: ";
    string qr;
    cin >> qr;

    Cigar cigar(cigar_string, cigar_reversed, rf_start, qr_start);

    cout << cigar << "\n";

    Read_Chunk chunk;
    Mutation_Cont mut_cont;
    tie(chunk, mut_cont) = Read_Chunk::make_chunk_from_cigar(cigar, qr);

    cout << chunk << "\n";
    cout << "Mutation_Cont:" << indent::inc;
    print_seq(cout, mut_cont, indent::nl, indent::nl);
    cout << indent::dec << indent::nl;

    Read_Chunk_Pos pos = chunk.get_start_pos();
    bool forward = true;
    while (forward or not (pos == chunk.get_start_pos()))
    {
        if (forward and pos == chunk.get_end_pos())
        {
            cout << "reached end; going backwards\n";
            forward = false;
        }
        cout << "position: " << pos << "\n";
        cout << "direction: " << (forward? "forward" : "backward") << "\n";

        cout << "breakpoint: ";
        Size_Type brk;
        cin >> brk;
        bool on_contig = true;
        if (brk > 0)
        {
            cout << "on_contig: ";
            cin >> on_contig;
        }
        if (forward)
            chunk.increment_pos(pos, brk, on_contig);
        else
            chunk.decrement_pos(pos, brk, on_contig);
    }

    return 0;
}
