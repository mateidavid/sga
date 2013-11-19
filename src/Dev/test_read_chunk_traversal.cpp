#include <iostream>
#include <string>
#include <tuple>
#include "Read_Chunk.hpp"
#include "print_seq.hpp"

using namespace std;
using namespace MAC;

auto inc_tab = indent::inc;
auto dec_tab = indent::dec;
using indent::nl;


void test_traversal(const Read_Chunk& chunk, Size_Type brk, bool on_contig)
{
    cout << nl << "testing brk=" << brk << ",on_contig=" << on_contig << inc_tab;

    vector< Read_Chunk_Pos > v;
    Read_Chunk_Pos pos = chunk.get_start_pos();
    bool reached_breakpoint = false;
    bool forward = true;
    while (true)
    {
        if (forward and pos == chunk.get_start_pos())
        {
            cout << dec_tab << nl << "at start, going forward" << inc_tab;
        }

        if (not forward and pos == chunk.get_end_pos())
        {
            cout << dec_tab << nl << "at end, going backward" << inc_tab;
        }

        cout << nl << "pos: " << pos;

        if (not reached_breakpoint and ((on_contig and pos.c_pos == brk) or (not on_contig and pos.r_pos == brk)))
        {
            cout << nl << "breakpoint reached";
            reached_breakpoint = true;
        }

        if (forward)
        {
            v.push_back(pos);
        }
        else
        {
            if (not (pos == v.back()))
            {
                cerr << "\nmismatched positions:\nforward: " << v.back() << "\nbackward: " << pos << "\n";
                exit(1);
            }
            v.resize(v.size() - 1);
        }

        if (not forward and pos == chunk.get_start_pos())
        {
            cout << dec_tab << nl << "at start, going backward: done";
            break;
        }

        if (forward and pos == chunk.get_end_pos())
        {
            cout << dec_tab << nl << "at end, forward: switching direction" << inc_tab;
            forward = false;
            reached_breakpoint = false;
            continue;
        }

        if (forward)
            chunk.increment_pos(pos, reached_breakpoint? 0 : brk, on_contig);
        else
            chunk.decrement_pos(pos, reached_breakpoint? 0 : brk, on_contig);
    }
    if (v.size() > 0)
    {
        cerr << "\nposition vector not empty:";
        print_seq(cerr, v, nl, nl);
        exit(1);
    }
}


int main()
{
    string cigar_string;
    Size_Type rf_start;
    Size_Type qr_start;

    while (cin >> cigar_string >> rf_start >> qr_start)
    {
        for (int cigar_reversed = 0; cigar_reversed < 2; ++cigar_reversed)
        {
            Cigar cigar(cigar_string, cigar_reversed, rf_start, qr_start);
            Read_Chunk chunk;
            shared_ptr< Mutation_Cont > mut_cont_sptr;
            std::tie(chunk, mut_cont_sptr) = Read_Chunk::make_chunk_from_cigar(cigar);

            cout << chunk << nl;
            cout << "Mutation_Cont:" << inc_tab;
            print_seq(cout, *mut_cont_sptr, nl, nl);
            cout << dec_tab << nl << "tests:" << inc_tab;

            test_traversal(chunk, 0, true);

            Size_Type brk;
            for (brk = cigar.get_rf_start(); brk <= cigar.get_rf_start() + cigar.get_rf_len(); ++brk)
            {
                test_traversal(chunk, brk, true);
            }
            for (brk = cigar.get_qr_start(); brk <= cigar.get_qr_start() + cigar.get_qr_len(); ++brk)
            {
                test_traversal(chunk, brk, false);
            }
            cout << dec_tab << nl << "done" << nl;
        }
    }

    return 0;
}
