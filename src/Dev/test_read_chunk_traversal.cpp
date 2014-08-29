#include <iostream>
#include <string>
#include <tuple>
#include "Read_Chunk.hpp"
#include "Contig_Entry.hpp"
#include "Read_Entry.hpp"

using namespace std;
using namespace MAC;

typedef std::list< std::vector< size_t > > Pos_List;

void __test_traversal(Read_Chunk_BPtr chunk_bptr, Size_Type brk, bool on_contig)
{
    cout << "testing brk=" << brk << ",on_contig=" << on_contig << "\n";

    list< Read_Chunk_Pos > l;
    Read_Chunk_Pos pos = chunk_bptr->get_start_pos();
    bool reached_breakpoint = false;
    bool forward = true;
    while (true)
    {
        if (forward and pos == chunk_bptr->get_start_pos())
        {
            cout << "at start, going forward\n";
        }

        if (not forward and pos == chunk_bptr->get_end_pos())
        {
            cout << "at end, going backward\n";
        }

        cout << "pos: " << pos.to_ptree();

        if (not reached_breakpoint
            and ((on_contig and pos.c_pos == brk) or (not on_contig and pos.r_pos == brk)))
        {
            cout << "breakpoint reached\n";
            reached_breakpoint = true;
        }

        if (forward)
        {
            l.push_back(pos);
        }
        else
        {
            if (not (pos == l.back()))
            {
                cerr << "mismatched positions:\n"
                     << "forward: " << l.back().to_ptree()
                     << "backward: " << pos.to_ptree();
                abort();
            }
            l.erase(--l.end());
        }

        if (not forward and pos == chunk_bptr->get_start_pos())
        {
            cout << "at start, going backward: done\n";
            break;
        }

        if (forward and pos == chunk_bptr->get_end_pos())
        {
            cout << "at end, forward: switching direction\n";
            forward = false;
            reached_breakpoint = false;
            continue;
        }

        if (forward)
            pos.increment(reached_breakpoint? 0 : brk, on_contig);
        else
            pos.decrement(reached_breakpoint? 0 : brk, on_contig);
    }
    if (l.size() > 0)
    {
        cerr << "position vector not empty:\n";
        for (const auto& e : l) cerr << e.to_ptree();
        abort();
    }
}

void test_traversal(Read_Chunk_BPtr rc_bptr, bool forward, Size_Type brk, bool on_contig,
                    const string& cigar, const Pos_List& res)
{
    auto tmp_ptree = ptree()
        .put("cigar", cigar)
        .put("forward", forward)
        .put("brk", brk)
        .put("on_contig", on_contig);
    Read_Chunk_Pos pos = forward? rc_bptr->get_start_pos() : rc_bptr->get_end_pos();
    Read_Chunk_Pos end_pos = forward? rc_bptr->get_end_pos() : rc_bptr->get_start_pos();
    auto it = res.begin();
    auto mca_cit = forward? rc_bptr->mut_ptr_cont().begin() : rc_bptr->mut_ptr_cont().end();
    bool bp_reached = false;
    size_t i = 0;
    while (true)
    {
        ostringstream oss;
        oss.str(string());
        oss << tmp_ptree.put("i", i);
        global::assert_message = oss.str();
        ASSERT(it != res.end());
        ASSERT(pos.c_pos == (*it)[0]);
        ASSERT(pos.r_pos == (*it)[1]);
        if ((*it)[2])
        {
            forward? ++mca_cit : --mca_cit;
        }
        ASSERT(pos.mca_cit == mca_cit);
        ASSERT(pos.mut_offset == (*it)[3]);
        clog << "checked position: " << i << "\n";
        if (pos == end_pos)
        {
            ASSERT(++it == res.end());
            break;
        }
        pos.advance(forward, not bp_reached? brk : 0, on_contig);
        bp_reached |= (on_contig? pos.c_pos == brk : pos.r_pos == brk);
        ++it;
        ++i;
    }
}

Pos_List get_expected_pos_list(const Cigar& cigar)
{
    Pos_List res;
    for (size_t i = 0; i <= cigar.get_n_ops(); ++i)
    {
        res.emplace_back(std::vector< size_t >({
            cigar.get_rf_offset(i),
            cigar.get_qr_offset(i),
            static_cast< size_t >(i > 0 and not cigar.is_match(i - 1)? 1 : 0),
            0}));
    }
    return res;
}

Pos_List reverse_pos_list(const Pos_List& pos_list)
{
    Pos_List res;
    size_t last_it_move = 0;
    for (auto it = pos_list.rbegin(); it != pos_list.rend(); ++it)
    {
        res.emplace_back(*it);
        std::swap(last_it_move, res.back()[2]);
    }
    return res;
}

std::ostream& operator << (std::ostream& os, const Pos_List& pos_list)
{
    os << "{";
    for (auto it = pos_list.begin(); it != pos_list.end(); ++it)
    {
        if (it != pos_list.begin())
        {
            os << ",";
        }
        os << "(" << (*it)[0] << "," << (*it)[1] << "," << (*it)[2] << "," << (*it)[3] << ")";
    }
    os << "}";
    return os;
}

int main()
{
    Mutation_Fact mut_fact;
    Mutation_Chunk_Adapter_Fact mca_fact;
    Read_Chunk_Fact rc_fact;
    Contig_Entry_Fact ce_fact;
    Read_Entry_Fact re_fact;

    string cigar_string;
    Size_Type rf_start = 0;
    Size_Type qr_start = 0;

    //while (cin >> cigar_string >> rf_start >> qr_start)
    {
        cigar_string = "3=2D3=2I3=";
        for (int cigar_reversed = 0; cigar_reversed < 2; ++cigar_reversed)
        {
            //int cigar_reversed = 0;
            Cigar cigar(cigar_string, cigar_reversed, rf_start, qr_start);
            Seq_Type rf_seq; for (size_t i = 0; i < cigar.get_rf_len(); ++i) rf_seq += "A";
            Seq_Type qr_seq; for (size_t i = 0; i < cigar.get_qr_len(); ++i) qr_seq += "A";
            Read_Chunk_BPtr chunk_bptr;
            Contig_Entry_BPtr ce_bptr;
            chunk_bptr = Read_Chunk::make_chunk_from_cigar(cigar, Seq_Type(rf_seq), qr_seq);
            ce_bptr = chunk_bptr->ce_bptr();

            clog << ptree("test")
                .put("chunk", chunk_bptr->to_ptree())
                .put("mut_cont", cont_to_ptree(ce_bptr->mut_cont()));

            Pos_List pos_list = get_expected_pos_list(cigar);
            cout << "pos_list: " << pos_list << "\n";
            auto rev_pos_list = reverse_pos_list(pos_list);
            cout << "rev_pos_list: " << rev_pos_list << "\n";

            test_traversal(chunk_bptr, true, 0, true, cigar_string, pos_list);
            test_traversal(chunk_bptr, false, 0, true, cigar_string, rev_pos_list);

            /*
            Size_Type brk;
            for (brk = cigar.get_rf_start(); brk <= cigar.get_rf_start() + cigar.get_rf_len(); ++brk)
            {
                test_traversal(chunk_bptr, brk, true);
            }
            for (brk = cigar.get_qr_start(); brk <= cigar.get_qr_start() + cigar.get_qr_len(); ++brk)
            {
                test_traversal(chunk_bptr, brk, false);
            }
            */
        }
    }

    return 0;
}
