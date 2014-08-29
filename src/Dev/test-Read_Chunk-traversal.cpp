#include <iostream>
#include <string>
#include <tuple>
#include "Read_Chunk.hpp"
#include "Contig_Entry.hpp"
#include "Read_Entry.hpp"

using namespace std;
using namespace MAC;

typedef std::list< std::vector< size_t > > Pos_List;

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

void test_traversal(Read_Chunk_BPtr rc_bptr, bool forward, Size_Type brk, bool on_contig,
                    const Cigar& cigar, const Pos_List& res)
{
    ostringstream oss;
    oss << res;
    auto tmp_ptree = ptree()
        .put("cigar", cigar.to_ptree())
        .put("forward", forward)
        .put("brk", brk)
        .put("on_contig", on_contig)
        .put("res", oss.str());

    clog << "checking: cigar=" << cigar.to_string() << " same_dir=" << not cigar.is_reversed()
         << " rf_start=" << cigar.get_rf_start() << " qr_start=" << cigar.get_qr_start()
         << " forward=" << forward << " brk=" << brk << " on_contig=" << on_contig
         << " res=" << res << "\n";

    Read_Chunk_Pos pos = forward? rc_bptr->get_start_pos() : rc_bptr->get_end_pos();
    Read_Chunk_Pos end_pos = forward? rc_bptr->get_end_pos() : rc_bptr->get_start_pos();
    auto it = res.begin();
    auto mca_cit = forward? rc_bptr->mut_ptr_cont().begin() : rc_bptr->mut_ptr_cont().end();
    bool bp_reached = false;
    size_t i = 0;
    while (true)
    {
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
        //clog << "checked position: " << i << "\n";
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

Pos_List add_brk_to_fwd_pos_list(const Pos_List& pos_list, Size_Type brk, bool on_contig)
{
    // must be forward traversal
    assert(pos_list.front()[0] <= pos_list.back()[0]);
    // brk==0 means no breakpoint
    if (brk == 0)
    {
        return pos_list;
    }
    // brk must not be reached already
    assert(not on_contig or brk > pos_list.front()[0]);
    bool same_dir = pos_list.front()[1] <= pos_list.back()[1];
    assert(on_contig or (same_dir? brk > pos_list.front()[1] : brk < pos_list.front()[1]));
    // construct new list
    Pos_List res;
    res.emplace_back(pos_list.front());
    for (auto it = ++pos_list.begin(); it != pos_list.end(); ++it)
    {
        // lengths of the current mapped stretch
        Size_Type c_len = (*it)[0] - res.back()[0];
        Size_Type r_len = same_dir? (*it)[1] - res.back()[1] : res.back()[1] - (*it)[1];
        if (on_contig and (*it)[0] >= brk)
        {
            // just passed breakpoint on contig
            if ((*it)[0] > brk)
            {
                assert(res.back()[0] < brk and brk < (*it)[0]);
                size_t c_delta = brk - res.back()[0];
                size_t r_delta = std::min(r_len, c_delta);
                res.emplace_back(std::vector< size_t >({
                    brk,
                    same_dir? res.back()[1] + r_delta : res.back()[1] - r_delta,
                    0,
                    (*it)[2] > 0? c_delta : 0
                }));
            }
            brk = pos_list.back()[0] + 1; // unreachable
        }
        else if (not on_contig and same_dir and (*it)[1] >= brk)
        {
            // just passed breakpoint on read, going in same dir
            if ((*it)[1] > brk)
            {
                assert(res.back()[1] < brk and brk < (*it)[1]);
                size_t r_delta = brk - res.back()[1];
                size_t c_delta = std::min(c_len, r_delta);
                res.emplace_back(std::vector< size_t >({
                    res.back()[0] + c_delta,
                    brk,
                    0,
                    (*it)[2] > 0? r_delta : 0
                }));
            }
            brk = pos_list.back()[1] + 1; // unreachable
        }
        else if (not on_contig and not same_dir and brk >= (*it)[1])
        {
            // just passed breakpoint on read, going in opposite dir
            if (brk > (*it)[1])
            {
                assert((*it)[1] < brk and brk < res.back()[1]);
                size_t r_delta = res.back()[1] - brk;
                size_t c_delta = std::min(c_len, r_delta);
                res.emplace_back(std::vector< size_t >({
                    res.back()[0] + c_delta,
                    brk,
                    0,
                    (*it)[2] > 0? r_delta : 0
                }));
            }
            brk = 0; // unreachable
        }
        res.emplace_back(*it);
    }
    return res;
}



int main()
{
    Mutation_Fact mut_fact;
    Mutation_Chunk_Adapter_Fact mca_fact;
    Read_Chunk_Fact rc_fact;
    Contig_Entry_Fact ce_fact;
    Read_Entry_Fact re_fact;

    //string cigar_string;

    //while (cin >> cigar_string)
    for (auto cigar_string : std::list< std::string >({"1=", "1D", "1I", "1D1=", "1I1=", "3=2D3=2I3="}))
    {
        for (int cigar_reversed = 0; cigar_reversed < 2; ++cigar_reversed)
        {
            for (const auto& p : std::vector< std::pair< Size_Type, Size_Type > >({
                std::make_pair(0, 0),
                std::make_pair(0, 10),
                std::make_pair(10, 0),
                std::make_pair(10, 10)
            }))
            {
                Size_Type rf_start;
                Size_Type qr_start;
                std::tie(rf_start, qr_start) = p;
                Cigar cigar(cigar_string, cigar_reversed, rf_start, qr_start);
                Seq_Type rf_seq; for (size_t i = 0; i < cigar.get_rf_len(); ++i) rf_seq += "A";
                Seq_Type qr_seq; for (size_t i = 0; i < cigar.get_qr_len(); ++i) qr_seq += "A";
                Read_Chunk_BPtr chunk_bptr;
                Contig_Entry_BPtr ce_bptr;
                chunk_bptr = Read_Chunk::make_chunk_from_cigar(cigar, Seq_Type(rf_seq), qr_seq);
                ce_bptr = chunk_bptr->ce_bptr();

                /*
                clog << ptree("test")
                    .put("chunk", chunk_bptr->to_ptree())
                    .put("mut_cont", cont_to_ptree(ce_bptr->mut_cont()));
                */

                // check traversal with no breakpoints
                Pos_List fwd_pos_list = get_expected_pos_list(cigar);
                test_traversal(chunk_bptr, true, 0, true, cigar, fwd_pos_list);
                auto rev_pos_list = reverse_pos_list(fwd_pos_list);
                test_traversal(chunk_bptr, false, 0, true, cigar, rev_pos_list);

                // check breakpoints on contig
                for (Size_Type c_brk = cigar.get_rf_start() + 1; c_brk < cigar.get_rf_end(); ++c_brk)
                {
                    auto pos_list = add_brk_to_fwd_pos_list(fwd_pos_list, c_brk, true);
                    test_traversal(chunk_bptr, true, c_brk, true, cigar, pos_list);
                    rev_pos_list = reverse_pos_list(pos_list);
                    test_traversal(chunk_bptr, false, c_brk, true, cigar, rev_pos_list);
                }

                // check breakpoints on read
                for (Size_Type r_brk = not cigar_reversed? cigar.get_qr_start() + 1 : (cigar.get_qr_end() > 0? cigar.get_qr_end() - 1 : 0);
                     not cigar_reversed? r_brk < cigar.get_qr_end() : r_brk > cigar.get_qr_start();
                     not cigar_reversed? ++r_brk : --r_brk)
                {
                    auto pos_list = add_brk_to_fwd_pos_list(fwd_pos_list, r_brk, false);
                    test_traversal(chunk_bptr, true, r_brk, false, cigar, pos_list);
                    rev_pos_list = reverse_pos_list(pos_list);
                    test_traversal(chunk_bptr, false, r_brk, false, cigar, rev_pos_list);
                }
            }
        }
    }
    clog << "all good\n";
    return 0;
}
