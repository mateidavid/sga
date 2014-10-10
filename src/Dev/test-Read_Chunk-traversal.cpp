#include <iostream>
#include <string>
#include <tuple>
#include "Read_Chunk.hpp"
#include "Contig_Entry.hpp"
#include "Read_Entry.hpp"
#include "fstr.hpp"

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

    clog << "checking: cigar=" << cigar.to_string() << " same_dir=" << not cigar.reversed()
         << " rf_start=" << cigar.rf_start() << " qr_start=" << cigar.qr_start()
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
        global::assert_message() = oss.str();
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

Pos_List get_expected_pos_list(Read_Chunk_CBPtr rc_cbptr)
{
    Pos_List res;
    bool same_dir = not rc_cbptr->get_rc();
    Size_Type c_pos = rc_cbptr->get_c_start();
    Size_Type r_pos = same_dir? rc_cbptr->get_r_start() : rc_cbptr->get_r_end();
    if (rc_cbptr->mut_ptr_cont().empty()
        or rc_cbptr->mut_ptr_cont().begin()->mut_cbptr()->get_start() > rc_cbptr->get_c_start())
    {
        res.emplace_back(std::vector< size_t >({ c_pos, r_pos, 0, 0 }));
    }
    for (auto mca_cbref : rc_cbptr->mut_ptr_cont())
    {
        Mutation_CBPtr mut_cbptr = (&mca_cbref)->mut_cbptr();
        assert(mut_cbptr->get_start() >= c_pos);
        Size_Type skip_len = mut_cbptr->get_start() - c_pos;
        assert(not res.empty() or skip_len == 0);
        c_pos += skip_len;
        r_pos = same_dir? r_pos + skip_len : r_pos - skip_len;
        res.emplace_back(std::vector< size_t >({ c_pos, r_pos, 0, 0 }));
        c_pos += mut_cbptr->get_len();
        r_pos = same_dir? r_pos + mut_cbptr->get_seq_len() : r_pos - mut_cbptr->get_seq_len();
        res.emplace_back(std::vector< size_t >({ c_pos, r_pos, 1, 0 }));
    }
    if (rc_cbptr->mut_ptr_cont().empty()
        or rc_cbptr->mut_ptr_cont().rbegin()->mut_cbptr()->get_end() < rc_cbptr->get_c_end())
    {
        Size_Type skip_len = rc_cbptr->get_c_end() - c_pos;
        assert(skip_len > 0);
        c_pos += skip_len;
        r_pos = same_dir? r_pos + skip_len : r_pos - skip_len;
        res.emplace_back(std::vector< size_t >({ c_pos, r_pos, 0, 0 }));
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

void test_cigar(const std::string& cigar_string,
                const std::map< std::string, std::vector< size_t > >& init_settings_map = std::map< std::string, std::vector< size_t > >())
{
    std::map< std::string, std::vector< size_t > > settings_map(init_settings_map);
    if (settings_map.count("cigar_orientation_v") == 0)
    {
        settings_map.emplace("cigar_orientation_v", std::vector< size_t >({ 0, 1 }) );
    }
    if (settings_map.count("rf_start_v") == 0)
    {
        settings_map.emplace("rf_start_v", std::vector< size_t >({ 0, 10 }) );
    }
    if (settings_map.count("qr_start_v") == 0)
    {
        settings_map.emplace("qr_start_v", std::vector< size_t >({ 0, 10 }) );
    }

    for (auto cigar_orientation : settings_map["cigar_orientation_v"])
    for (auto rf_start : settings_map["rf_start_v"])
    for (auto qr_start : settings_map["qr_start_v"])
    {
        // construct cigar
        Cigar cigar(cigar_string, cigar_orientation, rf_start, qr_start);
        Seq_Type rf_seq; for (size_t i = 0; i < cigar.rf_len(); ++i) rf_seq += "A";
        Seq_Type qr_seq; for (size_t i = 0; i < cigar.qr_len(); ++i) qr_seq += "A";
        // construct graph structures
        Read_Chunk_BPtr chunk_bptr;
        chunk_bptr = Read_Chunk::make_chunk_from_cigar(cigar, Seq_Type(rf_seq), qr_seq);
        // compute original (no-breakpoint, forward direction) position list
        Pos_List fwd_pos_list = get_expected_pos_list(chunk_bptr);
        // compute breakpoint positions if not specified
        if (settings_map.count("rel_c_brk_v") == 0)
        {
            settings_map.emplace("rel_c_brk_v", std::vector< size_t >({ 0 }) );
            for (Size_Type c_brk = 1; c_brk < cigar.rf_len(); ++c_brk)
            {
                settings_map["rel_c_brk_v"].push_back(c_brk);
            }
        }
        for (auto rel_c_brk : settings_map["rel_c_brk_v"])
        {
            Size_Type c_brk = rel_c_brk > 0? rf_start + rel_c_brk : 0;
            Pos_List pos_list = add_brk_to_fwd_pos_list(fwd_pos_list, c_brk, true);
            test_traversal(chunk_bptr, true, c_brk, true, cigar, pos_list);
            Pos_List rev_pos_list = reverse_pos_list(pos_list);
            test_traversal(chunk_bptr, false, c_brk, true, cigar, rev_pos_list);
        }
        if (settings_map.count("rel_r_brk_v") == 0)
        {
            for (Size_Type r_brk = 1; r_brk < cigar.qr_len(); ++r_brk)
            {
                settings_map["rel_r_brk_v"].push_back(r_brk);
            }
        }
        for (auto rel_r_brk : settings_map["rel_r_brk_v"])
        {
            Size_Type r_brk = rel_r_brk > 0? qr_start + rel_r_brk : 0;
            Pos_List pos_list = add_brk_to_fwd_pos_list(fwd_pos_list, r_brk, false);
            test_traversal(chunk_bptr, true, r_brk, false, cigar, pos_list);
            Pos_List rev_pos_list = reverse_pos_list(pos_list);
            test_traversal(chunk_bptr, false, r_brk, false, cigar, rev_pos_list);
        }
        // destroy graph structures
        Contig_Entry::dispose(chunk_bptr->ce_bptr());
    }
}

void usage(std::ostream& os, const std::string& prog_name)
{
    os << "use: " << prog_name << " [ <file> ]\n\n"
       << "For each cigar string, the program constructs a Read_Chunk object and checks its traversal using\n"
       << "Read_Chunk_Pos::increment() and decrement() methods. For each cigar the following combinations\n"
       << "of parameters are tested:\n"
       << "- same or reversed cigar orientation\n"
       << "- rf_start in {0, 10}\n"
       << "- qr_start in {0, 10}\n"
       << "- all possible contig breakpoints, including no breakpoint\n"
       << "- all possible read breakpoints\n"
       << "- foward or backward traversal\n";
}

int main(int argc, char* argv[])
{
    global::program_name() = argv[0];

    Mutation_Fact mut_fact;
    Mutation_Chunk_Adapter_Fact mca_fact;
    Read_Chunk_Fact rc_fact;
    Contig_Entry_Fact ce_fact;
    Read_Entry_Fact re_fact;

    if (argc > 2)
    {
        usage(std::cerr, argv[0]);
        abort();
    }
    if (argc == 2)
    {
        string filename(argv[1]);
        if (filename == "-?" or filename == "--help")
        {
            usage(std::cout, argv[0]);
            exit(EXIT_SUCCESS);
        }
        fstr tmp_fs(filename);
        string cigar_string;
        while (tmp_fs >> cigar_string)
        {
            test_cigar(cigar_string);
        }
    }
    else
    {
        for (auto cigar_string : std::list< std::string >({
            "1=", "3=", "1D", "3D", "1I", "3I",
            "1D1=", "3D1=", "1D3=", "3D3=",
            "1I1=", "3I1=", "1I3=", "3I3=",
            "1D1I", "3D1I", "1D3I", "3D3I",
            "1=1X", "1=3X", "3=1X", "3=3X",
            "1=1X1D1I1X1D1I1D1I1=1=1X",
            "3=2D3=2I3="}))
        {
            test_cigar(cigar_string);
        }
    }
    clog << "all good\n";
}
