#ifdef NDEBUG
#undef NDEBUG
#endif

#include <string>
#include <iostream>
#include <tclap/CmdLine.h>
#include "RC_Sequence.hpp"
#include "Util.h"

using namespace std;

namespace global
{
    string prog_desc =
        "Test RC Sequence type. "
        "The program repeats the following random test several times:\n"
        " - pick a DNA sequence\n"
        " - pick a substring position\n"
        " - optionally do a reverse operation\n"
        " - optionally do a complement operation\n"
        " - pick order of substring vs rev&comp\n"
        " - perform the operations using the old library\n"
        " - perform the operations using the new library\n\n";
    TCLAP::CmdLine cmd_parser(prog_desc);

    TCLAP::ValueArg< long > seed("", "seed", "Random seed (-1: use time).", false, -1, "int", cmd_parser);
    TCLAP::ValueArg< unsigned > n_ops("", "n-ops", "Number of operations.", false, 100000, "int", cmd_parser);
    TCLAP::ValueArg< unsigned > max_len("", "max-len", "Maximum string length.", false, 100, "int", cmd_parser);
    TCLAP::SwitchArg rc_only("", "rc-only", "Only id or rc operations.", cmd_parser, false);
} // namespace global

typedef rc_sequence::Sequence< std::string > sequence_type;
typedef sequence_type::proxy_type sequence_proxy_type;

size_t rand_int(size_t max)
{
    return static_cast< size_t >(drand48() * max);
}

void real_main()
{
    cout << "parameters: seed=" << global::seed << " max_len=" << global::max_len
         << " n_ops=" << global::n_ops << " rc_only=" << global::rc_only << "\n";
    size_t old_lib_us = 0;
    size_t new_lib_us = 0;
    vector< size_t > op_count(8, 0);
    const vector< string > op_tag = { "id + substr", "comp + substr", "rev + substr", "revcomp + substr",
                                      "substr + id", "substr + comp", "substr + rev", "substr + revcomp" };
    for (size_t i = 0; i < global::n_ops; ++i)
    {
        // pick 2 strings
        size_t s_len = rand_int(global::max_len - 1) + 1;
        size_t s_comp_len = rand_int(global::max_len - 1) + 1;
        string s_old_lib;
        string s_comp_old_lib;
        for (size_t j = 0; j < s_len; ++j)
        {
            static const string char_s = "ACGT";
            size_t char_idx = rand_int(4);
            s_old_lib += char_s[char_idx];
        }
        for (size_t j = 0; j < s_comp_len; ++j)
        {
            static const string char_s = "ACGT";
            size_t char_idx = rand_int(4);
            s_comp_old_lib += char_s[char_idx];
        }
        sequence_type s_new_lib(s_old_lib);
        sequence_type s_comp_new_lib(s_comp_old_lib);
        // pick operation order, rev/comp combo, and substr coordinates
        bool substr_first = rand_int(2);
        bool do_rev = rand_int(2);
        bool do_comp = global::rc_only? do_rev : rand_int(2);
        size_t cut_start = rand_int(s_len);
        size_t cut_len = rand_int(s_len - cut_start) + 1;
        size_t check_idx = rand_int(cut_len);
        ++op_count[ (static_cast< int >(substr_first) << 2) +
                    (static_cast< int >(do_rev) << 1) +
                    static_cast< int >(do_comp) ];
        // print operation
        clog << "s='" << s_old_lib << "' substr_first=" << substr_first << " do_rev=" << do_rev
             << " do_comp=" << do_comp << " cut_start=" << cut_start << " cut_len=" << cut_len
             << " check_idx=" << check_idx << " s_comp='" << s_comp_old_lib << "'\n";
        // apply operations using old library
        string res_old_lib;
        int res_comp_old_lib;
        {
            std::chrono::steady_clock::time_point start_t = std::chrono::steady_clock::now();
            if (substr_first and do_rev and do_comp)
            {
                res_old_lib = reverseComplement(s_old_lib.substr(cut_start, cut_len));
            }
            else if (substr_first and do_rev and not do_comp)
            {
                res_old_lib = reverse(s_old_lib.substr(cut_start, cut_len));
            }
            else if (substr_first and not do_rev and do_comp)
            {
                res_old_lib = complement(s_old_lib.substr(cut_start, cut_len));
            }
            else if (substr_first and not do_rev and not do_comp)
            {
                res_old_lib = s_old_lib.substr(cut_start, cut_len);
            }
            else if (not substr_first and do_rev and do_comp)
            {
                res_old_lib = reverseComplement(s_old_lib).substr(cut_start, cut_len);
            }
            else if (not substr_first and do_rev and not do_comp)
            {
                res_old_lib = reverse(s_old_lib).substr(cut_start, cut_len);
            }
            else if (not substr_first and not do_rev and do_comp)
            {
                res_old_lib = complement(s_old_lib).substr(cut_start, cut_len);
            }
            else if (not substr_first and not do_rev and not do_comp)
            {
                res_old_lib = s_old_lib.substr(cut_start, cut_len);
            }
            res_comp_old_lib = res_old_lib.compare(s_comp_old_lib);
            std::chrono::steady_clock::time_point end_t = std::chrono::steady_clock::now();
            old_lib_us += std::chrono::duration_cast<std::chrono::microseconds>(end_t - start_t).count();
        }
        // apply operations using new library
        sequence_proxy_type res_new_lib;
        int res_comp_new_lib;
        {
            std::chrono::steady_clock::time_point start_t = std::chrono::steady_clock::now();
            if (substr_first)
            {
                res_new_lib = std::move(s_new_lib.substr(cut_start, cut_len).rev(do_rev).comp(do_comp));
            }
            else
            {
                res_new_lib = std::move(s_new_lib.rev(do_rev).comp(do_comp).substr(cut_start, cut_len));
            }
            res_comp_new_lib = res_new_lib.compare(s_comp_new_lib);
            std::chrono::steady_clock::time_point end_t = std::chrono::steady_clock::now();
            new_lib_us += std::chrono::duration_cast<std::chrono::microseconds>(end_t - start_t).count();
        }
        assert(res_new_lib == res_old_lib);
        assert(res_new_lib[check_idx] == res_old_lib[check_idx]);
        assert((res_comp_new_lib < 0) == (res_comp_old_lib < 0));
        assert((res_comp_new_lib > 0) == (res_comp_old_lib > 0));
    }
    cout << "op counts:\n";
    for (size_t i = 0; i < 8; ++i)
    {
        cout << "  " << op_tag[i] << ": " << op_count[i] << "\n";
    }
    cout << "old lib time: " << old_lib_us << "\n";
    cout << "new lib time: " << new_lib_us << "\n";
}

int main(int argc, char * argv[])
{
    global::cmd_parser.parse(argc, argv);
    if (global::seed < 0)
    {
        global::seed.get() = time(nullptr);
    }
    srand48(global::seed);
    real_main();
}
