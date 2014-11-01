#include <string>
#include <iostream>
#include <boost/program_options.hpp>
#include "DNA_Sequence.hpp"
#include "Util.h"

using namespace std;
namespace bo = boost::program_options;

typedef dnasequence::Sequence< std::string > sequence_type;
typedef sequence_type::proxy_type sequence_proxy_type;

struct Program_Options
{
    size_t max_len;
    size_t n_ops;
    size_t seed;
    bool rc_only;
};

size_t rand_int(size_t max)
{
    return static_cast< size_t >(drand48() * max);
}

void real_main(const Program_Options& po)
{
    cout << "parameters: seed=" << po.seed << " max_len=" << po.max_len
         << " n_ops=" << po.n_ops << " rc_only=" << po.rc_only << "\n";
    size_t old_lib_us = 0;
    size_t new_lib_us = 0;
    vector< size_t > op_count(8, 0);
    const vector< string > op_tag = { "id + substr", "comp + substr", "rev + substr", "revcomp + substr",
                                      "substr + id", "substr + comp", "substr + rev", "substr + revcomp" };
    for (size_t i = 0; i < po.n_ops; ++i)
    {
        // pick 2 strings
        size_t s_len = rand_int(po.max_len - 1) + 1;
        size_t s_comp_len = rand_int(po.max_len - 1) + 1;
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
        bool do_comp = po.rc_only? do_rev : rand_int(2);
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

int main(int argc, char* argv[])
{
    Program_Options po;

    bo::options_description generic_opts_desc("Generic options");
    bo::options_description config_opts_desc("Configuration options");
    bo::options_description hidden_opts_desc("Hidden options");
    bo::options_description cmdline_opts_desc;
    bo::options_description visible_opts_desc("Allowed options");
    generic_opts_desc.add_options()
        ("help,?", "produce help message")
        ;
    config_opts_desc.add_options()
        ("max-len", bo::value<size_t>(&po.max_len)->default_value(100), "maximum string length")
        ("n-ops", bo::value<size_t>(&po.n_ops)->default_value(100000), "number of operations")
        ("seed", bo::value<size_t>(&po.seed)->default_value(0), "random number generator seed")
        ("rc-only", bo::bool_switch(&po.rc_only), "only id or rc operations")
        ;
    cmdline_opts_desc.add(generic_opts_desc).add(config_opts_desc).add(hidden_opts_desc);
    visible_opts_desc.add(generic_opts_desc).add(config_opts_desc);
    bo::variables_map vm;
    store(bo::command_line_parser(argc, argv).options(cmdline_opts_desc).run(), vm);
    notify(vm);
    if (vm.count("help"))
    {
        cout << "Test DNA Sequence type.\n"
            "The program repeats the following random test several times:\n"
            " - pick a DNA sequence\n"
            " - pick a substring position\n"
            " - optionally do a reverse operation\n"
            " - optionally do a complement operation\n"
            " - pick order of substring vs rev&comp\n"
            " - perform the operations using the old library\n"
            " - perform the operations using the new library\n\n";
        cout << visible_opts_desc;
        exit(EXIT_SUCCESS);
    }
    if (po.seed == 0)
    {
        po.seed = time(NULL);
    }
    real_main(po);
}
