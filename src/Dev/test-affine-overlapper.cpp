#include <iostream>
#include <time.h>
#include <tclap/CmdLine.h>

#include "overlapper.h"
#include "Cigar.hpp"

using namespace std;

namespace global
{
    string prog_desc =
        "Test affine overlapper with random queries.";

    TCLAP::CmdLine cmd_parser(prog_desc);

    TCLAP::ValueArg< long > seed("s", "seed", "Random seed (-1: use time).", false, -1, "int", cmd_parser);
    TCLAP::ValueArg< unsigned > n_ops("n", "n-ops", "Number of operations.", false, 1000, "int", cmd_parser);
    TCLAP::SwitchArg check_alt_custom("c", "check-alt-custom", "Check ALT_CUSTOM alignment.", cmd_parser, false);
    TCLAP::SwitchArg use_random_scores("r", "use-random-scores", "Use random scores instead of defaults.", cmd_parser, false);
} // namespace global

typedef cigar::Cigar<> Cigar;

void real_main()
{
    cout << "using seed=" << global::seed << "\n";
    for (size_t i = 0; i < global::n_ops; ++i)
    {
        // construct random strings
        int s_len[2];
        string s[2];
        for (int k = 0; k < 2; ++k)
        {
            s_len[k] = 1 + static_cast< size_t >(drand48() * 19);
            for (int j = 0; j < s_len[k]; ++j)
            {
                static const string bases = "ACGT";
                s[k].push_back(bases[static_cast< size_t >(drand48() * 4)]);
            }
        }
        // pick random alignment type
        static const vector< AffineAlignmentType > alt_types = { ALT_OVERLAP, ALT_GLOBAL, ALT_CONTAINMENT, ALT_CUSTOM };
        static const vector< string > alt_names = { "ALT_OVERLAP", "ALT_GLOBAL", "ALT_CONTAINMENT", "ALT_CUSTOM" };
        OverlapperParams params = affine_default_params;
        size_t alt_type_idx = static_cast< size_t >(drand48() * (3 + static_cast< int >(global::check_alt_custom)));
        params.type = static_cast< AffineAlignmentType >(alt_type_idx);
        if (params.type == ALT_CUSTOM)
        {
            params.gap_s1_start = static_cast< int >(drand48() * 2);
            params.gap_s1_end = static_cast< int >(drand48() * 2);
            params.gap_s2_start = static_cast< int >(drand48() * 2);
            params.gap_s2_end = static_cast< int >(drand48() * 2);
        }
        if (global::use_random_scores)
        {
            params.match_score = static_cast< int >(drand48() * 20);
            params.mismatch_penalty = - static_cast< int >(drand48() * 20);
            params.gap_penalty = - static_cast< int >(drand48() * 20);
            params.gap_ext_penalty = params.gap_penalty + static_cast< int >(drand48() * (- params.gap_penalty));
        }
        // print settings
        cout << "checking alignment: s[0]='" << s[0] << "' s[1]='" << s[1] << "' type="
                  << alt_names[alt_type_idx];
        if (params.type == ALT_CUSTOM)
        {
            cout << "("
                 << params.gap_s1_start << ","
                 << params.gap_s1_end << ","
                 << params.gap_s2_start << ","
                 << params.gap_s2_end << ")";
        }
        if (global::use_random_scores)
        {
            cout << " scores=("
                 << params.match_score << ","
                 << params.mismatch_penalty << ","
                 << params.gap_penalty << ","
                 << params.gap_ext_penalty << ")";
        }
        cout << "\n";

        // construct alignment
        auto output = Overlapper::computeAlignmentAffine(s[0], s[1], params);
        cout << "cigar='" << output.cigar << "' score=" << output.score << "\n";

        // check start&end positions
        if (params.type == ALT_OVERLAP)
        {
            assert(output.match[0].start == 0 or output.match[1].start == 0);
            assert(output.match[0].end == s_len[0] - 1 or output.match[1].end == s_len[1] - 1);
        }
        else if (params.type == ALT_GLOBAL)
        {
            assert(output.match[0].start == 0);
            assert(output.match[0].end == s_len[0] - 1);
            assert(output.match[1].start == 0);
            assert(output.match[1].end == s_len[1] - 1);
        }
        else if (params.type == ALT_CONTAINMENT)
        {
            assert(output.match[1].start == 0);
            assert(output.match[1].end == s_len[1] - 1);
        }
        else if (params.type == ALT_CUSTOM)
        {
            assert(params.gap_s1_start or output.match[0].start == 0);
            assert(params.gap_s1_end or output.match[0].end == s_len[0] - 1);
            assert(params.gap_s2_start or output.match[1].start == 0);
            assert(params.gap_s2_end or output.match[1].end == s_len[1] - 1);
        }
        else
        {
            abort();
        }

        // check cigar length
        Cigar cigar(output.cigar);
        assert(static_cast< long >(cigar.rf_len()) == output.match[0].end - output.match[0].start + 1);
        assert(static_cast< long >(cigar.qr_len()) == output.match[1].end - output.match[1].start + 1);

        // check scores
        cigar.disambiguate(s[0].substr(output.match[0].start, output.match[0].end - output.match[0].start + 1),
                           s[1].substr(output.match[1].start, output.match[1].end - output.match[1].start + 1));
        int score = 0;
        for (size_t j = 0; j < cigar.n_ops(); ++j)
        {
            if (cigar.op_type(j) == '=')
            {
                score += cigar.op_len(j) * params.match_score;
            }
            else if (cigar.op_type(j) == 'X')
            {
                score += cigar.op_len(j) * params.mismatch_penalty;
            }
            else if (cigar.op_type(j) == 'I' || cigar.op_type(j) == 'D')
            {
                score += params.gap_penalty + cigar.op_len(j) * params.gap_ext_penalty;
            }
        }
        assert(score == output.score);
    }
    cout << "all good\n";
}

int main(int argc, char * argv[])
{
    global::cmd_parser.parse(argc, argv);
    global_assert::prog_name() = global::cmd_parser.getProgramName();
    if (global::seed < 0)
    {
        global::seed.get() = time(nullptr);
    }
    srand48(global::seed);
    real_main();
}
