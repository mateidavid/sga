#include <iostream>
#include <time.h>
#include <boost/program_options.hpp>

#include "overlapper.h"
#include "Cigar.hpp"

namespace bo = boost::program_options;


struct Program_Options
{
    size_t n_ops;
    size_t seed;
    bool check_alt_custom;
    bool use_random_scores;
};

void real_main(const Program_Options& po)
{
    srand48(po.seed);
    std::cout << "using seed=" << po.seed << "\n";
    for (size_t i = 0; i < po.n_ops; ++i)
    {
        // construct random strings
        int s_len[2];
        std::string s[2];
        for (int k = 0; k < 2; ++k)
        {
            s_len[k] = 1 + static_cast< size_t >(drand48() * 19);
            for (int j = 0; j < s_len[k]; ++j)
            {
                static const std::string bases = "ACGT";
                s[k].push_back(bases[static_cast< size_t >(drand48() * 4)]);
            }
        }
        // pick random alignment type
        static const std::vector< AffineAlignmentType > alt_types = { ALT_OVERLAP, ALT_GLOBAL, ALT_CONTAINMENT, ALT_CUSTOM };
        static const std::vector< std::string > alt_names = { "ALT_OVERLAP", "ALT_GLOBAL", "ALT_CONTAINMENT", "ALT_CUSTOM" };
        OverlapperParams params = affine_default_params;
        size_t alt_type_idx = static_cast< size_t >(drand48() * (3 + static_cast< int >(po.check_alt_custom)));
        params.type = static_cast< AffineAlignmentType >(alt_type_idx);
        if (params.type == ALT_CUSTOM)
        {
            params.gap_s1_start = static_cast< int >(drand48() * 2);
            params.gap_s1_end = static_cast< int >(drand48() * 2);
            params.gap_s2_start = static_cast< int >(drand48() * 2);
            params.gap_s2_end = static_cast< int >(drand48() * 2);
        }
        if (po.use_random_scores)
        {
            params.match_score = static_cast< int >(drand48() * 20);
            params.mismatch_penalty = - static_cast< int >(drand48() * 20);
            params.gap_penalty = - static_cast< int >(drand48() * 20);
            params.gap_ext_penalty = params.gap_penalty + static_cast< int >(drand48() * (- params.gap_penalty));
        }
        // print settings
        std::cout << "checking alignment: s[0]='" << s[0] << "' s[1]='" << s[1] << "' type="
                  << alt_names[alt_type_idx];
        if (params.type == ALT_CUSTOM)
        {
            std::cout << "("
                      << params.gap_s1_start << ","
                      << params.gap_s1_end << ","
                      << params.gap_s2_start << ","
                      << params.gap_s2_end << ")";
        }
        if (po.use_random_scores)
        {
            std::cout << " scores=("
                      << params.match_score << ","
                      << params.mismatch_penalty << ","
                      << params.gap_penalty << ","
                      << params.gap_ext_penalty << ")";
        }
        std::cout << "\n";

        // construct alignment
        auto output = Overlapper::computeAlignmentAffine(s[0], s[1], params);
        std::cout << "cigar='" << output.cigar << "' score=" << output.score << "\n";

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
        MAC::Cigar cigar(output.cigar);
        assert(static_cast< long >(cigar.get_rf_len()) == output.match[0].end - output.match[0].start + 1);
        assert(static_cast< long >(cigar.get_qr_len()) == output.match[1].end - output.match[1].start + 1);

        // check scores
        cigar.disambiguate(s[0].substr(output.match[0].start, output.match[0].end - output.match[0].start + 1),
                           s[1].substr(output.match[1].start, output.match[1].end - output.match[1].start + 1));
        int score = 0;
        for (size_t j = 0; j < cigar.get_n_ops(); ++j)
        {
            if (cigar.get_op(j) == '=')
            {
                score += cigar.get_op_len(j) * params.match_score;
            }
            else if (cigar.get_op(j) == 'X')
            {
                score += cigar.get_op_len(j) * params.mismatch_penalty;
            }
            else if (cigar.get_op(j) == 'I' || cigar.get_op(j) == 'D')
            {
                score += params.gap_penalty + cigar.get_op_len(j) * params.gap_ext_penalty;
            }
        }
        assert(score == output.score);
    }
    std::cout << "all good\n";
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
        ("n-ops,n", bo::value<size_t>(&po.n_ops)->default_value(1000), "number of operations")
        ("seed,s", bo::value<size_t>(&po.seed)->default_value(0), "random number generator seed")
        ("check-alt-custom,c", bo::bool_switch(&po.check_alt_custom), "check ALT_CUSTOM alignment")
        ("use-random-scores,r", bo::bool_switch(&po.use_random_scores), "use random scores instead of defaults")
        ;
    cmdline_opts_desc.add(generic_opts_desc).add(config_opts_desc).add(hidden_opts_desc);
    visible_opts_desc.add(generic_opts_desc).add(config_opts_desc);
    bo::variables_map vm;
    store(bo::command_line_parser(argc, argv).options(cmdline_opts_desc).run(), vm);
    notify(vm);
    if (vm.count("help"))
    {
        std::cout <<
            "Test affine overlapper with random queries.\n"
            "\n";
        std::cout << visible_opts_desc;
        exit(EXIT_SUCCESS);
    }
    if (po.seed == 0)
    {
        po.seed = time(NULL);
    }
    real_main(po);
}
