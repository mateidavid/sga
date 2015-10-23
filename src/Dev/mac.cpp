#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <list>
#include <tclap/CmdLine.h>

#include "version.hpp"
#include "zstr.hpp"
#include "Graph.hpp"
//#include "Hap_Map.hpp"
#include "logger.hpp"
#include "CLI.hpp"

#include "Variation_Cleaner.hpp"
#include "Validate_Variations.hpp"
#include "Read_Merger.hpp"

#include "BWT.h"
#include "BWTAlgorithms.h"

using namespace std;
using namespace MAC;

/// Program options.
/// Global variables storing program options. Their default values are specified
/// in the option descriptions.
namespace opts
{
    using namespace TCLAP;
    string prog_desc =
        "Multi-Allelic Contig assembler";
    CmdLine cmd_parser(prog_desc, ' ', package_version);
    //
    // generic
    //
    MultiArg< string > log_level("d", "log-level", "Log level.", false, "string", cmd_parser);
    ValueArg< long > seed("", "seed", "Random seed (-1: use time).", false, -1, "int", cmd_parser);
    ValueArg< unsigned > progress_count("c", "progress-count", "Progress count.", false, 0, "int", cmd_parser);
    //
    // file-related
    //
    ValueArg< string > input_file("i", "input-file", "Load ASQG input file.", false, "", "file", cmd_parser);
    ValueArg< string > load_file("L", "load-file", "Load MAC graph from file.", false, "", "file", cmd_parser);
    ValueArg< string > save_file("S", "save-file", "Save MAC graph to file.", false, "", "file", cmd_parser);
    ValueArg< string > gfa_file("", "export-gfa", "Export graph to GFA file.", false, "", "file", cmd_parser);
    ValueArg< string > bwt_prefix("", "bwt-prefix", "Prefix of {.bwt,.ssa,.id.gz} files representing BWT index, sampled suffix array, and id list of the reads.", false, "", "file", cmd_parser);
    ValueArg< string > aux_bwt_file("", "aux-bwt-file", "BWT index of 2GS reads used to validate variations.", false, "", "file", cmd_parser);
    ValueArg< double > aux_coverage("", "aux-coverage", "Base coverage by 2GS reads.", false, -1, "float", cmd_parser);
    ValueArg< unsigned > aux_read_len("", "aux-read-len", "Read length of 2GS reads.", false, 1000, "int", cmd_parser);
    ValueArg< string > stats_file("", "stats-file", "Stats file.", false, "", "file", cmd_parser);
    ValueArg< string > supercontigs_stats_file("", "supercontigs-stats-file", "Supercontigs stats file.", false, "", "file", cmd_parser);
    ValueArg< string > supercontigs_file("", "supercontigs-file", "Supercontigs file.", false, "", "file", cmd_parser);
    ValueArg< string > unmappable_contigs_file("", "unmappable-contigs-file", "Unmappable contigs file.", false, "", "file", cmd_parser);
    ValueArg< string > terminal_reads_file("", "terminal-reads-file", "Terminal reads file.", false, "", "file", cmd_parser);
    //ValueArg< string > hapmap_stats_file("", "hapmap-stats-file", "Haplotype map stats file.", false, "", "file", cmd_parser);
    //
    // asqg loading
    //
    ValueArg< unsigned > unmap_trigger_len("u", "unmap-trigger-len", "Trigger length for unmapping.", false, 9, "int", cmd_parser);
    SwitchArg cat_at_step("s", "cat-at-step", "Catenate contigs after each operation.", cmd_parser, false);
    SwitchArg print_at_step("", "print-at-step", "Print graph after each operation.", cmd_parser, false);
    SwitchArg check_at_step("", "check-at-step", "Check graph after each operation.", cmd_parser, false);
    SwitchArg trim_tuc_step("", "trim-tuc-step", "Trim terminal unmappable chunks during unmapping.", cmd_parser, false);
    //
    // post-loading
    //
    SwitchArg cat_at_end("", "cat-at-end", "Catenate contigs after loading graph.", cmd_parser, false);
    SwitchArg trim_tuc_end("", "trim-tuc-end", "Trim terminal unmappable chunks at the end.", cmd_parser, false);
    SwitchArg print_at_end("", "print-at-end", "Print graph after loading graph.", cmd_parser, false);
    SwitchArg validate_variations("", "validate-variations", "Validate variations with auxiliary 2GS data.", cmd_parser, false);
    SwitchArg unmap_single_chunks("", "unmap-single-chunks", "Unmap read chunks not supported by other reads.", cmd_parser, false);
    SwitchArg unmap_single_terminal_regions("", "unmap-single-terminal-regions", "Unmap terminal read regions not supported by other reads.", cmd_parser, false);
    ValueArg< int > unmap_short_contigs("", "unmap-short-contigs", "Unmap contigs smaller than a given size.", false, 0, "int", cmd_parser);
    SwitchArg clean_variations("", "clean-variations", "Cleanup variations.", cmd_parser, false);
    SwitchArg unmap_homopolymer_indels("", "unmap-homopolymer-indels", "Unmap indels around homopolymers.", cmd_parser, false);
    SwitchArg copy_num("", "copy-num", "Compute copy numbers.", cmd_parser, false);
    SwitchArg merge_reads_1("", "merge-reads-1", "Merge reads.", cmd_parser, false);
    SwitchArg merge_reads_2("", "merge-reads-2", "Merge reads.", cmd_parser, false);
    ValueArg< unsigned > merged_weight("", "merged-weight", "Weight of merged reads.", false, 5, "int", cmd_parser);
    SwitchArg gfa_show_mutations("", "gfa-show-mutations", "Show mutations in GFA output.", cmd_parser, false);
    //
    // other
    //
    SwitchArg interactive("", "interactive", "Enter interactive mode.", cmd_parser, false);
    SwitchArg test_1("", "test-1", "Internal switch.", cmd_parser, false);
    SwitchArg test_2("", "test-2", "Internal switch.", cmd_parser, false);
    SwitchArg test_3("", "test-3", "Internal switch.", cmd_parser, false);
} // namespace opts


void load_asqg(std::istream& is, Graph& g)
{
    LOG("io", info) << ptree("begin");
    string line;
    size_t line_count = 0;
    while (getline(is, line))
    {
        LOG("mac", debug) << ptree("op").put("line", line);
        istringstream iss(line + "\n");
        string rec_type;
        iss >> rec_type;
        if (iss.eof() or rec_type[0] == '#')
        {
            // ignore blank lines & lines starting with '#'
        }
        if (rec_type == "HT")
        {
            // ignore header line
        }
        else if (rec_type == "VT")
        {
            string r_id;
            string r_seq;
            iss >> r_id >> r_seq;
            ASSERT(not iss.eof());
            global_assert::global_msg() = string("VT ") + r_id;
            g.add_read(std::move(r_id), std::move(r_seq));
        }
        else if (rec_type == "ED")
        {
            string r1_id;
            string r2_id;
            size_t r1_start, r1_end, r1_len;
            size_t r2_start, r2_end, r2_len;
            int rc, tmp;
            string sam_cigar;
            string sam_pi;
            iss >> r1_id >> r2_id
            >> r1_start >> r1_end >> r1_len
            >> r2_start >> r2_end >> r2_len
            >> rc >> tmp >> sam_cigar >> sam_pi;
            ASSERT(not iss.eof());
            // switch to open interval ends
            ++r1_end;
            ++r2_end;
            global_assert::global_msg() = string("ED ") + r1_id + " " + r2_id;
            g.add_overlap(r1_id, r2_id, r1_start, r1_end - r1_start, r2_start, r2_end - r2_start, rc, sam_cigar.substr(5));
        }
        if (opts::check_at_step)
        {
            g.check_all();
            g.check_leaks();
        }
        LOG("mac", debug2) << g.to_ptree();
        if (opts::progress_count > 0)
        {
            if ((++line_count % opts::progress_count) == 0)
            {
                LOG("mac", info) << ptree("progress_count").put("count", line_count);
            }
        }
    }
    LOG("io", info) << ptree("end");
    g.check_all();
    g.check_leaks();
}

int real_main()
{
    Graph g;
    //Hap_Map hm;
    // option validation
    if (opts::input_file.get().empty() == opts::load_file.get().empty())
    {
        LOG("mac", error) << "exactly 1 of input-file or load-file options must be specified" << endl;
        exit(EXIT_FAILURE);
    }
    if (not opts::load_file.get().empty())
    {
        /* TODO:
        vector< string > forbidden_options = { "unmap-trigger-len" };
        vector< string > forbidden_bool_switches = { "cat-at-step", "print-at-step", "check-at-step" };
        for (const auto& s : forbidden_options)
        {
            if (opts::cat_at_step)
            {
                LOG("mac", warning) << s << ": ignored when loading a mac graph";
            }
        }
        for (const auto& s : forbidden_bool_switches)
        {
            if (opts::vm.at(s).as< bool >())
            {
                LOG("mac", warning) << s << ": ignored when loading a mac graph";
            }
        }
        */
    }
    if (opts::validate_variations and opts::aux_bwt_file.get().empty())
    {
        LOG("mac", error) << "--aux-bwt-file must be specifed with --validate-variations" << endl;
        exit(EXIT_FAILURE);
    }
    // main work
    if (not opts::input_file.get().empty())
    {
        // load asqg graph
        LOG("mac", info) << ptree("load_asqg").put("file", opts::input_file.get());
        zstr::ifstream tmp_fs(opts::input_file);
        load_asqg(tmp_fs, g);
    }
    else
    {
        // load mac graph
        LOG("mac", info) << ptree("load_mac").put("file", opts::load_file.get());
        strict_fstream::fstream ifs(opts::load_file, ios_base::in | ios_base::binary);
        g.load(ifs);
    }
    g.unmap_trigger_len() = opts::unmap_trigger_len;
    g.cat_at_step() = opts::cat_at_step;
    g.trim_tuc_step() = opts::trim_tuc_step;
    LOG("mac", info) << ptree().put("factory_stats", g.factory_stats());
    if (opts::cat_at_end)
    {
        g.cat_all_read_contigs();
        g.check_all();
    }
    if (opts::trim_tuc_end)
    {
        g.trim_tucs();
        g.check_all();
    }
    if (not opts::bwt_prefix.get().empty())
    {
        g.load_bwt(opts::bwt_prefix);
    }
    if (not opts::aux_bwt_file.get().empty())
    {
        g.load_aux_bwt(opts::aux_bwt_file, opts::aux_read_len, opts::aux_coverage);
    }

    if (opts::validate_variations)
    {
        Validate_Variations{g}();
    }
    if (opts::unmap_single_terminal_regions)
    {
        g.unmapper().unmap_single_terminal_regions();
    }
    if (opts::clean_variations)
    {
        Variation_Cleaner{g}();
    }
    if (opts::unmap_single_chunks)
    {
        g.unmapper().unmap_single_chunks();
    }
    if (opts::unmap_short_contigs > 0)
    {
        g.unmapper().unmap_short_contigs(opts::unmap_short_contigs, 1);
    }
    if (opts::unmap_homopolymer_indels)
    {
        g.unmapper().unmap_homopolymer_indels();
        if (opts::clean_variations)
        {
            Variation_Cleaner{g}();
        }
    }
    if (opts::copy_num)
    {
        g.compute_unique_flank_lens(10, 100);
        g.compute_mutation_copy_num();
    }
    if (opts::merge_reads_1)
    {
        Read_Merger{g, 1, opts::merged_weight}.merge_haploid_alleles();
    }
    if (opts::merge_reads_2)
    {
        Read_Merger{g, 1, opts::merged_weight}.remove_contained();
    }
    if (opts::test_1)
    {
        g.test_mutation_allele_swapping();
    }
    if (opts::interactive)
    {
        //cli(std::cin, std::cout, g, hm);
        cli(std::cin, std::cout, g);
    }
    if (not opts::stats_file.get().empty())
    {
        LOG("mac", info) << ptree("print_stats").put("file", opts::stats_file.get());
        strict_fstream::fstream tmp_fs(opts::stats_file, ios_base::out);
        g.print_detailed_counts(tmp_fs);
    }
    else
    {
        g.print_basic_stats(clog);
    }
    if (not opts::supercontigs_stats_file.get().empty())
    {
        LOG("mac", info) << ptree("supercontigs_stats").put("file", opts::supercontigs_stats_file.get());
        strict_fstream::fstream tmp_fs(opts::supercontigs_stats_file, ios_base::out);
        g.print_supercontig_stats(tmp_fs);
    }
    if (not opts::supercontigs_file.get().empty())
    {
        LOG("mac", info) << ptree("supercontigs").put("file", opts::supercontigs_file.get());
        strict_fstream::fstream tmp_fs(opts::supercontigs_file, ios_base::out);
        g.print_supercontigs(tmp_fs, true, 1);
    }
    if (not opts::terminal_reads_file.get().empty())
    {
        LOG("mac", info) << ptree("print_terminal_reads").put("file", opts::terminal_reads_file.get());
        strict_fstream::fstream tmp_fs(opts::terminal_reads_file, ios_base::out);
        g.get_terminal_reads(tmp_fs);
    }
    if (opts::print_at_end)
    {
        cout << g.to_ptree();
    }
    if (not opts::unmappable_contigs_file.get().empty())
    {
        LOG("mac", info) << ptree("print_unmappable_contigs").put("file", opts::unmappable_contigs_file.get());
        strict_fstream::fstream tmp_fs(opts::unmappable_contigs_file, ios_base::out);
        g.print_unmappable_contigs(tmp_fs);
    }
    /*
    if (not opts::hapmap_stats_file.get().empty())
    {
        LOG("mac", info) << ptree("hapmap-stats").put("file", opts::hapmap_stats_file.get());
        strict_fstream::fstream tmp_fs(opts::hapmap_stats_file, ios_base::out);
        hm.dump_stats(tmp_fs, g);
    }
    */
    LOG("mac", info) << ptree("done_output");

    if (not opts::save_file.get().empty())
    {
        LOG("mac", info) << ptree("save_mac").put("file", opts::save_file.get());
        strict_fstream::fstream tmp_fs(opts::save_file, ios_base::out | ios_base::binary);
        g.save(tmp_fs);
    }
    if (not opts::gfa_file.get().empty())
    {
        LOG("mac", info) << ptree("export_gfa").put("file", opts::gfa_file.get());
        strict_fstream::fstream tmp_fs(opts::gfa_file, ios_base::out);
        g.export_gfa(tmp_fs, opts::gfa_show_mutations);
    }

    LOG("mac", info) << ptree().put("factory_stats", g.factory_stats());
    g.clear_and_dispose();
    //hm.clear_and_dispose();
    LOG("mac", info) << ptree("graph_cleared");
    return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
    opts::cmd_parser.parse(argc, argv);
    global_assert::prog_name() = opts::cmd_parser.getProgramName();
    // set log levels
    Logger::set_default_level(Logger::level::info);
    Logger::set_levels_from_options(opts::log_level, &clog);
    // print options
    LOG("main", info) << "program: " << opts::cmd_parser.getProgramName() << endl
                      << "version: " << opts::cmd_parser.getVersion() << endl
                      << "args: " << opts::cmd_parser.getOrigArgv() << endl;
    // set random seed
    if (opts::seed < 0)
    {
        opts::seed.get() = time(nullptr);
        LOG("main", info) << "seed: " << opts::seed << endl;
    }
    srand48(opts::seed);
    return real_main();
}
