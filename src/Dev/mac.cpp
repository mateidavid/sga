#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <list>
#include <boost/program_options.hpp>

#include "Graph.hpp"
#include "Hap_Map.hpp"
#include "ixstream.hpp"
#include "fstr.hpp"
#include "logger.hpp"
#include "variables_map_converter.hpp"
#include "CLI.hpp"

#include "Unmap_Mut_Clusters.hpp"
#include "Validate_Variations.hpp"

#include "BWT.h"
#include "BWTAlgorithms.h"

using namespace std;
using namespace MAC;
namespace bo = boost::program_options;

/// Program options.
/// Global variables storing program options. Their default values are specified
/// in the option descriptions.
namespace opts
{
    // main variable map
    bo::variables_map vm;
    //
    // generic
    //
    vector< string > log_level;
    long seed;
    unsigned progress_count;
    //
    // file-related
    //
    string input_file;
    string load_file;
    string save_file;
    string stats_file;
    string supercontigs_stats_file;
    string mutations_file;
    string unmappable_contigs_file;
    string terminal_reads_file;
    string hapmap_stats_file;
    string bwt_file;
    string aux_bwt_file;
    //
    // asqg loading
    //
    unsigned unmap_trigger_len;
    bool cat_at_step;
    bool print_at_step;
    bool check_at_step;
    //
    // post-loading
    //
    bool cat_at_end;
    bool print_at_end;
    bool unmap_read_ends;
    bool resolve_unmappable_regions;
    bool validate_variations;
    bool unmap_single_chunks;
    bool unmap_mut_clusters;
    bool build_hap_map;
    //
    // other
    //
    bool interactive;
} // namespace opts


void load_asqg(std::istream& is, Graph& g)
{
    LOG("io", info) << ptree("load_asqg_start");
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
            global::assert_message() = string("VT ") + r_id;
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
            global::assert_message() = string("ED ") + r1_id + " " + r2_id;
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
                LOG("mac", info) << ptree("progress").put("count", line_count);
            }
        }
    }
    LOG("io", info) << ptree("load_asqg_end");
    g.check_all();
    g.check_leaks();
}

int real_main()
{
    Graph g;
    Hap_Map hm;
    // option validation
    if (opts::input_file.empty() == opts::load_file.empty())
    {
        LOG("mac", error) << "exactly 1 of input-file or load-file options must be specified" << endl;
        abort();
    }
    if (not opts::load_file.empty())
    {
        vector< string > forbidden_options = { "unmap-trigger-len" };
        vector< string > forbidden_bool_switches = { "cat-at-step", "print-at-step", "check-at-step" };
        for (const auto& s : forbidden_options)
        {
            if (opts::vm.count(s))
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
    }
    if (opts::validate_variations and opts::aux_bwt_file.empty())
    {
        LOG("mac", error) << "--aux-bwt-file must be specifed with --validate-mutations" << endl;
        abort();
    }
    // main work
    if (not opts::input_file.empty())
    {
        // load asqg graph
        LOG("mac", info) << ptree("loading").put("input_file", opts::input_file);
        ixstream tmp_fs(opts::input_file);
        g.unmap_trigger_len() = opts::unmap_trigger_len;
        g.cat_at_step() = opts::cat_at_step;
        load_asqg(tmp_fs, g);
    }
    else
    {
        // load mac graph
        LOG("mac", info) << ptree("loading").put("load-file", opts::load_file);
        fstr ifs(opts::load_file, ios_base::in | ios_base::binary);
        g.load(ifs);
    }
    LOG("mac", info) << ptree("factory_stats", g.factory_stats());
    if (opts::cat_at_end)
    {
        LOG("mac", info) << ptree("cat_at_end_start");
        g.cat_all_read_contigs();
        g.check_all();
        LOG("mac", info) << ptree("cat_at_end_end");
    }
    if (opts::validate_variations)
    {
        BWTIndexSet index_set;
        LOG("mac", info) << ptree("validate_variations__load_index__start").put("aux_bwt_file", opts::aux_bwt_file);
        index_set.pBWT = new BWT(opts::aux_bwt_file);
        LOG("mac", info) << ptree("validate_variations__load_index__end");
        Validate_Variations()(g, index_set);
        LOG("mac", info) << ptree("validate_variations__delete_index__start");
        delete index_set.pBWT;
        LOG("mac", info) << ptree("validate_variations__delete_index__end");
    }
    if (opts::unmap_read_ends)
    {
        LOG("mac", info) << ptree("unmap_read_ends_start");
        g.unmap_read_ends();
        g.check_all();
        LOG("mac", info) << ptree("unmap_read_ends_end");
    }
    if (opts::unmap_mut_clusters)
    {
        LOG("mac", info) << ptree("unmap_mut_clusters_start");
        Unmap_Mut_Clusters()(g);
        g.check_all();
        LOG("mac", info) << ptree("unmap_mut_clusters_end");
    }
    if (opts::resolve_unmappable_regions)
    {
        LOG("mac", info) << ptree("resolve_unmappable_regions_start");
        g.resolve_unmappable_regions();
        g.check_all();
        LOG("mac", info) << ptree("resolve_unmappable_regions_end");
    }
    if (opts::unmap_single_chunks)
    {
        LOG("mac", info) << ptree("unmap_single_chunks_start");
        g.unmap_single_chunks();
        g.check_all();
        LOG("mac", info) << ptree("unmap_single_chunks_end");
    }
    if (opts::build_hap_map)
    {
        LOG("mac", info) << ptree("hap_map_start");
        hm.build(g);
        LOG("mac", info) << ptree("hap_map_end");
    }
    if (opts::interactive)
    {
        cli(std::cin, std::cout, g, hm);
    }
    if (not opts::stats_file.empty())
    {
        LOG("mac", info) << ptree("print_stats").put("stats_file", opts::stats_file);
        fstr tmp_fs(opts::stats_file, ios_base::out);
        g.dump_detailed_counts(tmp_fs);
    }
    if (not opts::supercontigs_stats_file.empty())
    {
        LOG("mac", info) << ptree("supercontigs_stats").put("file", opts::supercontigs_stats_file);
        fstr tmp_fs(opts::supercontigs_stats_file, ios_base::out);
        g.print_supercontig_stats(tmp_fs);
    }
    if (not opts::mutations_file.empty())
    {
        LOG("mac", info) << ptree("print_mutations").put("file", opts::mutations_file);
        fstr tmp_fs(opts::mutations_file, ios_base::out);
        g.print_mutations(tmp_fs);
    }
    if (not opts::terminal_reads_file.empty())
    {
        LOG("mac", info) << ptree("print_terminal_reads").put("file", opts::terminal_reads_file);
        fstr tmp_fs(opts::terminal_reads_file, ios_base::out);
        g.get_terminal_reads(tmp_fs);
    }
    if (opts::print_at_end)
    {
        cout << g.to_ptree();
    }
    if (not opts::unmappable_contigs_file.empty())
    {
        LOG("mac", info) << ptree("print_unmappable_contigs").put("file", opts::unmappable_contigs_file);
        fstr tmp_fs(opts::unmappable_contigs_file, ios_base::out);
        g.print_unmappable_contigs(tmp_fs);
    }
    if (not opts::hapmap_stats_file.empty())
    {
        LOG("mac", info) << ptree("hapmap-stats").put("file", opts::hapmap_stats_file);
        fstr tmp_fs(opts::hapmap_stats_file, ios_base::out);
        hm.dump_stats(tmp_fs, g);
    }
    LOG("mac", info) << ptree("done_output");

    if (not opts::save_file.empty())
    {
        LOG("mac", info) << ptree("saving").put("save_file", opts::save_file);
        fstr tmp_fs(opts::save_file, ios_base::out | ios_base::binary);
        g.save(tmp_fs);
    }
    LOG("mac", info) << ptree("factory_stats", g.factory_stats());
    g.clear_and_dispose();
    hm.clear_and_dispose();
    LOG("mac", info) << ptree("graph_cleared");
    return EXIT_SUCCESS;
}

int main(int argc, char* argv[])
{
    global::program_name() = argv[0];

    bo::options_description generic_opts_desc("Generic options");
    bo::options_description config_opts_desc("Configuration options");
    bo::options_description hidden_opts_desc("Hidden options");
    bo::options_description cmdline_opts_desc;
    bo::options_description visible_opts_desc("Allowed options");
    generic_opts_desc.add_options()
        ("help,?", "produce help message")
        // hack, see: http://lists.boost.org/boost-users/2010/01/55054.php
        ("log-level,d",
             bo::value(&opts::log_level)->default_value(vector< string >(), ""),
                 "log level")
        ("seed",
             bo::value(&opts::seed)->default_value(time(nullptr), "use time"),
                 "random seed")
        ("progress-count,c",
             bo::value(&opts::progress_count)->default_value(0),
                 "progress count")
        ;
    config_opts_desc.add_options()
        //
        // file-related options
        //
        ("input-file,i",
             bo::value(&opts::input_file),
                 "input file")
        ("load-file,L",
             bo::value(&opts::load_file),
                 "load file")
        ("save-file,S",
             bo::value(&opts::save_file),
                 "save file")
        ("stats-file",
             bo::value(&opts::stats_file),
                 "stats file")
        ("supercontigs-stats-file",
             bo::value(&opts::supercontigs_stats_file),
                 "supercontig stats file")
        ("mutations-file",
             bo::value(&opts::mutations_file),
                 "mutations file")
        ("unmappable-contigs-file",
             bo::value(&opts::unmappable_contigs_file),
                 "unmappable contigs file")
        ("terminal-reads-file",
             bo::value(&opts::terminal_reads_file),
                 "terminal reads file")
        ("hapmap-stats-file",
             bo::value(&opts::hapmap_stats_file),
                 "haplotype map stats file")
        ("aux-bwt-file",
             bo::value(&opts::aux_bwt_file),
                 "BWT index of 2GS reads used to validate variations")
        //
        // asqg loading options
        //
        ("unmap-trigger-len,u",
             bo::value(&opts::unmap_trigger_len)->default_value(9),
                 "unmap trigger len")
        ("cat-at-step,s",
             bo::bool_switch(&opts::cat_at_step),
                 "cat contigs at each step")
        ("print-at-step",
             bo::bool_switch(&opts::print_at_step),
                 "print graph at each step")
        ("check-at-step",
             bo::bool_switch(&opts::check_at_step),
                 "check graph at each step")
        //
        // post-loading options
        //
        ("cat-at-end",
             bo::bool_switch(&opts::cat_at_end),
                 "cat contigs at end")
        ("print-at-end",
             bo::bool_switch(&opts::print_at_end),
                 "print graph at end")
        ("unmap-read-ends",
             bo::bool_switch(&opts::unmap_read_ends),
                 "unmap read ends")
        ("resolve-unmappable-regions",
             bo::bool_switch(&opts::resolve_unmappable_regions),
                 "resolve unmappable regions")
        ("validate-variations",
             bo::bool_switch(&opts::validate_variations),
                 "validate variations with auxiliary 2GS data")
        ("unmap-single-chunks",
             bo::bool_switch(&opts::unmap_single_chunks),
                 "unmap single chunks")
        ("unmap-mut-clusters",
             bo::bool_switch(&opts::unmap_mut_clusters),
                 "unmap mutation clusters")
        ("hap-map",
             bo::bool_switch(&opts::build_hap_map),
                 "build haplotype map")
        //
        // other
        //
        ("interactive",
             bo::bool_switch(&opts::interactive),
                 "run interactive commands")
        ;
    any_converter ac;
    ac.add_string_converter< string >();
    ac.add_string_converter< bool >();
    ac.add_string_converter< unsigned >();
    ac.add_string_converter< long >();
    ac.add_converter(&cont_to_ptree< vector< string > >);
    cmdline_opts_desc.add(generic_opts_desc).add(config_opts_desc).add(hidden_opts_desc);
    visible_opts_desc.add(generic_opts_desc).add(config_opts_desc);
    try
    {
        store(bo::command_line_parser(argc, argv).options(cmdline_opts_desc).run(), opts::vm);
        // if help requested, print it and stop
        if (opts::vm.count("help"))
        {
            //usage(cout);
            cout << visible_opts_desc;
            exit(EXIT_SUCCESS);
        }
        notify(opts::vm);
    }
    catch(bo::error& e) 
    { 
        cerr << "ERROR: " << e.what() << endl << endl;
        //usage(cerr);
        cerr << visible_opts_desc;
        exit(EXIT_FAILURE);
    }
    // set log levels
    Logger::set_levels_from_options(opts::log_level, &clog);
    // set random seed
    srand48(opts::seed);
    // print options
    LOG("mac", info) << variables_map_converter::to_ptree(opts::vm, ac);
    return real_main();
}
