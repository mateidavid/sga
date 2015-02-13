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
#include "Validate_Mutations.hpp"

#include "BWT.h"
#include "BWTAlgorithms.h"

using namespace std;
using namespace MAC;
namespace bo = boost::program_options;


void load_asqg(std::istream& is, const bo::variables_map& vm, Graph& g)
{
    LOG("io", info) << ptree("load_asqg_start");
    string line;
    size_t line_count = 0;
    size_t progress_count = vm.at("progress-count").as< unsigned >();
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
        if (vm.at("check-at-step").as< bool >())
        {
            g.check_all();
            g.check_leaks();
        }
        LOG("mac", debug2) << g.to_ptree();
        if (progress_count > 0)
        {
            if ((++line_count % progress_count) == 0)
            {
                LOG("mac", info) << ptree("progress").put("count", line_count);
            }
        }
    }
    LOG("io", info) << ptree("load_asqg_end");
    g.check_all();
    g.check_leaks();
}


int real_main(const bo::variables_map& vm)
{
    Graph g;
    Hap_Map hm;

    // option validation
    if (vm.count("input-file") == vm.count("load-file"))
    {
        LOG("mac", error) << "exactly 1 of input-file or load-file options must be specified" << endl;
        abort();
    }
    if (vm.count("load-file"))
    {
        vector< string > forbidden_options = { "unmap-trigger-len" };
        vector< string > forbidden_bool_switches = { "cat-at-step", "print-at-step", "check-at-step" };
        for (const auto& s : forbidden_options)
        {
            if (vm.count(s))
            {
                LOG("mac", warning) << s << ": ignored when loading a mac graph";
            }
        }
        for (const auto& s : forbidden_bool_switches)
        {
            if (vm.at(s).as< bool >())
            {
                LOG("mac", warning) << s << ": ignored when loading a mac graph";
            }
        }
    }
    if (vm.count("2gs-validate-mutations") and not vm.count("2gs-bwt-file"))
    {
        LOG("mac", error) << "--2gs-bwt-file must be specifed with --2gs-validate-mutations" << endl;
        abort();
    }
    // main work
    if (vm.count("input-file"))
    {
        // load asqg graph
        const string& fn = vm.at("input-file").as< string >();
        LOG("mac", info) << ptree("loading").put("input_file", fn);
        ixstream tmp_fs(fn);
        if (vm.count("unmap-trigger-len"))
        {
            g.unmap_trigger_len() = vm.at("unmap-trigger-len").as< unsigned >();
        }
        g.cat_at_step() = vm.at("cat-at-step").as< bool >();
        load_asqg(tmp_fs, vm, g);
    }
    else
    {
        // load mac graph
        const string& fn = vm.at("load-file").as< string >();
        LOG("mac", info) << ptree("loading").put("load-file", fn);
        fstr ifs(fn, ios_base::in | ios_base::binary);
        g.load(ifs);
    }
    LOG("mac", info) << ptree("factory_stats", g.factory_stats());
    if (vm.at("cat-at-end").as< bool >())
    {
        LOG("mac", info) << ptree("cat_at_end_start");
        g.cat_all_read_contigs();
        g.check_all();
        LOG("mac", info) << ptree("cat_at_end_end");
    }
    if (vm.count("2gs-validate-mutations"))
    {
        const string& fn = vm.at("2gs-bwt-file").as< string >();
        LOG("mac", info) << ptree("validate_mutations_start").put("bwt_file", fn);
        BWTIndexSet index_set;
        index_set.pBWT = new BWT(fn);
        Validate_Mutations()(g, index_set);
        delete index_set.pBWT;
        LOG("mac", info) << ptree("validate_mutations_end");
    }
    if (vm.at("unmap-read-ends").as< bool >())
    {
        LOG("mac", info) << ptree("unmap_read_ends_start");
        g.unmap_read_ends();
        g.check_all();
        LOG("mac", info) << ptree("unmap_read_ends_end");
    }
    if (vm.at("unmap-mut-clusters").as< bool >())
    {
        LOG("mac", info) << ptree("unmap_mut_clusters_start");
        Unmap_Mut_Clusters()(g);
        g.check_all();
        LOG("mac", info) << ptree("unmap_mut_clusters_end");
    }
    if (vm.at("resolve-unmappable-regions").as< bool >())
    {
        LOG("mac", info) << ptree("resolve_unmappable_regions_start");
        g.resolve_unmappable_regions();
        g.check_all();
        LOG("mac", info) << ptree("resolve_unmappable_regions_end");
    }
    if (vm.at("unmap-single-chunks").as< bool >())
    {
        LOG("mac", info) << ptree("unmap_single_chunks_start");
        g.unmap_single_chunks();
        g.check_all();
        LOG("mac", info) << ptree("unmap_single_chunks_end");
    }
    if (vm.at("hap-map").as< bool >())
    {
        LOG("mac", info) << ptree("hap_map_start");
        hm.build(g);
        LOG("mac", info) << ptree("hap_map_end");
    }
    if (vm.at("interactive").as< bool >())
    {
        cli(std::cin, std::cout, g, hm);
    }
    if (vm.count("stats-file"))
    {
        const string& fn = vm.at("stats-file").as< string >();
        LOG("mac", info) << ptree("print_stats").put("stats_file", fn);
        fstr tmp_fs(fn, ios_base::out);
        g.dump_detailed_counts(tmp_fs);
    }
    if (vm.count("supercontig-stats-file"))
    {
        const string& fn = vm.at("supercontig-stats-file").as< string >();
        LOG("mac", info) << ptree("supercontig_stats").put("file", fn);
        fstr tmp_fs(fn, ios_base::out);
        g.print_supercontig_stats(tmp_fs);
    }
    if (vm.count("mutations-file"))
    {
        const string& fn = vm.at("mutations-file").as< string >();
        LOG("mac", info) << ptree("print_mutations").put("file", fn);
        fstr tmp_fs(fn, ios_base::out);
        g.print_mutations(tmp_fs);
    }
    if (vm.count("terminal-reads-file"))
    {
        const string& fn = vm.at("terminal-reads-file").as< string >();
        LOG("mac", info) << ptree("print_terminal_reads").put("file", fn);
        fstr tmp_fs(fn, ios_base::out);
        g.get_terminal_reads(tmp_fs);
    }
    if (vm.at("print-at-end").as< bool >())
    {
        cout << g.to_ptree();
    }
    if (vm.count("unmappable-contigs-file"))
    {
        const string& fn = vm.at("unmappable-contigs-file").as< string >();
        LOG("mac", info) << ptree("print_unmappable_contigs").put("file", fn);
        fstr tmp_fs(fn, ios_base::out);
        g.print_unmappable_contigs(tmp_fs);
    }
    if (vm.count("hapmap-stats-file"))
    {
        const string& fn = vm.at("hapmap-stats-file").as< string >();
        LOG("mac", info) << ptree("hapmap-stats").put("file", fn);
        fstr tmp_fs(fn, ios_base::out);
        hm.dump_stats(tmp_fs, g);
    }
    LOG("mac", info) << ptree("done_output");

    if (vm.count("save-file"))
    {
        const string& fn = vm.at("save-file").as< string >();
        LOG("mac", info) << ptree("saving").put("save_file", fn);
        fstr tmp_fs(fn, ios_base::out | ios_base::binary);
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
        ("log-level,d", bo::value< vector< string > >()->default_value(vector< string >(), ""), "log level")
        ("seed", bo::value< long >()->default_value(time(nullptr), "use time"), "random seed")
        ("progress-count,c", bo::value< unsigned >()->default_value(0), "progress count")
        ;
    config_opts_desc.add_options()
        //
        // file-related options
        //
        ("input-file,i", bo::value< string >(), "input file")
        ("load-file,L", bo::value< string >(), "load file")
        ("save-file,S", bo::value< string >(), "save file")
        ("stats-file", bo::value< string >(), "stats file")
        ("supercontig-stats-file", bo::value< string >(), "supercontig stats file")
        ("mutations-file", bo::value< string >(), "mutations file")
        ("unmappable-contigs-file", bo::value< string >(), "unmappable contigs file")
        ("terminal-reads-file", bo::value< string >(), "terminal reads file")
        ("hapmap-stats-file", bo::value< string >(), "haplotype map stats file")
        ("2gs-bwt-file", bo::value< string >(), "BWT index of 2GS reads that can be used to validate mutations")
        //
        // asqg loading options
        //
        ("unmap-trigger-len,u", bo::value< unsigned >(), "unmap trigger len")
        ("cat-at-step,s", bo::bool_switch(), "cat contigs at each step")
        ("print-at-step", bo::bool_switch(), "print graph at each step")
        ("check-at-step", bo::bool_switch(), "check graph at each step")
        //
        // post-loading options
        //
        ("cat-at-end", bo::bool_switch(), "cat contigs at end")
        ("print-at-end", bo::bool_switch(), "print graph at end")
        ("unmap-read-ends", bo::bool_switch(), "unmap read ends")
        ("resolve-unmappable-regions", bo::bool_switch(), "resolve unmappable regions")
        ("2gs-validate-mutations", bo::bool_switch(), "validate mutations using 2GS data")
        ("unmap-single-chunks", bo::bool_switch(), "unmap single chunks")
        ("unmap-mut-clusters", bo::bool_switch(), "unmap mutation clusters")
        ("hap-map", bo::bool_switch(), "build haplotype map")
        ("interactive", bo::bool_switch(), "run interactive commands")
        ;
    any_converter ac;
    ac.add_string_converter< string >();
    ac.add_string_converter< bool >();
    ac.add_string_converter< unsigned >();
    ac.add_string_converter< long >();
    ac.add_converter(&cont_to_ptree< vector< string > >);
    cmdline_opts_desc.add(generic_opts_desc).add(config_opts_desc).add(hidden_opts_desc);
    visible_opts_desc.add(generic_opts_desc).add(config_opts_desc);
    bo::variables_map vm;
    store(bo::command_line_parser(argc, argv).options(cmdline_opts_desc).run(), vm);
    notify(vm);

    // if help requested, print it and stop
    if (vm.count("help"))
    {
        cout << visible_opts_desc;
        exit(EXIT_SUCCESS);
    }
    // set log levels
    for (const auto& l : vm.at("log-level").as< vector< string > >())
    {
        size_t i = l.find(':');
        if (i == string::npos)
        {
            Logger::set_default_level(l);
            clog << "set default log level to: "
                 << static_cast< int >(Logger::get_default_level()) << "\n";
        }
        else
        {
            Logger::set_facility_level(l.substr(0, i), l.substr(i + 1));
            clog << "set log level of '" << l.substr(0, i) << "' to: "
                 << static_cast< int >(Logger::get_facility_level(l.substr(0, i))) << "\n";
        }
    }
    // set random seed
    srand48(vm.at("seed").as< long >());

    // print options
    LOG("mac", info) << variables_map_converter::to_ptree(vm, ac);

    return real_main(vm);
}
