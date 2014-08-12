#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <list>
#include <boost/program_options.hpp>

#include "globals.hpp"
#include "global_assert.hpp"
#include "shortcuts.hpp"
#include "Graph.hpp"
#include "indent.hpp"
#include "ixstream.hpp"
#include "logger.hpp"

#define LOG_FACILITY "main"


using namespace std;
using namespace MAC;
namespace bo = boost::program_options;

/*
auto inc_tab = indent::inc;
auto dec_tab = indent::dec;
using indent::nl;
using indent::tab;
*/

struct Program_Options
{
    string input_file;
    string stats_file_1;
    string stats_file_2;
    string supercontig_lengths_file;
    string mutations_file;
    string unmappable_contigs_file;
    size_t progress_count;
    size_t unmap_trigger_len;
    size_t default_log_level;
    bool cat_at_step;
    bool cat_at_end;
    bool print_at_step;
    bool print_at_end;

    boost::property_tree::ptree to_ptree() const
    {
        return ptree().put("input file", input_file)
                      .put("stats file 1", stats_file_1)
                      .put("stats file 2", stats_file_2)
                      .put("supercontig lengths file", supercontig_lengths_file)
                      .put("mutations file", mutations_file)
                      .put("unmappable contigs file", unmappable_contigs_file)
                      .put("progress count", progress_count)
                      .put("unmap trigger length", unmap_trigger_len)
                      .put("default log level", default_log_level)
                      .put("cat contigs at each step", cat_at_step)
                      .put("cat contigs at end", cat_at_end)
                      .put("print graph at each step", print_at_step)
                      .put("print graph at end", print_at_end)
                      ;
    }
};

int real_main(const Program_Options& po)
{
    log_l(info) << ptree().put("main() settings", po.to_ptree());

    if (po.input_file.empty())
    {
        cerr << "no input file\n";
        exit(EXIT_FAILURE);
    }

    Graph g;
    ixstream ixs(po.input_file);
    if (not ixs)
    {
        cerr << "error opening file: " << po.input_file << '\n';
        return EXIT_FAILURE;
    }

    string line;
    size_t line_count = 0;
    while (getline(ixs, line))
    {
        log_l(debug) << ptree().put("main() op", line);

        istringstream iss(line + "\n");
        string rec_type;
        iss >> rec_type;
        if (iss.eof() or rec_type[0] == '#')
        {
            // ignore blank lines & lines starting with '#'
        }
        if (rec_type == "HT")
        {
            // ignore header line for now
        }
        else if (rec_type == "VT")
        {
            string r_id;
            string r_seq;
            iss >> r_id >> r_seq;
            ASSERT(not iss.eof());
            global::assert_message = string("VT ") + r_id;
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
            global::assert_message = string("ED ") + r1_id + " " + r2_id;
            g.add_overlap(r1_id, r2_id, r1_start, r1_end - r1_start, r2_start, r2_end - r2_start, rc, sam_cigar.substr(5), po.cat_at_step);
        }

        log_l(debug2) << g.to_ptree();

        if (po.progress_count > 0)
        {
            if ((++line_count % po.progress_count) == 0)
            {
                log_l(info) << ptree().put("count", line_count);
            }
        }
        //ASSERT(g.check_all());
    }
    log_l(info) << ptree().put("main() done loop", "");
    g.unmap_single_chunks();
    g.set_contig_ids();
    ASSERT(g.check_all());
    log_l(info) << ptree().put("main() done postprocessing", "");
    if (po.cat_at_end)
    {
        if (not po.stats_file_2.empty())
        {
            ofstream stats_2_ofs(po.stats_file_2);
            g.dump_detailed_counts(stats_2_ofs);
        }
        g.cat_all_read_contigs();
        ASSERT(g.check_all());
    }
    if (not po.stats_file_1.empty())
    {
        ofstream stats_1_ofs(po.stats_file_1);
        g.dump_detailed_counts(stats_1_ofs);
    }
    if (not po.supercontig_lengths_file.empty())
    {
        //TODO
        //ofstream lengths_file(po.supercontig_lengths_file);
        //g.print_supercontig_lengths_2(lengths_file);
        //ASSERT(g.check_colours());
    }
    if (not po.mutations_file.empty())
    {
        //TODO
        //ofstream mutations_ofs(po.mutations_file);
        //g.print_separated_het_mutations(mutations_ofs, 2, 20);
    }
    if (po.print_at_end)
    {
        cout << g.to_ptree();
    }
    if (not po.unmappable_contigs_file.empty())
    {
        //TODO
        //ofstream unmappable_contigs_ofs(po.unmappable_contigs_file);
        //g.print_unmappable_contigs(unmappable_contigs_ofs);
    }
    log_l(info) << ptree().put("main() done output", "");

    g.clear_and_dispose();
    log_l(info) << ptree().put("main() graph cleared", "");

    return EXIT_SUCCESS;
}

/*
    string input_file;
    string stats_file_1;
    string stats_file_2;
    string supercontig_lengths_file;
    string mutations_file;
    string unmappable_contigs_file;
    size_t progress_count;
    size_t unmap_trigger_len;
    bool merge_contigs_at_each_step;
    bool merge_contigs_at_end;
    bool print_graph;
    bool print_graph_each_step;
*/

int main(int argc, char* argv[])
{
    global::program_name = argv[0];

    Program_Options po;
    try
    {
        bo::options_description generic_opts_desc("Generic options");
        bo::options_description config_opts_desc("Configuration options");
        bo::options_description hidden_opts_desc("Hidden options");
        bo::options_description cmdline_opts_desc;
        bo::options_description visible_opts_desc("Allowed options");
        generic_opts_desc.add_options()
        ("help,?", "produce help message")
        ;
        config_opts_desc.add_options()
        ("input-file,i", bo::value< string >(&po.input_file), "input file")
        ("stats-file-1,x", bo::value< string >(&po.stats_file_1), "stats file 1")
        ("stats-file-2,y", bo::value< string >(&po.stats_file_2), "stats file 2")
        ("supercontig-lengths-file,l", bo::value< string >(&po.supercontig_lengths_file), "supercontig lengths file")
        ("mutations-file,M", bo::value< string >(&po.mutations_file), "mutations file")
        ("unmappable-contigs-file,U", bo::value< string >(&po.unmappable_contigs_file), "unmappable contigs file")
        ("progress-count,c", bo::value< size_t >(&po.progress_count)->default_value(0), "progress count")
        ("unmap-trigger-len,u", bo::value< size_t >(&po.unmap_trigger_len)->default_value(9), "unmap trigger len")
        ("cat-each-step,s", bo::bool_switch(&po.cat_at_step), "cat contigs at each step")
        ("cat-end,e", bo::bool_switch(&po.cat_at_end), "cat contigs at end")
        ("print-at-step,G", bo::bool_switch(&po.print_at_step), "print graph at each step")
        ("print-at-end,g", bo::bool_switch(&po.print_at_end), "print graph at end")
        ("default-log-level,d", bo::value< size_t >(&po.default_log_level)->default_value(0), "default log level")
        ;
        cmdline_opts_desc.add(generic_opts_desc).add(config_opts_desc).add(hidden_opts_desc);
        visible_opts_desc.add(generic_opts_desc).add(config_opts_desc);
        bo::variables_map vm;
        store(bo::command_line_parser(argc, argv).options(cmdline_opts_desc).run(), vm);
        notify(vm);

        if (vm.count("help"))
        {
            cout << visible_opts_desc;
            exit(EXIT_SUCCESS);
        }
        // set default log level
        logger::log::default_level() = logger::level(po.default_log_level);
        /*
        if (po.seed == 0)
        {
            po.seed = time(NULL);
        }
        */
    }
    catch(exception& e)
    {
        cout << e.what() << "\n";
        return EXIT_FAILURE;
    }
    return real_main(po);
}
