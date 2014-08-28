#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <list>
#include <boost/program_options.hpp>

#include "Graph.hpp"
#include "ixstream.hpp"
#include "logger.hpp"

#define LOG_FACILITY "mac"


using namespace std;
using namespace MAC;
namespace bo = boost::program_options;


struct Program_Options
{
    string input_file;
    string stats_file;
    string supercontig_lengths_file;
    string mutations_file;
    string unmappable_contigs_file;
    string terminal_reads_file;
    string save_file;
    string load_file;
    size_t progress_count;
    size_t unmap_trigger_len;
    long seed;
    vector< string > log_level;
    bool cat_at_step;
    bool cat_at_end;
    bool unmap_read_ends;
    bool print_at_step;
    bool print_at_end;

    boost::property_tree::ptree to_ptree() const
    {
        return ptree().put("input file", input_file)
                      .put("load file", load_file)
                      .put("save file", save_file)
                      .put("stats file", stats_file)
                      .put("supercontig lengths file", supercontig_lengths_file)
                      .put("mutations file", mutations_file)
                      .put("unmappable contigs file", unmappable_contigs_file)
                      .put("terminal reads file", terminal_reads_file)
                      .put("progress count", progress_count)
                      .put("unmap trigger length", unmap_trigger_len)
                      .put("log levels", cont_to_ptree< vector<string>, string >(log_level, [] (const string& s) { return boost::property_tree::ptree(s); }))
                      .put("cat contigs at each step", cat_at_step)
                      .put("cat contigs at end", cat_at_end)
                      .put("unmap read ends", unmap_read_ends)
                      .put("print graph at each step", print_at_step)
                      .put("print graph at end", print_at_end)
                      .put("seed", seed)
                      ;
    }
};

void load_asqg(std::istream& is, const Program_Options& po, Graph& g)
{
    logger("io", info) << ptree("load_asqg_start");
    string line;
    size_t line_count = 0;
    while (getline(is, line))
    {
        logger(debug) << ptree("op").put("line", line);
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
        logger(debug2) << g.to_ptree();
        if (po.progress_count > 0)
        {
            if ((++line_count % po.progress_count) == 0)
            {
                logger(info) << ptree("progress").put("count", line_count);
            }
        }
        //ASSERT(g.check_all());
    }
    logger("io", info) << ptree("load_asqg_end");
    ASSERT(g.check_all());
}


int real_main(const Program_Options& po)
{
    Graph g;

    logger(info) << ptree("settings", po.to_ptree());

    if (po.input_file.empty() == po.load_file.empty())
    {
        cerr << "exactly 1 of input file or load file must be specified\n";
        exit(EXIT_FAILURE);
    }
    if (not po.input_file.empty())
    {
        ixstream ixs(po.input_file);
        if (not ixs)
        {
            cerr << "error opening file: " << po.input_file << '\n';
            return EXIT_FAILURE;
        }
        logger(info) << ptree("loading").put("input_file", po.input_file);
        load_asqg(ixs, po, g);
    }
    else
    {
        ifstream ifs(po.load_file, ios_base::in | ios_base::binary);
        if (not ifs)
        {
            cerr << "error opening file: " << po.load_file << '\n';
            return EXIT_FAILURE;
        }
        logger(info) << ptree("loading").put("load_file", po.load_file);
        g.load(ifs);
    }
    logger("mac", info) << ptree("factory_stats", g.factory_stats());

    if (po.unmap_read_ends)
    {
        g.unmap_read_ends();
        ASSERT(g.check_all());
    }
    if (po.cat_at_end)
    {
        g.cat_all_read_contigs();
        ASSERT(g.check_all());
    }
    if (not po.stats_file.empty())
    {
        ofstream stats_ofs(po.stats_file);
        g.dump_detailed_counts(stats_ofs);
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
    if (not po.terminal_reads_file.empty())
    {
        ofstream ofs(po.terminal_reads_file);
        if (not ofs)
        {
            cerr << "error opening file: " << po.terminal_reads_file << "\n";
            exit(EXIT_FAILURE);
        }
        g.get_scontig_terminal_reads(ofs);
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
    logger(info) << ptree("done_output");

    if (not po.save_file.empty())
    {
        ofstream ofs(po.save_file.c_str(), ios_base::out | ios_base::binary);
        if (not ofs)
        {
            cerr << "error opening file: " << po.save_file << '\n';
            return EXIT_FAILURE;
        }
        logger(info) << ptree("saving").put("save_file", po.save_file);
        g.save(ofs);
    }

    logger("mac", info) << ptree("factory_stats", g.factory_stats());
    g.clear_and_dispose();
    logger(info) << ptree("graph_cleared");

    return EXIT_SUCCESS;
}


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
        ("stats-file,x", bo::value< string >(&po.stats_file), "stats file")
        ("supercontig-lengths-file,l", bo::value< string >(&po.supercontig_lengths_file), "supercontig lengths file")
        ("mutations-file,M", bo::value< string >(&po.mutations_file), "mutations file")
        ("unmappable-contigs-file,U", bo::value< string >(&po.unmappable_contigs_file), "unmappable contigs file")
        ("terminal-reads", bo::value< string >(&po.terminal_reads_file), "terminal reads file")
        ("progress-count,c", bo::value< size_t >(&po.progress_count)->default_value(0), "progress count")
        ("unmap-trigger-len,u", bo::value< size_t >(&po.unmap_trigger_len)->default_value(9), "unmap trigger len")
        ("cat-each-step,s", bo::bool_switch(&po.cat_at_step), "cat contigs at each step")
        ("cat-end,e", bo::bool_switch(&po.cat_at_end), "cat contigs at end")
        ("unmap-read-ends", bo::bool_switch(&po.unmap_read_ends), "unmap read ends")
        ("print-at-step,G", bo::bool_switch(&po.print_at_step), "print graph at each step")
        ("print-at-end,g", bo::bool_switch(&po.print_at_end), "print graph at end")
        ("log-level,d", bo::value< vector< string > >(&po.log_level)->composing(), "default log level")
        ("save,S", bo::value< string >(&po.save_file), "save file")
        ("load,L", bo::value< string >(&po.load_file), "load file")
        ("seed", bo::value< long >(&po.seed), "RNG seed")
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

        // set log levels
        for (const auto& l : po.log_level)
        {
           size_t i = l.find(':');
           if (i == string::npos)
           {
              Logger::set_default_level(Logger::level_from_string(l));
              clog << "set default log level to: " << static_cast< int >(Logger::get_default_level()) << "\n";
           }
           else
           {
              Logger::set_facility_level(l.substr(0, i), Logger::level_from_string(l.substr(i + 1)));
              clog << "set log level of '" << l.substr(0, i) << "' to: "
                   << static_cast< int >(Logger::get_facility_level(l.substr(0, i))) << "\n";
           }
        }
        if (po.seed == 0)
        {
            po.seed = time(nullptr);
        }
    }
    catch(exception& e)
    {
        cout << e.what() << "\n";
        return EXIT_FAILURE;
    }
    return real_main(po);
}
