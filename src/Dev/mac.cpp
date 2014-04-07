#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <list>

#include "globals.hpp"
#include "global_assert.hpp"
#include "shortcuts.hpp"
#include "Graph.hpp"
#include "indent.hpp"
#include "ixstream.hpp"

using namespace std;
using namespace MAC;

auto inc_tab = indent::inc;
auto dec_tab = indent::dec;
using indent::nl;
using indent::tab;


int main(int argc, char* argv[])
{
    global::program_name = argv[0];
    char c;
    while ((c = getopt(argc, argv, "i:x:y:segGpc:u:l:M:U:")) != -1)
    {
        istringstream optarg_s(optarg != NULL? optarg : "");
        switch (c) {
            case 'i':
                global::input_file = optarg;
                break;
            case 'x':
                global::stats_file_1 = optarg;
                break;
            case 'y':
                global::stats_file_2 = optarg;
                break;
            case 's':
                global::merge_contigs_at_each_step = true;
                break;
            case 'e':
                global::merge_contigs_at_end = true;
                break;
            case 'g':
                global::print_graph = true;
                break;
            case 'G':
                global::print_graph_each_step = true;
                break;
            case 'p':
                global::progress_graph_op = true;
                break;
            case 'c':
                optarg_s >> global::progress_count;
                break;
            case 'u':
                optarg_s >> global::unmap_trigger_len;
                break;
            case 'l':
                global::supercontig_lengths_file = optarg;
                break;
            case 'M':
                global::mutations_file = optarg;
                break;
            case 'U':
                global::unmappable_contigs_file = optarg;
                break;
            default:
                cerr << "unrecognized option: " << c << endl;
                exit(EXIT_FAILURE);
        }
    }
    if (global::input_file.empty())
    {
        cerr << "no input file\n";
        exit(EXIT_FAILURE);
    }
    clog << "merging contigs at each step: " << (global::merge_contigs_at_each_step? "yes" : "no") << '\n';
    clog << "merging contigs at end: " << (global::merge_contigs_at_end? "yes" : "no") << '\n';
    clog << "using unmap_trigger_len: " << global::unmap_trigger_len << '\n';

    Graph g;
    ixstream ixs(global::input_file);
    if (not ixs)
    {
        cerr << "error opening file: " << argv[1] << '\n';
        return EXIT_FAILURE;
    }

    string line;
    size_t line_count = 0;
    while (getline(ixs, line))
    {
        if (global::progress_graph_op)
        {
            cerr << line << '\n';
        }
        istringstream iss(line + "\n");
        string rec_type;
        iss >> rec_type;
        if (rec_type == "HT")
        {
            // ignore header line for now
        }
        else if (rec_type == "VT")
        {
            string* r_id_ptr = new string();
            string* r_seq_ptr = new string();
            iss >> *r_id_ptr >> *r_seq_ptr;
            ASSERT(not iss.eof());
            global::assert_message = string("VT ") + *r_id_ptr;
            g.add_read(r_id_ptr, r_seq_ptr);
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
            //g.add_overlap(r1_id, r2_id, r1_start, r1_end - r1_start, r2_start, r2_end - r2_start, rc, sam_cigar.substr(5));
        }
        if (global::print_graph_each_step)
        {
            cerr << g;
        }
        if (global::progress_count > 0)
        {
            if ((++line_count % global::progress_count) == 0)
            {
                cerr << line_count << '\n';
            }
        }
        //ASSERT(g.check_all());
    }
    //g.unmap_single_chunks();
    //g.set_contig_ids();
    ASSERT(g.check_all());
    if (global::merge_contigs_at_end)
    {
        if (not global::stats_file_2.empty())
        {
            ofstream stats_2_ofs(global::stats_file_2);
            //g.dump_detailed_counts(stats_2_ofs);
        }
        //g.merge_all_read_contigs();
        ASSERT(g.check_all());
    }
    if (not global::stats_file_1.empty())
    {
        ofstream stats_1_ofs(global::stats_file_1);
        //g.dump_detailed_counts(stats_1_ofs);
    }
    if (not global::supercontig_lengths_file.empty())
    {
        ofstream lengths_file(global::supercontig_lengths_file);
        //g.print_supercontig_lengths_2(lengths_file);
        //ASSERT(g.check_colours());
    }
    if (not global::mutations_file.empty())
    {
        ofstream mutations_ofs(global::mutations_file);
        //g.print_separated_het_mutations(mutations_ofs, 2, 20);
    }
    if (global::print_graph)
    {
        cout << g;
    }
    if (not global::unmappable_contigs_file.empty())
    {
        ofstream unmappable_contigs_ofs(global::unmappable_contigs_file);
        //g.print_unmappable_contigs(unmappable_contigs_ofs);
    }
    cerr << "success\n";

    return EXIT_SUCCESS;
}
