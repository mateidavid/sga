//
// Construct a Bloom Filter from a set of (normal/G2) reads,
// then go through another set of (long/moleculo/G3) reads,
// and report for each of those,
// trim the ends which do not appear in the first set
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <time.h>
#include <atomic>

#include "version.h"
#include <tclap/CmdLine.h>

#include "logger.hpp"
#include "zstr.hpp"
#include "RC_Sequence.hpp"
#include "BloomFilter.h"
#include "pfor.hpp"
#include "kmer_gen.hpp"

using namespace std;

namespace global
{
    string program_name;

    // main Bloom Filter
    BloomFilter* bf_p;
    // run stats
    atomic< long long > rf_kmers_total;
    atomic< long long > rf_kmers_used;
    atomic< long long > qr_kmers_total;
    atomic< long long > qr_kmers_used;
    atomic< long long > qr_kmers_in_rf;
    const string allowed_chars("ACGT");

    TCLAP::CmdLine cmd_parser("Bloom Filter trimmer", ' ', PACKAGE_VERSION);
    //
    // general parameters
    //
    TCLAP::MultiArg< string > log_level("d", "log-level", "Log level.", false, "string", cmd_parser);
    TCLAP::ValueArg< long > seed("", "seed", "Random seed (0: use time).", false, -1, "int", cmd_parser);
    TCLAP::ValueArg< unsigned > num_threads("t", "threads",   "Number of threads.", false, 1, "int", cmd_parser);
    TCLAP::ValueArg< unsigned > progress("", "progress",   "Progress count.", false, 0, "int", cmd_parser);
    TCLAP::ValueArg< unsigned > chunk_size("", "chunk-size", "Thread chunk size.", false, 100, "int", cmd_parser);
    //
    // i/o-related parameters
    //
    TCLAP::ValueArg< string > rf_reads_fn("r", "rf-reads", "Reference reads.", false, "", "file", cmd_parser);
    TCLAP::ValueArg< string > qr_reads_fn("q", "qr-reads", "Query reads.", false, "", "file", cmd_parser);
    TCLAP::ValueArg< string > bf_load_fn("l", "bf-load", "BF to load.", false, "", "file", cmd_parser);
    TCLAP::ValueArg< string > bf_save_fn("s", "bf-save", "BF to save.", false, "", "file", cmd_parser);
    TCLAP::SwitchArg fasta_input("f", "fasta-input", "Reads in fasta input, not fastq.", cmd_parser, false);
    TCLAP::ValueArg< unsigned > phred_offset("", "phred-offset", "PHRED offset.", false, 33, "int", cmd_parser);
    //
    // core parameters
    //
    TCLAP::ValueArg< unsigned > k("k", "kmer-size", "Kmer size.", false, 61, "int", cmd_parser);
    TCLAP::ValueArg< unsigned > min_qv("m", "min-qv", "Minimum quality value.", false, 20, "int", cmd_parser);
    TCLAP::ValueArg< unsigned long > bf_size("", "bf-size", "BF size, in bits.", false, 20llu * 2llu * 3100000000llu, "int", cmd_parser);
    TCLAP::ValueArg< unsigned > bf_hashes("", "bf-hashes", "BF hashes.", false, 10, "int", cmd_parser);
}

// no frills sequence holder, along with istream reader
class SimpleRead
{
public:
    string id;
    rc_sequence::Sequence< string > sq[2];
    rc_sequence::Sequence< string > qv[2];

    bool get_from(istream& is, bool is_fasta = false)
    {
        id.clear();
        sq[0].clear();
        sq[1].clear();
        qv[0].clear();
        qv[1].clear();
        string tmp;
        if (not is_fasta)
        {
            return getline(is, id)
                and getline(is, sq[0])
                and getline(is, tmp)
                and getline(is, qv[0]);
        }
        else
        {
            if (not getline(is, id))
                return false;
            while (true)
            {
                int c = is.peek();
                if (c == '>' or c == std::char_traits<char>::eof())
                    break;
                if (not getline(is, tmp))
                    break;
                sq[0] += tmp;
            }
            return true;
        }
    }
}; // class SimpleRead

struct Thread_Output_Storage
{
    Thread_Output_Storage() {}
    ~Thread_Output_Storage()
    {
        cout << os.str() << flush;
        clog << es.str() << flush;
    }
    ostringstream os;
    ostringstream es;
}; // struct Thread_Output_Storage

bool is_valid_base(const string& sq, const string& qv, size_t i)
{
    return (global::allowed_chars.find(sq[i]) != string::npos
            and qv[i] - global::phred_offset >= global::min_qv);
}

void add_to_rf_graph(SimpleRead& r, Thread_Output_Storage& tos)
{
    if (r.sq[0].size() < (size_t)global::k)
        // no kmers here
        return;

    for (size_t i = 0; i < r.sq[0].size(); ++i)
        r.sq[0][i] = toupper(r.sq[0][i]);
    r.sq[1] = r.sq[0].revcomp();
    if (not global::fasta_input)
        r.qv[1] = r.qv[0].rev();

    global::rf_kmers_total += 2 * (r.sq[0].size() - global::k + 1);

    for (int st = 0; st < 2; ++st)
    {
        auto kg = kmer_gen(
            r.sq[st].size(), global::k,
            [&] (size_t i) { return is_valid_base(r.sq[st], r.qv[st], i); });

        for (auto it = kg.begin(); it != kg.end(); ++it)
        {
            // have full kmer
            LOG("main", debug1, tos.es) << "kmer: " << r.sq[0].substr(*it, global::k) << "\n";

            global::bf_p->add(&r.sq[st][*it], global::k);
            ++global::rf_kmers_used;
        }
    }
}

void build_rf_graph(istream& is)
{
    LOG("main", info)
        << "allocating Bloom Filter" << endl;
    global::bf_p = new BloomFilter(global::bf_size, global::bf_hashes);
    LOG("main", info)
        << "starting build threads" << endl;

    pfor< SimpleRead, Thread_Output_Storage >(
        nullptr,
        [&is] (SimpleRead& r)
        {
            return r.get_from(is, global::fasta_input);
        },
        &add_to_rf_graph,
        nullptr,
        global::num_threads,
        (global::fasta_input? 1 : global::chunk_size),
        1);

    LOG("main", info)
        << "ending build threads" << endl
        << "reference kmers total: " << global::rf_kmers_total << endl
        << "reference kmers used: " << global::rf_kmers_used << endl;
}

// Potentially add current interval to vector of results
//
// min_i : position (0-based) of 1st kmer that was ok
// i : position (0-based) of 1st kmer that is no longer ok
//
// The interval potentially added to v consists of 1-based positions
// *completely* spanned by (k) good kmers.
void add_int_to_res(vector<pair<size_t,size_t>>& v, size_t min_i, size_t i)
{
    if (min_i + global::k <= i)
        v.push_back(pair<size_t,size_t>(min_i + global::k - 1 + 1, i - 1 + 1));
}


void process_qr_read(SimpleRead& r, Thread_Output_Storage& tos)
{
    if (r.sq[0].size() < (size_t)global::k)
        // no kmers here
        return;

    vector< pair< size_t, size_t > > res;

    global::qr_kmers_total += r.sq[0].size() - global::k + 1;

    auto kg = kmer_gen(
        r.sq[0].size(), global::k,
        [&] (size_t i) { return is_valid_base(r.sq[0], r.qv[0], i); });

    size_t next_i = 0;
    size_t min_i = 0;
    for (auto it = kg.begin(); it != kg.end(); ++it)
    {
        // new kmer starting at i
        size_t i = *it;
        ++global::qr_kmers_used;

        if (i == next_i and global::bf_p->test(&r.sq[0][i], global::k))
        {
            LOG("main", debug1, tos.es) << "found kmer: " << r.sq[0].substr(i, global::k) << "\n";
            // contiguous with previous
            ++global::qr_kmers_in_rf;
            next_i = i + 1;
        }
        else
        {
            if (i == next_i)
            {
                LOG("main", debug1, tos.es) << "missing kmer: " << r.sq[0].substr(i, global::k) << "\n";
            }
            else
            {
                LOG("main", debug1, tos.es) << "skipping from next_i=" << next_i << " to i=" << i << "\n";
            }
            // jump
            add_int_to_res(res, min_i, next_i);
            min_i = i + 1;
            next_i = i + 1;
        }
    }
    add_int_to_res(res, min_i, next_i);

    size_t first_bp;
    size_t last_bp;
    if (res.size() == 0)
    {
        first_bp = r.sq[0].size() + 1;
        last_bp = r.sq[0].size();
    }
    else
    {
        first_bp = res[0].first - global::k + 1;
	last_bp = res[res.size() - 1].second + global::k - 1;
    }

    tos.os
        << r.id << "\n"
        << r.sq[0].substr(first_bp - 1, (last_bp + 1) - first_bp) << "\n"
        << "+" << r.sq[0].substr(0, first_bp - 1) << " "
        << r.qv[0].substr(0, first_bp - 1) << " "
        << r.sq[0].substr(last_bp) << " "
        << r.qv[0].substr(last_bp) << "\n"
        << r.qv[0].substr(first_bp - 1, (last_bp + 1) - first_bp) << "\n";

    if (Logger::get_facility_level("main") >= level_wrapper::debug1)
    {
        tos.es << r.id << "\t" << r.sq[0].size() << "\t";
        for (size_t l = 0; l < res.size(); ++l)
        {
            if (l > 0) tos.es << ",";
            tos.es << res[l].first - global::k + 1 << "-"
                   << res[l].second + global::k - 1;
        }
        tos.es << "\n";
    }
}

void process_qr_reads(istream& is)
{
    LOG("main", info)
        << "starting processing threads" << endl;

    pfor< SimpleRead, Thread_Output_Storage >
        (nullptr,
         [&is] (SimpleRead& r)
         {
             return r.get_from(is);
         },
         &process_qr_read,
         nullptr,
         global::num_threads,
         global::chunk_size,
         1);

    LOG("main", info)
        << "ending processing threads" << endl
        << "query kmers total: " << global::qr_kmers_total << endl
        << "query kmers used: " << global::qr_kmers_used << endl
        << "query kmers in reference: " << global::qr_kmers_in_rf << endl;
}

void usage(ostream & os)
{
    os << "Usage:" << endl
       << "  " << global::program_name << " OPTS -r <reference_reads> -q <query_reads>" << endl
       << "  " << global::program_name << " OPTS -r <reference_reads> -s <save_bf_file>" << endl
       << "  " << global::program_name << " OPTS -l <load_bf_file> -q <query_reads>" << endl
       << endl
       << "Synopsis:" << endl
       << "  Construct a Bloom Filter using kmers from a reference read set," << endl
       << "  then trim the read ends from a query read set." << endl
       << endl;
}

int main(int argc, char* argv[])
{
    global::program_name = argv[0];

    /*
    bo::options_description generic_opts_desc("Generic options");
    bo::options_description config_opts_desc("Configuration options");
    bo::options_description hidden_opts_desc("Hidden options");
    bo::options_description cmdline_opts_desc;
    bo::options_description visible_opts_desc("Allowed options");
    generic_opts_desc.add_options()
        ("help,?", "produce help message")
        // hack, see: http://lists.boost.org/boost-users/2010/01/55054.php
        ("log-level,d", bo::value(&global::log_level)->default_value(vector< string >(), ""), "log level")
        ("seed", bo::value(&global::seed)->default_value(time(nullptr), "use time"), "random seed")
        ("threads,t", bo::value(&global::num_threads)->default_value(1), "number of threads")
        ("progress", bo::value(&global::progress)->default_value(0), "progress count")
        ("chunk-size", bo::value(&global::chunk_size)->default_value(100), "progress count")
        ;
    config_opts_desc.add_options()
        //
        // file-related options
        //
        ("rf-reads,r", bo::value(&global::rf_reads_fn), "reference reads file")
        ("qr-reads,q", bo::value(&global::qr_reads_fn), "query reads file")
        ("bf-load,l", bo::value(&global::bf_load_fn), "load Bloom Filter from file")
        ("bf-save,s", bo::value(&global::bf_save_fn), "save Bloom Filter to file")
        ("fasta-input,f", bo::bool_switch(&global::fasta_input)->default_value(false), "input reads in fasta (not fastq) format")
        ("phred-offset", bo::value(&global::phred_offset)->default_value(33), "phred offset")
        //
        // general parameters
        //
        ("kmer-size,k", bo::value(&global::k)->default_value(61), "kmer size")
        ("min-qv,m", bo::value(&global::min_qv)->default_value(20), "minimum quality value")
        //
        // Bloom Filter parameters
        //
        ("bf-size", bo::value(&global::bf_size)->default_value(20ll * 2ll * 3100000000ll), "Bloom Filter size, in bits")
        ("bf-hashes", bo::value(&global::bf_hashes)->default_value(10), "Bloom Filter number of hashes")
        ;
    any_converter ac;
    ac.add_string_converter< string >();
    ac.add_string_converter< bool >();
    ac.add_string_converter< unsigned >();
    ac.add_string_converter< long >();
    ac.add_string_converter< size_t >();
    ac.add_converter(&cont_to_ptree< vector< string > >);
    cmdline_opts_desc.add(generic_opts_desc).add(config_opts_desc).add(hidden_opts_desc);
    visible_opts_desc.add(generic_opts_desc).add(config_opts_desc);
    bo::variables_map vm;
    try
    {
        store(bo::command_line_parser(argc, argv).options(cmdline_opts_desc).run(), vm);
        // if help requested, print it and stop
        if (vm.count("help"))
        {
            usage(cout);
            cout << visible_opts_desc;
            exit(EXIT_SUCCESS);
        }
        notify(vm);
    }
    catch(bo::error& e) 
    { 
        cerr << "ERROR: " << e.what() << endl << endl;
        usage(cerr);
        cerr << visible_opts_desc;
        exit(EXIT_FAILURE);
    }
    */

    global::cmd_parser.parse(argc, argv);

    // validate command-line options
    if (global::rf_reads_fn.getValue().empty() == global::bf_load_fn.getValue().empty()
        or (not global::bf_load_fn.getValue().empty() and not global::bf_save_fn.getValue().empty())
        or (global::qr_reads_fn.getValue().empty() and global::bf_save_fn.getValue().empty()))
    {
        usage(cerr);
        //cerr << visible_opts_desc;
        exit(EXIT_FAILURE);
    }

    // set log levels
    Logger::set_levels_from_options(global::log_level, &clog);
    // set random seed
    if (global::seed < 0)
    {
        global::seed.getValue() = time(nullptr);
    }
    srand48(global::seed);

    // print options
    LOG("main", info) << "program: " << global::cmd_parser.getProgramName() << endl;
    LOG("main", info) << "version: " << global::cmd_parser.getVersion() << endl;
    if (Logger::get_facility_level("main") <= level_wrapper::info)
    {
        ostringstream oss;
        for (int i = 1; i < argc; ++i)
        {
            oss << (i > 1? " " : "") << argv[i];
        }
        LOG("main", info) << "args: " << oss.str() << endl;
    }
    //LOG("main", info) << variables_map_converter::to_ptree(vm, ac);

    if (not global::rf_reads_fn.getValue().empty())
    {
        // create Bloom Filter from ref reads
        unique_ptr< std::istream > rf_reads_is_p;
        if (global::rf_reads_fn.getValue() == "-")
        {
            rf_reads_is_p = unique_ptr< std::istream >(new zstr::istream(std::cin));
        }
        else
        {
            rf_reads_is_p = unique_ptr< std::istream >(new zstr::ifstream(global::rf_reads_fn));
        }
        build_rf_graph(*rf_reads_is_p);
        // save Bloom Filter
        if (not global::bf_save_fn.getValue().empty())
        {
            strict_fstream::fstream bf_save_os(global::bf_save_fn, ios_base::out | ios_base::binary);
            LOG("main", info) << "saving Bloom Filter to file: " << global::bf_save_fn.getValue() << endl;
            global::bf_p->save(bf_save_os);
            LOG("main", info) << "done saving Bloom Filter" << endl;
        }
    }
    else
    {
        // load BF from file
        strict_fstream::fstream bf_load_is(global::bf_load_fn, ios_base::in | ios_base::binary);
        global::bf_p = new BloomFilter();
        LOG("main", info) << "loading Bloom Filter from file: " << global::bf_load_fn.getValue() << endl;
        global::bf_p->load(bf_load_is);
        LOG("main", info) << "done loading Bloom Filter" << endl;
    }
    if (not global::qr_reads_fn.getValue().empty())
    {
        unique_ptr< std::istream > qr_reads_is_p;
        if (global::qr_reads_fn.getValue() == "-")
        {
            qr_reads_is_p = unique_ptr< std::istream >(new zstr::istream(std::cin));
        }
        else
        {
            qr_reads_is_p = unique_ptr< std::istream >(new zstr::ifstream(global::qr_reads_fn));
        }
        process_qr_reads(*qr_reads_is_p);
    }
    if (global::bf_p) delete global::bf_p;
}
