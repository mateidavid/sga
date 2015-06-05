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
    // main Bloom Filter
    BloomFilter* bf_p;
    // run stats
    atomic< long long > rf_kmers_total;
    atomic< long long > rf_kmers_used;
    atomic< long long > qr_kmers_total;
    atomic< long long > qr_kmers_used;
    atomic< long long > qr_kmers_in_rf;
    const string allowed_chars("ACGT");

    string description =
        "Construct a Bloom Filter using kmers from a reference read set, "
        "then trim the read ends from a query read set.\n"
        "Either -r or -l must be used. "
        "Either -q or -s must be used. "
        "-l and -s cannot be used together.";

    TCLAP::CmdLine cmd_parser(description, ' ', PACKAGE_VERSION);
    //
    // general parameters
    //
    TCLAP::MultiArg< string > log_level("d", "log-level", "Log level.", false, "string", cmd_parser);
    TCLAP::ValueArg< long > seed("", "seed", "Random seed (-1: use time).", false, -1, "int", cmd_parser);
    TCLAP::ValueArg< unsigned > progress("", "progress",   "Progress count.", false, 0, "int", cmd_parser);
    TCLAP::ValueArg< unsigned > num_threads("t", "threads",   "Number of threads.", false, 1, "int", cmd_parser);
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
        // avoid trimming reads to size 0: keep initial base only
        first_bp = 1;
        last_bp = 1;
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

int main(int argc, char* argv[])
{
    global::cmd_parser.parse(argc, argv);

    // validate command-line options
    if (global::rf_reads_fn.get().empty() == global::bf_load_fn.get().empty()
        or (not global::bf_load_fn.get().empty() and not global::bf_save_fn.get().empty())
        or (global::qr_reads_fn.get().empty() and global::bf_save_fn.get().empty()))
    {
        global::cmd_parser.usage();
        exit(EXIT_FAILURE);
    }

    // set log levels
    Logger::set_levels_from_options(global::log_level, &clog);
    // set random seed
    if (global::seed < 0)
    {
        global::seed.get() = time(nullptr);
    }
    srand48(global::seed);

    // print options
    LOG("main", info) << "program: " << global::cmd_parser.getProgramName() << endl;
    LOG("main", info) << "version: " << global::cmd_parser.getVersion() << endl;
    LOG("main", info) << "args: " << global::cmd_parser.getOrigArgv() << endl;

    if (not global::rf_reads_fn.get().empty())
    {
        // create Bloom Filter from ref reads
        unique_ptr< std::istream > rf_reads_is_p;
        if (global::rf_reads_fn.get() == "-")
        {
            rf_reads_is_p = unique_ptr< std::istream >(new zstr::istream(std::cin));
        }
        else
        {
            rf_reads_is_p = unique_ptr< std::istream >(new zstr::ifstream(global::rf_reads_fn));
        }
        build_rf_graph(*rf_reads_is_p);
        // save Bloom Filter
        if (not global::bf_save_fn.get().empty())
        {
            strict_fstream::fstream bf_save_os(global::bf_save_fn, ios_base::out | ios_base::binary);
            LOG("main", info) << "saving Bloom Filter to file: " << global::bf_save_fn.get() << endl;
            global::bf_p->save(bf_save_os);
            LOG("main", info) << "done saving Bloom Filter" << endl;
        }
    }
    else
    {
        // load BF from file
        strict_fstream::fstream bf_load_is(global::bf_load_fn, ios_base::in | ios_base::binary);
        global::bf_p = new BloomFilter();
        LOG("main", info) << "loading Bloom Filter from file: " << global::bf_load_fn.get() << endl;
        global::bf_p->load(bf_load_is);
        LOG("main", info) << "done loading Bloom Filter" << endl;
    }
    if (not global::qr_reads_fn.get().empty())
    {
        unique_ptr< std::istream > qr_reads_is_p;
        if (global::qr_reads_fn.get() == "-")
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
