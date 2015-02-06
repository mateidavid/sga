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
#include <omp.h>
#include <time.h>

#include <boost/program_options.hpp>
#include "variables_map_converter.hpp"

#include "logger.hpp"
#include "ixstream.hpp"
#include "RC_Sequence.hpp"
#include "BloomFilter.h"
#include "pfor-omp.hpp"
#include "kmer_gen.hpp"

using namespace std;
namespace bo = boost::program_options;


namespace global
{
    string program_name;
    bo::variables_map vm;
    unsigned k;
    BloomFilter* bf_p;
    long long rf_kmers_total;
    long long rf_kmers_used;
    long long qr_kmers_total;
    long long qr_kmers_used;
    long long qr_kmers_in_rf;
    const string allowed_chars("ACGT");
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
};

bool is_valid_base(const string& sq, const string& qv, size_t i)
{
    return (global::allowed_chars.find(sq[i]) != string::npos
            and qv[i] - global::vm.at("phred-offset").as< unsigned >() >= global::vm.at("min-qv").as< unsigned >());
}

void add_to_rf_graph(TLS_pfor& tls, SimpleRead& r)
{
    if (r.sq[0].size() < (size_t)global::k)
        // no kmers here
        return;

    for (size_t i = 0; i < r.sq[0].size(); ++i)
        r.sq[0][i] = toupper(r.sq[0][i]);
    r.sq[1] = r.sq[0].revcomp();
    if (not global::vm.at("fasta-input").as< bool >())
        r.qv[1] = r.qv[0].rev();

#pragma omp atomic
    global::rf_kmers_total += 2 * (r.sq[0].size() - global::k + 1);

    for (int st = 0; st < 2; ++st)
    {
        auto kg = kmer_gen(
            r.sq[st].size(), global::k,
            [&] (size_t i) { return is_valid_base(r.sq[st], r.qv[st], i); });

        for (auto it = kg.begin(); it != kg.end(); ++it)
        {
            // have full kmer
            if (Logger::get_facility_level("main") >= level_wrapper::debug1)
                *tls.log_str_p << "kmer: " << r.sq[0].substr(*it, global::k) << "\n";

            global::bf_p->add(&r.sq[st][*it], global::k);
#pragma omp atomic
            ++global::rf_kmers_used;
        }
    }
}

void build_rf_graph(istream& is)
{
    LOG("main", info)
        << "allocating Bloom Filter" << endl;
    global::bf_p = new BloomFilter(global::vm.at("bf-size").as< size_t >(), global::vm.at("bf-hashes").as< unsigned >());

    LOG("main", info)
        << "starting build threads" << endl;
    pfor<SimpleRead,TLS_pfor>(NULL,
                              [&is] (TLS_pfor& tls, SimpleRead& r)
                              {
                                  (void)tls;
                                  return r.get_from(is, global::vm.at("fasta-input").as< bool >());
                              },
                              &add_to_rf_graph,
                              NULL,
                              global::vm.at("threads").as< unsigned >(),
                              (global::vm.at("fasta-input").as< bool >()? 1 : global::vm.at("chunk-size").as< unsigned >()),
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


void process_qr_read(TLS_pfor& tls, SimpleRead& r)
{
    if (r.sq[0].size() < (size_t)global::k)
        // no kmers here
        return;

    vector<pair<size_t,size_t>> res;

#pragma omp atomic
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
#pragma omp atomic
        ++global::qr_kmers_used;


        if (Logger::get_facility_level("main") >= level_wrapper::debug2)
            *tls.log_str_p << r.sq[0].substr(i, global::k) << "\t";

        if (i == next_i and global::bf_p->test(&r.sq[0][i], global::k))
        {
            // contiguous with previous
            if (Logger::get_facility_level("main") >= level_wrapper::debug2)
                *tls.log_str_p << "yes\n";

#pragma omp atomic
            ++global::qr_kmers_in_rf;

            next_i = i + 1;
        }
        else
        {
            // jump
            if (Logger::get_facility_level("main") >= level_wrapper::debug2)
                *tls.log_str_p << "no\n";

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

    *(tls.out_str_p)
      << r.id << "\n"
      << r.sq[0].substr(first_bp - 1, (last_bp + 1) - first_bp) << "\n"
      << "+" << r.sq[0].substr(0, first_bp - 1) << " "
      << r.qv[0].substr(0, first_bp - 1) << " "
      << r.sq[0].substr(last_bp) << " "
      << r.qv[0].substr(last_bp) << "\n"
      << r.qv[0].substr(first_bp - 1, (last_bp + 1) - first_bp) << "\n";

    if (Logger::get_facility_level("main") >= level_wrapper::debug1)
    {
        *tls.log_str_p << r.id << "\t" << r.sq[0].size() << "\t";
        for (size_t l = 0; l < res.size(); ++l)
        {
            if (l > 0) *tls.out_str_p << ",";
            *tls.log_str_p << res[l].first - global::k + 1 << "-"
                           << res[l].second + global::k - 1;
        }
        *tls.log_str_p << "\n";
    }
}

void process_qr_reads(istream& is)
{
    LOG("main", info)
        << "starting processing threads" << endl;

    pfor<SimpleRead,TLS_pfor>(NULL,
                              [&is] (TLS_pfor& tls, SimpleRead& r)
                              {
                                  (void)tls;
                                  return r.get_from(is);
                              },
                              &process_qr_read,
                              NULL,
                              global::vm.at("threads").as< unsigned >(),
                              global::vm.at("chunk-size").as< unsigned >(),
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

    bo::options_description generic_opts_desc("Generic options");
    bo::options_description config_opts_desc("Configuration options");
    bo::options_description hidden_opts_desc("Hidden options");
    bo::options_description cmdline_opts_desc;
    bo::options_description visible_opts_desc("Allowed options");
    generic_opts_desc.add_options()
        ("help,?", "produce help message")
        // hack, see: http://lists.boost.org/boost-users/2010/01/55054.php
        ("log-level,d", bo::value< vector< string > >()->default_value(vector< string >(), ""), "log level")
        ("threads,t", bo::value< unsigned >()->default_value(1), "number of threads")
        ("seed", bo::value< unsigned >()->default_value(time(nullptr), "use time"), "random seed")
        ("progress", bo::value< unsigned >()->default_value(0), "progress count")
        ("chunk-size", bo::value< unsigned >()->default_value(100), "progress count")
        ;
    config_opts_desc.add_options()
        //
        // file-related options
        //
        ("rf-reads,r", bo::value< string >()->default_value(""), "reference reads file")
        ("qr-reads,q", bo::value< string >()->default_value(""), "query reads file")
        ("bf-load,l", bo::value< string >()->default_value(""), "load Bloom Filter from file")
        ("bf-save,s", bo::value< string >()->default_value(""), "save Bloom Filter to file")
        ("fasta-input,f", bo::bool_switch()->default_value(false), "input reads in fasta (not fastq) format")
        ("phred-offset", bo::value< unsigned >()->default_value(33), "phred offset")
        //
        // general parameters
        //
        ("kmer-size,k", bo::value< unsigned >()->default_value(61), "kmer size")
        ("min-qv,m", bo::value< unsigned >()->default_value(20), "minimum quality value")
        //
        // Bloom Filter parameters
        //
        ("bf-size", bo::value< size_t >()->default_value(20ll * 2ll * 3100000000ll), "size of Bloom Filter in bits")
        ("bf-hashes", bo::value< unsigned >()->default_value(10), "number of hashes in Bloom Filter")
        ;
    any_converter ac;
    ac.add_string_converter< string >();
    ac.add_string_converter< bool >();
    ac.add_string_converter< unsigned >();
    ac.add_string_converter< size_t >();
    ac.add_converter(&cont_to_ptree< vector< string > >);
    cmdline_opts_desc.add(generic_opts_desc).add(config_opts_desc).add(hidden_opts_desc);
    visible_opts_desc.add(generic_opts_desc).add(config_opts_desc);
    try
    {
        store(bo::command_line_parser(argc, argv).options(cmdline_opts_desc).run(), global::vm);
        // if help requested, print it and stop
        if (global::vm.count("help"))
        {
            usage(cout);
            cout << visible_opts_desc;
            exit(EXIT_SUCCESS);
        }
        notify(global::vm);
    }
    catch(bo::error& e) 
    { 
        cerr << "ERROR: " << e.what() << endl << endl;
        usage(cerr);
        cerr << visible_opts_desc;
        exit(EXIT_FAILURE);
    }

    // validate command-line options
    const string& rf_reads_fn = global::vm.at("rf-reads").as< string >();
    const string& qr_reads_fn = global::vm.at("qr-reads").as< string >();
    const string& bf_load_fn = global::vm.at("bf-load").as< string >();
    const string& bf_save_fn = global::vm.at("bf-save").as< string >();
    if (rf_reads_fn.empty() == bf_load_fn.empty()
        or (not bf_load_fn.empty() and not bf_save_fn.empty())
        or (qr_reads_fn.empty() and bf_save_fn.empty()))
    {
        usage(cerr);
        cerr << visible_opts_desc;
        exit(EXIT_FAILURE);
    }

    // set log levels
    for (const auto& l : global::vm.at("log-level").as< vector< string > >())
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
    srand48(global::vm.at("seed").as< unsigned >());

    // print options
    LOG("main", info) << variables_map_converter::to_ptree(global::vm, ac);

    // save k as global variable for easier access
    global::k = global::vm.at("kmer-size").as< unsigned >();

    if (not rf_reads_fn.empty())
    {
        // create Bloom Filter from ref reads
        ixstream rf_reads_is(rf_reads_fn);
        build_rf_graph(rf_reads_is);
        // save Bloom Filter
        if (not bf_save_fn.empty())
        {
            fstr bf_save_os(bf_save_fn, ios_base::out | ios_base::binary);
            LOG("main", info) << "saving Bloom Filter to file: " << bf_save_fn << endl;
            global::bf_p->save(bf_save_os);
            LOG("main", info) << "done saving Bloom Filter" << endl;
        }
    }
    else
    {
        // load BF from file
        fstr bf_load_is(bf_load_fn, ios_base::in | ios_base::binary);
        global::bf_p = new BloomFilter();
        LOG("main", info) << "loading Bloom Filter from file: " << bf_load_fn << endl;
        global::bf_p->load(bf_load_is);
        LOG("main", info) << "done loading Bloom Filter" << endl;
    }
    if (not qr_reads_fn.empty())
    {
        ixstream qr_reads_is(qr_reads_fn);
        process_qr_reads(qr_reads_is);
    }
}
