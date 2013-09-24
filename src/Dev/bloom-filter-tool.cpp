//
// Construct a deBruijn graph from a set of (normal) reads,
// then go through another set of (long/moleculo) reads,
// and report for each of those,
// trim the ends which do not appear in the first set
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <omp.h>
#include <time.h>

#include "ixstream.hpp"
#include "Util.h"
#include "BloomFilter.h"
#include "pfor.hpp"
#include "kmer_gen.hpp"

using namespace std;


namespace global
{
    string prog_name;
    int verbosity;
    int num_threads = 1;
    size_t chunk_size = 100;

    BloomFilter* bf_p;
    int k = 61;
    size_t bf_size = 20ll * 2ll * 3100000000ll;
    int bf_hashes = 10;
    int min_qv = 20;
    int phred_offset = 33;
    long long progress = 100000;
    long seed;
    bool fasta_input = false;

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
    string sq[2];
    string qv[2];

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

void add_to_rf_graph(TLS_pfor& tls, SimpleRead& r)
{
    if (r.sq[0].size() < (size_t)global::k)
        // no kmers here
        return;

    for (size_t i = 0; i < r.sq[0].size(); ++i)
        r.sq[0][i] = toupper(r.sq[0][i]);
    r.sq[1] = reverseComplementIUPAC(r.sq[0]);
    if (not global::fasta_input)
        r.qv[1] = reverse(r.qv[0]);

#pragma omp atomic
    global::rf_kmers_total += 2 * (r.sq[0].size() - global::k + 1);

    for (int st = 0; st < 2; ++st)
    {
        Kmer_gen g(&r.sq[st], global::k, &global::allowed_chars,
                   &r.qv[st], global::min_qv, global::phred_offset);
        size_t i = 0;
        while ((i = g.get_next()) != string::npos)
        {
            // have full kmer
            if (global::verbosity > 2)
                *tls.log_str_p << "kmer: " << r.sq[0].substr(i, global::k) << "\n";

            global::bf_p->add(&r.sq[st][i], global::k);
#pragma omp atomic
            ++global::rf_kmers_used;
        }
    }
}

void build_rf_graph(istream& is)
{
    if (global::verbosity > 0) clog << "allocating Bloom Filter\n";
    global::bf_p = new BloomFilter(global::bf_size, global::bf_hashes);

    if (global::verbosity > 0) clog << "starting build threads\n";
    pfor<SimpleRead,TLS_pfor>(NULL,
                              [&is] (TLS_pfor& tls, SimpleRead& r)
                              {
                                  return r.get_from(is, global::fasta_input);
                              },
                              &add_to_rf_graph,
                              NULL,
                              global::num_threads,
                              (global::fasta_input? 1 : global::chunk_size),
                              1);
    if (global::verbosity > 0) clog << "ending build threads\n";

    if (global::verbosity > 0)
        clog << "reference kmers total: " << global::rf_kmers_total
             << "\nreference kmers used: " << global::rf_kmers_used << "\n";
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

    Kmer_gen g(&r.sq[0], global::k, &global::allowed_chars,
               &r.qv[0], global::min_qv, global::phred_offset);

    size_t i = 0;
    size_t next_i = 0;
    size_t min_i = 0;
    while ((i = g.get_next()) != string::npos)
    {
        // new kmer starting at i
#pragma omp atomic
        ++global::qr_kmers_used;

        if (global::verbosity > 2)
            *tls.log_str_p << r.sq[0].substr(i, global::k) << "\t";

        if (i == next_i and global::bf_p->test(&r.sq[0][i], global::k))
        {
            // contiguous with previous
            if (global::verbosity > 2)
                *tls.log_str_p << "yes\n";

#pragma omp atomic
            ++global::qr_kmers_in_rf;

            next_i = i + 1;
        }
        else
        {
            // jump
            if (global::verbosity > 2)
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

    if (global::verbosity > 1)
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
    if (global::verbosity > 0) clog << "starting processing threads\n";

    pfor<SimpleRead,TLS_pfor>(NULL,
                              [&is] (TLS_pfor& tls, SimpleRead& r)
                              {
                                  return r.get_from(is);
                              },
                              &process_qr_read,
                              NULL,
                              global::num_threads,
                              global::chunk_size,
                              1);

    if (global::verbosity > 0) clog << "ending processing threads\n";

    clog << "query kmers total: " << global::qr_kmers_total
         << "\nquery kmers used: " << global::qr_kmers_used
         << "\nquery kmers in reference: " << global::qr_kmers_in_rf << "\n";
}


void usage(ostream & os)
{
    os << "basic usage: " << global::prog_name << " -r <reference_reads> -q <query_reads>\n";
}


int main(int argc, char * argv[])
{
    global::prog_name = argv[0];

    string rf_reads_file;
    string qr_reads_file;
    string load_bf_file;
    string save_bf_file;

    char c;
    while ((c = getopt(argc, argv, "fr:q:l:s:k:m:z:n:S:vh")) != -1)
    {
        istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
        case 'f':
            global::fasta_input = true;
            break;
        case 'r':
            rf_reads_file = optarg;
            break;
        case 'q':
            qr_reads_file = optarg;
            break;
        case 'l':
            load_bf_file = optarg;
            break;
        case 's':
            save_bf_file = optarg;
            break;
        case 'k':
            arg >> global::k;
            break;
        case 'm':
            arg >> global::min_qv;
            break;
        case 'z':
            arg >> global::bf_size;
            break;
        case 'n':
            arg >> global::num_threads;
            break;
        case 'S':
            arg >> global::seed;
            break;
        case 'v':
            global::verbosity++;
            break;
        case 'h':
            usage(cout);
            exit(EXIT_SUCCESS);
        default:
            cerr << "unrecognized option: " << c << "\n";
            usage(cerr);
            exit(EXIT_FAILURE);
        }
    }
    if (optind != argc)
    {
        usage(cerr);
        exit(EXIT_FAILURE);
    }

    if ((rf_reads_file == "") == (load_bf_file == ""))
    {
        cerr << "Exactly one of reference reads file and BloomFilter file must be given\n";
        exit(EXIT_FAILURE);
    }

    if (load_bf_file != "" and save_bf_file != "")
    {
        cerr << "When loading a BloomFilter file, save cannot be given\n";
        exit(EXIT_FAILURE);
    }

    if (qr_reads_file == "" and save_bf_file == "")
    {
        cerr << "Either a query reads file or a BloomFilter save file must be given\n";
        exit(EXIT_FAILURE);
    }

    while (global::seed == 0)
    {
        global::seed = long(time(NULL));
    }
    srand(global::seed);

    if (global::verbosity > 0)
    {
        clog << "num_threads = " << global::num_threads << "\n";
        clog << "seed = " << global::seed << "\n";
        clog << "k = " << global::k << "\n";
        clog << "min_qv = " << global::min_qv << "\n";
        clog << "bf_size = " << global::bf_size << "\n";
        clog << "bf_hashes = " << global::bf_hashes << "\n";
    }

    ixstream* rf_reads_is_p = NULL;
    ixstream* qr_reads_is_p = NULL;
    ifstream* load_bf_is_p = NULL;
    ofstream* save_bf_os_p = NULL;
    if (rf_reads_file != "") rf_reads_is_p = new ixstream(rf_reads_file);
    if (qr_reads_file != "") qr_reads_is_p = new ixstream(qr_reads_file);
    if (load_bf_file != "") load_bf_is_p = new ifstream(load_bf_file);
    if (save_bf_file != "") save_bf_os_p = new ofstream(save_bf_file);

    if (rf_reads_file != "")
    {
        // create deBruijn graph from ref reads
        build_rf_graph(*rf_reads_is_p);

        if (save_bf_file != "")
        {
            global::bf_p->save(*save_bf_os_p);
        }
    }
    else
    {
        // load BF from file
        global::bf_p = new BloomFilter();
        global::bf_p->load(*load_bf_is_p);
        clog << "loaded Bloom Filter\n";
    }

    if (qr_reads_file != "")
    {
        process_qr_reads(*qr_reads_is_p);
    }

    return 0;
}
