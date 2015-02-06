//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

// pfor
//   General implementation of a parallel for loop.
//
//   Given a number of threads and a chunk size, the total work is split into
//   chunks, each chunk is processed independently by a thread, and the outputs
//   are reordered and sent to stdout in the same order they were read.
//
//
// TLS_pfor
//   Thread local storage during the loop.
//
// void do_before(TLS&)
// void do_after(TLS&)
//   Functions to execute before and after parallel section.
//
// bool get_item(TLS&, T&)
//   Read one item; return false if no item was retrieved. Before the call, the
//   item placeholder is in the state left after the previous loop (or after
//   default constructor).
//
// void process_item(TLS&, T&)
//   Process item, possibly using out&log streams in TLS.
//
//
// Sample call: (using lambda expressions, gcc >= 4.6)
//   pfor<int, TLS_pfor>(NULL,
//                       [&cin] (TLS_pfor& tls, int& e) { return bool(cin >> e); },
//                       [] (TLS_pfor& tls, int& e) { *tls.out_str_p << e << "\n"; },
//                       NULL,
//                       4,
//                       10);
//
//
// NOTE: might produce 'unused parameter (TLS)' warnings. Disable with
// -Wno-unused-parameter.


#ifndef __PFOR_HPP
#define __PFOR_HPP

#include <iostream>
#include <sstream>
#include <functional>
#include <vector>
#include <algorithm>
#include <queue>
#include <omp.h>
#include <assert.h>


// Thread Local Storage during pfor loop
class TLS_pfor
{
public:
    int tid;                    // thread id
    size_t cid;                 // chunk id
    std::ostream* out_str_p;    // output stream
    std::ostream* log_str_p;    // log stream
};

// Output Stream Chunk (used internally)
class OSC_pfor
{
public:
    int tid;
    size_t cid;
    std::ostringstream* out_str_p;
    std::ostringstream* log_str_p;
};

class OSC_pfor_comparator
{
public:
    bool operator () (const OSC_pfor& lhs, const OSC_pfor& rhs)
    {
        return lhs.cid > rhs.cid;
    }
};


template <typename T, typename TLS>
void pfor(std::function<void(TLS&)> do_before,
          std::function<bool(TLS&, T&)> get_item,
          std::function<void(TLS&, T&)> process_item,
          std::function<void(TLS&)> do_after,
          int num_threads,
          size_t chunk_size,
          size_t chunk_progress = 1)
{
    using namespace std;

    priority_queue<OSC_pfor, vector<OSC_pfor>, OSC_pfor_comparator> h;
    size_t cid_in = 0;
    size_t cid_out = 0;
    bool done = false;

#pragma omp parallel num_threads(num_threads)
    {
        TLS tls;
        tls.tid = omp_get_thread_num();
        vector<T> buff(chunk_size);
        size_t load = 0;

        if (do_before) do_before(tls);

        do
        {
            load = 0;

#pragma omp critical (pfor_input)
            {
                if (not done)
                {
                    while (load < chunk_size and get_item(tls, buff[load]))
                        ++load;
                    done = (load < chunk_size);
                    if (load > 0)
                    {
                        tls.cid = cid_in++;
                        if (cid_in % chunk_progress == 0)
                        {
#pragma omp critical (pfor_out)
                            {
                                clog << (cid_in - 1) * chunk_size + load << "..." << endl;
                            }
                        }
                    }
                }
            }
// end omp critical (pfor_input)

            if (load == 0)
                break;

            // real work
            tls.out_str_p = new ostringstream();
            tls.log_str_p = new ostringstream();
            for (size_t i = 0; i < load; ++i)
                process_item(tls, buff[i]);

            // output part
            {
                OSC_pfor osc;
                osc.tid = tls.tid;
                osc.cid = tls.cid;
                osc.out_str_p = static_cast<ostringstream*>(tls.out_str_p);
                osc.log_str_p = static_cast<ostringstream*>(tls.log_str_p);

#pragma omp critical (pfor_output)
                {
                    h.push(osc);
                    while (h.size() > 0)
                    {
                        osc = h.top();
                        assert(osc.cid >= cid_out);
                        if (osc.cid > cid_out)
                            break;
#pragma omp critical (pfor_out)
                        {
                            //clog << "cid:" << osc.cid << " work_thread:" << osc.tid << " output_thread:" << tls.cid << endl;
                            cout << osc.out_str_p->str();
                            cout.flush();
                        }
                        delete osc.out_str_p;
#pragma omp critical (pfor_out)
                        {
                            clog << osc.log_str_p->str();
                            clog.flush();
                        }
                        delete osc.log_str_p;
                        h.pop();
                        ++cid_out;
                    }
                }
// end omp critical (pfor_output)
            }

        } while (not done);

        if (do_after) do_after(tls);
    }
// end omp parallel
}


#endif
