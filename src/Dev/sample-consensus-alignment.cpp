#include <iostream>
#include <string>
#include <vector>

#include <seqan/store.h>
#include <seqan/consensus.h>
#include <seqan/align.h>
#include <seqan/graph_msa.h>

using namespace std;

void do_alignment(const vector< string >& v)
{
    seqan::FragmentStore<> store;
    seqan::AlignedReadLayout layout;
    seqan::ConsensusAlignmentOptions options;
    //SimpleScore score;

    // Resize contigStore and contigNameStore (required for printing the first layout).
    //resize(store.contigStore, 1);
    //appendValue(store.contigNameStore, "ref");

    // Append reads (includes small errors).
    //appendRead(store, "AATGGATGGCAAAATAGTTGTTCCATGAATACATCTCTAAAGAGCTTT");
    //appendRead(store, "AAAGTAGTTGTTCCATGAATACATCTCTAAAGAGCTTTGATGCTAATTT");
    //appendRead(store, "AGTTGTCCATGAATACATCTCTAAAGAGCTTTGATGCTAATTTAGTCAATTTTCAATACTGTA");
    //appendRead(store, "ACATCTCTTAAAGAGCTTTGATGCTAATTTAGTCAAATT");
    //appendRead(store, "AGAGCTTTGATGCTAATTTAGTCAAATTTTCAATACTGTACAATCTTCTCTAG");

    // Actual read layout.
    //
    // AATGGATGGCAAAATAGTTGTTCCATGAATACATCTCTAAAGAGCTTT
    //           AAAATAGTTGTTCCATGAATACATCTCTAAAGAGCTTTGATGCTAATTT
    //                AGTTGTTCCATGAATACATCTCTAAAGAGCTTTGATGCTAATTTAGTCAAATTTTCAATACTGTA
    //                               ACATCTCTAAAGAGCTTTGATGCTAATTTAGTCAAATT
    //                                         AGAGCTTTGATGCTAATTTAGTCAAATTTTCAATACTGTACAATCTTCTCTAG

    // The position used in the following are only approximate and would
    // not lead to the read layout above.
    size_t max_read_len = 0;
    cout << "Reads:" << endl;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        cout << v[i] << endl;
        seqan::appendRead(store, v[i].data());
        seqan::appendAlignedRead(store, i, 0, 0, 0 + (int)length(store.readSeqStore[i]));
        max_read_len = max(max_read_len, v[i].size());
    }
    //convertMatchesToGlobalAlignment(store, score, true);

    // Print the (wrong) alignment.
    //cout << "Initial alignment" << endl << endl;
    //layoutAlignment(layout, store);
    //printAlignment(cout, layout, store, /*contigID=*/ 0, /*beginPos=*/ 0, /*endPos=*/ 100, 0, 30);

    options.useContigID = true;
    options.usePositions = true;
    seqan::consensusAlignment(store, options);

    cout << endl << "Consensus alignment:" << endl;
    seqan::layoutAlignment(layout, store);
    seqan::printAlignment(cout, layout, store, /*contigID=*/ 0, /*beginPos=*/ 0, /*endPos=*/ 100, 0, 30);

    cout << endl << "Contig:" << endl << store.contigStore[0].seq << endl;

    seqan::reAlignment(store, 0, 1, 2*max_read_len, false);
    cout << endl << "Consensus after realignment:" << endl;
    seqan::layoutAlignment(layout, store);
    seqan::printAlignment(cout, layout, store, /*contigID=*/ 0, /*beginPos=*/ 0, /*endPos=*/ 100, 0, 30);

    cout << endl << "Contig:" << endl << store.contigStore[0].seq << endl;

}

void do_msa(const vector< string >& v)
{
    seqan::Align< seqan::DnaString > align;
    resize(rows(align), v.size());
    for (unsigned i = 0; i < v.size(); ++i)
    {
        seqan::assignSource(row(align, i), v[i]);
    }
    seqan::globalMsaAlignment(align, seqan::SimpleScore(5, -3, -1, -3));
    cout << "Align:" << endl << align << endl;
}

int main()
{
    vector< string > v;
    string s;
    while (getline(cin, s))
    {
        v.emplace_back(move(s));
    }
    if (v.empty())
    {
        v = { "AGAGCTTTTGATGCTAATTTAGTGAAAT",
              "AGAGCTTTGATGCTAATTAGTGAAAT",
              "AGAGCTTTGATCCTAATTTAGTGAAAT",
              "AGAGCTTTGATCCTAATTTAGTGAAAT",
              "AGAGCTTTGATGCTAATTTAGTGAAAT" };
    }
    do_alignment(v);
    do_msa(v);
}
