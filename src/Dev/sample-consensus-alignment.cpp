#include <iostream>

#include <seqan/store.h>
#include <seqan/consensus.h>

using namespace std;
using namespace seqan;

int main()
{
    FragmentStore<> store;
    AlignedReadLayout layout;
    ConsensusAlignmentOptions options;
    //SimpleScore score;

    // Resize contigStore and contigNameStore (required for printing the first layout).
    //resize(store.contigStore, 1);
    //appendValue(store.contigNameStore, "ref");

    // Actual read layout.
    //
    // AATGGATGGCAAAATAGTTGTTCCATGAATACATCTCTAAAGAGCTTT
    //           AAAATAGTTGTTCCATGAATACATCTCTAAAGAGCTTTGATGCTAATTT
    //                AGTTGTTCCATGAATACATCTCTAAAGAGCTTTGATGCTAATTTAGTCAAATTTTCAATACTGTA
    //                               ACATCTCTAAAGAGCTTTGATGCTAATTTAGTCAAATT
    //                                         AGAGCTTTGATGCTAATTTAGTCAAATTTTCAATACTGTACAATCTTCTCTAG

    // Append reads (includes small errors).
    //appendRead(store, "AATGGATGGCAAAATAGTTGTTCCATGAATACATCTCTAAAGAGCTTT");
    //appendRead(store, "AAAGTAGTTGTTCCATGAATACATCTCTAAAGAGCTTTGATGCTAATTT");
    //appendRead(store, "AGTTGTCCATGAATACATCTCTAAAGAGCTTTGATGCTAATTTAGTCAATTTTCAATACTGTA");
    //appendRead(store, "ACATCTCTTAAAGAGCTTTGATGCTAATTTAGTCAAATT");
    //appendRead(store, "AGAGCTTTGATGCTAATTTAGTCAAATTTTCAATACTGTACAATCTTCTCTAG");

    appendRead(store, "AGAGCTTTTGATGCTAATTTAGTGAAAT");
    appendRead(store, "AGAGCTTTGATGCTAATTAGTGAAAT");
    appendRead(store, "AGAGCTTTGATCCTAATTTAGTGAAAT");
    appendRead(store, "AGAGCTTTGATCCTAATTTAGTGAAAT");
    appendRead(store, "AGAGCTTTGATGCTAATTTAGTGAAAT");

    // The position used in the following are only approximate and would
    // not lead to the read layout above.
    appendAlignedRead(store, 0, 0, 0, 0 + (int)length(store.readSeqStore[0]));
    appendAlignedRead(store, 1, 0, 0, 0 + (int)length(store.readSeqStore[1]));
    appendAlignedRead(store, 2, 0, 0, 0 + (int)length(store.readSeqStore[2]));
    appendAlignedRead(store, 3, 0, 0, 0 + (int)length(store.readSeqStore[3]));
    appendAlignedRead(store, 4, 0, 0, 0 + (int)length(store.readSeqStore[4]));
    //convertMatchesToGlobalAlignment(store, score, true);

    // Print the (wrong) alignment.
    //cout << "Initial alignment" << endl << endl;
    //layoutAlignment(layout, store);
    //printAlignment(cout, layout, store, /*contigID=*/ 0, /*beginPos=*/ 0, /*endPos=*/ 100, 0, 30);

    options.useContigID = true;
    consensusAlignment(store, options);

    cout << "Consensus alignment:" << endl << endl;
    layoutAlignment(layout, store);
    printAlignment(cout, layout, store, /*contigID=*/ 0, /*beginPos=*/ 0, /*endPos=*/ 100, 0, 30);

    cout << "Contig:" << endl
         << store.contigStore[0].seq << endl;

    return 0;
}
