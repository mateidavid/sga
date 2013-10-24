#include "MAC.hpp"

using namespace std;


namespace MAC
{
    void add_read(const string& name, const Seq_Type& seq, Read_Entry* re_ptr, Contig_Entry* ce_ptr)
    {
        // first, reset the read entry
        re_ptr->name = name;
        re_ptr->chunk_cont.clear();

        // reset contig entry
        ce_ptr->seq = seq;
        ce_ptr->mut_cont.clear();
        ce_ptr->chunk_ptr_cont.clear();

        // create new read chunk
        Read_Chunk c(re_ptr, ce_ptr, seq.size());
        re_ptr->chunk_cont.insert(c);

        // place read chunk into contig entry list
        const Read_Chunk* c_p = &(*re_ptr->chunk_cont.begin());
        ce_ptr->chunk_ptr_cont.insert(c_p);
    }
}
