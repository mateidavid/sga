#include "MAC.hpp"

using namespace std;


namespace MAC
{
    void add_read(const string& name, const Seq_Type& seq, Read_Entry_Desc re_desc, Contig_Entry_Desc ce_desc)
    {
        // first, reset the read entry
        re_desc->name = name;
        re_desc->chunk_list.clear();

        // reset contig entry
        ce_desc->seq = seq;
        ce_desc->mut_list.clear();
        ce_desc->chunk_dlist.clear();

        // create new read chunk
        re_desc->chunk_list.emplace_back(re_desc, ce_desc, seq.size());

        // place read chunk into contig entry list
        ce_desc->chunk_dlist.push_back(&(*re_desc->chunk_list.begin()));
    }
}
