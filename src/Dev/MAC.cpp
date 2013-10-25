#include "MAC.hpp"

using namespace std;
using boost::tie;


namespace MAC
{
    void add_read(Read_Entry_Cont& re_cont, Contig_Entry_Cont& ce_cont,
                  const string* name_ptr, const Seq_Type* seq_ptr)
    {
        // first, create read entry & contig entry
        Read_Entry re(name_ptr, seq_ptr->size());
        Contig_Entry ce(seq_ptr);
        ce.add_chunk(re.get_cptr_first_chunk());

        // insert them in their containers
        Read_Entry_Cont::iterator re_it;
        bool success;
        tie(re_it, success) = re_cont.insert(re);
        assert(success);
        Contig_Entry_Cont::iterator ce_it;
        tie(ce_it, success) = ce_cont.insert(ce_cont.end(), ce);
        assert(success);

        // set the contig entry pointer inside the single read chunk
        auto rc_modifier = [&] (Read_Chunk& rc) { rc.assign_to_contig(&(*ce_it), 0, seq_ptr->size(), false, vector<const Mutation*>()); };
        auto re_modifier = [&] (Read_Entry& re) { re.modify_read_chunk(re_it->get_cptr_first_chunk(), rc_modifier); };
        re_cont.modify(re_it, re_modifier);
    }

    void add_overlap(Read_Entry_Cont& re_cont, Contig_Entry_Cont& ce_cont,
                     const string& r1_name, const string& r2_name,
                     Size_Type r1_start, Size_Type r1_len,
                     Size_Type r2_start, Size_Type r2_len, bool r2_rc,
                     const string& cigar)
    {
        (void)ce_cont;
        (void)r1_start;
        (void)r1_len;
        (void)r2_start;
        (void)r2_len;
        (void)r2_rc;
        (void)cigar;

        Read_Entry_Cont::iterator re1_it = re_cont.find(r1_name);
        assert(re1_it != re_cont.end());
        Read_Entry_Cont::iterator re2_it = re_cont.find(r2_name);
        assert(re2_it != re_cont.end());
    }
}
