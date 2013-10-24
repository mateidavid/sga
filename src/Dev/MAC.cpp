#include "MAC.hpp"

using namespace std;
using boost::tie;


namespace MAC
{
    void add_read(const string* name_ptr, const Seq_Type* seq_ptr, Read_Entry_Cont& re_cont, Contig_Entry_Cont& ce_cont)
    {
        // first, create read entry & contig entry
        Read_Entry re(name_ptr, seq_ptr);
        Contig_Entry ce(seq_ptr);
        ce.add_chunk(re.get_ptr_first_chunk());

        // insert them in their containers
        Read_Entry_Cont::iterator re_it;
        bool success;
        tie(re_it, success) = re_cont.insert(re);
        ce_cont.push_back(ce);

        // set the contig entry pointer inside the single read chunk
        auto f = bind(&Read_Entry::set_ce_ptr, placeholders::_1, re_it->get_ptr_first_chunk(), &ce_cont.back());
        re_cont.modify(re_it, f);
    }
}
