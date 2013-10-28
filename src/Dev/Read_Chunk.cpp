//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Read_Chunk.hpp"

#include <boost/tuple/tuple.hpp>

#include "Mutation.hpp"
#include "print_seq.hpp"

using namespace std;
using boost::tie;


namespace MAC
{
    /** Create a set of mutations to a reference string based on a cigar object.
     * @param rf_start Index in reference of the start of the match (0-based).
     * @param cigar Cigar object describing the match.
     * @return Container of Mutation objects.
     */
    Mutation_Cont make_mutations_from_cigar(Size_Type rf_start, const Cigar& cigar)
    {
        Mutation_Cont res;
        for (size_t i = 0; i < cigar.get_n_ops(); ++i)
        {
            if (not Cigar::is_match_op(cigar.get_op(i)))
            {
                Mutation_Cont::iterator it;
                bool success;
                tie(it, success) = res.insert(Mutation(rf_start + cigar.get_rf_offset(i), cigar.get_rf_op_len(i), cigar.get_qr_op_len(i)));
                assert(success);
            }
        }
        return res;
    }

    vector< pair< Read_Chunk, Mutation_Cont > > Read_Chunk::make_chunks_from_cigar(
        Size_Type r1_start, Size_Type r2_start, const Cigar& cigar)
    {
        // create objects with default constructor
        vector< pair< Read_Chunk, Mutation_Cont > > res;

        // fix lengths and rc flags
        Size_Type r1_len = cigar.get_rf_len();
        Size_Type r2_len = cigar.get_qr_len();
        res[0].first._r_start = r2_start;
        res[0].first._r_len = r2_len;
        res[0].first._c_start = r1_start;
        res[0].first._c_len = r1_len;
        res[0].first._rc = cigar.is_reversed();
        res[1].first._r_start = r1_start;
        res[1].first._r_len = r1_len;
        res[1].first._c_start = r2_start;
        res[1].first._c_len = r2_len;
        res[1].first._rc = cigar.is_reversed();

        // construct r1->r2 mutations
        res[0].second = make_mutations_from_cigar(r1_start, cigar);
        for (auto it = res[0].second.begin(); it != res[0].second.end(); ++it)
            res[0].first._mut_ptr_cont.push_back(&(*it));

        // construct r2->r1 mutations
        res[1].second = make_mutations_from_cigar(r2_start, cigar.complement());
        for (auto it = res[1].second.begin(); it != res[1].second.end(); ++it)
            res[1].first._mut_ptr_cont.push_back(&(*it));

        return res;
    }

    ostream& operator << (ostream& os, const Read_Chunk& rhs)
    {
        os << "(r_start=" << rhs.get_r_start() << ",r_len=" << rhs.get_r_len()
        << ",c_start=" << rhs.get_c_start() << ",c_len=" << rhs.get_c_len() << ",rc=" << (int)rhs.get_rc() << ",mut_list=(";
        print_ptr_seq(os, rhs._mut_ptr_cont);
        os << "))";
        return os;
    }
}
