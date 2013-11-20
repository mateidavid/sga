//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Mutation.hpp"

#include <iostream>
#include "../Util/Util.h"
#include "indent.hpp"

using namespace std;


namespace MAC
{
    Mutation Mutation::cut(Size_Type base_offset, Size_Type alt_offset)
    {
        Mutation res;

        assert(base_offset <= _len);
        assert(alt_offset <= _seq_len);

        if (have_seq())
        {
            res = Mutation(_start + base_offset, _len - base_offset, _seq.substr(alt_offset));
            _len = base_offset;
            _seq_len = alt_offset;
            _seq.erase(alt_offset);
        }
        else
        {
            res = Mutation(_start + base_offset, _len - base_offset, _seq_len - alt_offset);
            _len = base_offset;
            _seq_len = alt_offset;
        }
        return res;
    }

    shared_ptr< Mutation_Cont > make_mutations_from_cigar(const Cigar& cigar, const string& qr)
    {
        shared_ptr< Mutation_Cont > res(new Mutation_Cont());
        Mutation m;
        for (size_t i = 0; i <= cigar.get_n_ops(); ++i)
        {
            if (i == cigar.get_n_ops() or cigar.get_op(i) == '=')
            {
                // construct mutation and add it to container
                if (m.get_len() > 0 or m.get_seq_len() > 0)
                {
                    Mutation_Cont::iterator it;
                    bool success;
                    tie(it, success) = res->insert(m);
                    assert(success);
                }
                m = Mutation(cigar.get_rf_offset(i) + cigar.get_rf_op_len(i), 0); // ok to call even with i == n_ops
            }
            else
            {
                // disallow ambigous 'M' operation
                assert(cigar.get_op(i) != 'M');
                Mutation tmp;
                if (qr.size() > 0)
                {
                    // accept either the matched part of the query (length = cigar.qr_len)
                    // or entire query (length > cigar.qr_len)
                    assert(qr.size() >= cigar.get_qr_len());
                    Size_Type qr_offset = (not cigar.is_reversed()? cigar.get_qr_offset(i) : cigar.get_qr_offset(i) - cigar.get_qr_op_len(i));
                    if (qr.size() == cigar.get_qr_len())
                    {
                        // assume the query sequence before qr_start is missing
                        qr_offset -= cigar.get_qr_start();
                    }
                    tmp = Mutation(cigar.get_rf_offset(i), cigar.get_rf_op_len(i),
                                   (not cigar.is_reversed()?
                                    qr.substr(qr_offset, cigar.get_qr_op_len(i))
                                    : reverseComplement(qr.substr(qr_offset, cigar.get_qr_op_len(i)))));
                }
                else
                {
                    tmp = Mutation(cigar.get_rf_offset(i), cigar.get_rf_op_len(i), cigar.get_qr_op_len(i));
                }
                m.merge(tmp);
            }
        }
        return res;
    }

    ostream& operator << (ostream& os, const Mutation& rhs)
    {
        os << "(Mutation &=" << (void*)&rhs
           << indent::inc << indent::nl << "start=" << (size_t)rhs.get_start()
           << ",len=" << (size_t)rhs.get_len()
           << ",seq_len=" << rhs.get_seq_len()
           << ",seq=" << rhs.get_seq()
           << indent::dec << indent::nl << ")";
        return os;
    }
}
