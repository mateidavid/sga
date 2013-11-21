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

    void Mutation::simplify(const Seq_Type& rf)
    {
        assert(rf.size() == _len);
        if (not have_seq())
            return;
        while (_len > 0 and _seq_len > 0 and rf[_len - 1] == _seq[_len - 1])
        {
            --_len;
            --_seq_len;
            _seq.resize(_seq.size() - 1);
        }
        size_t i = 0;
        while (i < _len and i < _seq_len and rf[i] == _seq[i])
        {
            ++i;
        }
        if (i > 0)
        {
            _start += i;
            _len -= i;
            _seq_len -= i;
            _seq = _seq.substr(i);
        }
    }

    bool Mutation::operator == (const Mutation& rhs) const
    {
        return (_start == rhs._start
                and _len == rhs._len
                and _seq_len == rhs._seq_len
                and have_seq() == rhs.have_seq()
                and (not have_seq() or _seq == rhs._seq));
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

    Mutation_CPtr add_mut_to_cont(Mutation_Cont& mut_cont, const Mutation& mut)
    {
        Mutation_Cont::iterator it;
        Mutation_Cont::iterator it_end;
        for (std::tie(it, it_end) = mut_cont.equal_range(mut.get_key()); it != it_end; ++it)
        {
            if (mut == *it)
            {
                return &*it;
            }
        }
        bool success;
        tie(it, success) = mut_cont.insert(mut);
        assert(success);
        return &*it;
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
