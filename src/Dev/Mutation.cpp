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

    ASSERT(base_offset <= _len);
    ASSERT(alt_offset <= _seq_len);

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
    ASSERT(rf.size() == _len);
    if (not have_seq())
    {
        return;
    }
    while (_len > 0 and _seq_len > 0 and rf[_len - 1] == _seq[_seq_len - 1])
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

bool operator == (const Mutation& lhs, const Mutation& rhs)
{
    return (lhs._start == rhs._start
            and lhs._len == rhs._len
            and lhs._seq_len == rhs._seq_len
            and lhs.have_seq() == rhs.have_seq()
            and (not lhs.have_seq() or lhs._seq == rhs._seq));
}

Mutation_Cont::Mutation_Cont(const Cigar& cigar, const string& qr)
{
    Mutation_BPtr m_ptr = Mutation_Fact::new_elem(cigar.get_rf_start(), 0);
    for (size_t i = 0; i <= cigar.get_n_ops(); ++i)
    {
        ASSERT(m_ptr);
        if (i == cigar.get_n_ops() or cigar.get_op(i) == '=')
        {
            // construct mutation and add it to container
            if (m_ptr->get_len() > 0 or m_ptr->get_seq_len() > 0)
            {
                insert(m_ptr);
            }
            else
            {
                Mutation_Fact::del_elem(m_ptr);
            }
            m_ptr = Mutation_Fact::new_elem(cigar.get_rf_offset(i) + cigar.get_rf_op_len(i), 0); // ok to call even with i == n_ops
        }
        else
        {
            // disallow ambigous 'M' operation
            ASSERT(cigar.get_op(i) != 'M');
            Mutation_BPtr tmp_m_ptr;
            if (qr.size() > 0)
            {
                // accept either the matched part of the query (length = cigar.qr_len)
                // or entire query (length > cigar.qr_len)
                ASSERT(qr.size() >= cigar.get_qr_len());
                Size_Type qr_offset = (not cigar.is_reversed()? cigar.get_qr_offset(i) : cigar.get_qr_offset(i) - cigar.get_qr_op_len(i));
                if (qr.size() == cigar.get_qr_len())
                {
                    // assume the query sequence before qr_start is missing
                    qr_offset -= cigar.get_qr_start();
                }
                tmp_m_ptr = Mutation_Fact::new_elem(cigar.get_rf_offset(i), cigar.get_rf_op_len(i),
                                                    (not cigar.is_reversed()?
                                                     qr.substr(qr_offset, cigar.get_qr_op_len(i))
                                                     : reverseComplement(qr.substr(qr_offset, cigar.get_qr_op_len(i)))));
            }
            else
            {
                tmp_m_ptr = Mutation_Fact::new_elem(cigar.get_rf_offset(i), cigar.get_rf_op_len(i), cigar.get_qr_op_len(i));
            }
            m_ptr->merge(std::move(static_cast< Mutation& >(*tmp_m_ptr)));
            Mutation_Fact::del_elem(tmp_m_ptr);
        }
    }
}

/*
Mutation_BPtr Mutation_Cont::add_mut(Mutation_BPtr mut_bptr)
{
    Base::iterator it;
    Base::iterator it_end;
    for (std::tie(it, it_end) = this->equal_range(*mut_bptr); it != it_end; ++it)
    {
        if (*mut_bptr == *it)
        {
            return &*it;
        }
    }
    this->insert(*mut_bptr);
    return mut_bptr;
}
*/

ostream& operator << (ostream& os, const Mutation& rhs)
{
    os << "(Mutation &=" << (void*)&rhs << ", start=" << (size_t)rhs.get_start() << ",len=" << (size_t)rhs.get_len()
       << ",seq_len=" << rhs.get_seq_len() << ",seq=" << rhs.get_seq() << ")";
    return os;
}

ostream& operator << (ostream& os, const Mutation_Ptr_Node& rhs)
{
    os << "(Mutation_Ptr_Node &=" << (void*)&rhs << ", _m_bptr=" << rhs._m_bptr << ")";
    return os;
}

} // namespace MAC
