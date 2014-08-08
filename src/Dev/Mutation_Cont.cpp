//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Mutation_Cont.hpp"

using namespace std;


namespace MAC
{

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
            m_ptr->extend(tmp_m_ptr);
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

Mutation_CBPtr Mutation_Cont::find_span_pos(Size_Type c_pos) const
{
    auto it_range = Base::iintersect(c_pos, c_pos);
    for (auto it = it_range.begin(); it != it_range.end(); ++it)
    {
        if (it->get_start() < c_pos and c_pos < it->get_end())
        {
            return &*it;
        }
    }
    return nullptr;
}

Mutation_Cont Mutation_Cont::split(Size_Type c_brk, Mutation_CBPtr mut_left_cbptr)
{
    ASSERT(not mut_left_cbptr or (mut_left_cbptr->get_start() == c_brk and mut_left_cbptr->is_ins()));
    Mutation_Cont rhs_cont;
    auto it = begin();
    while (it != end())
    {
        Mutation_BPtr mut_bptr = &*it;
        ++it;
        // no Mutation may span c_brk
        ASSERT(not (mut_bptr->get_start() < c_brk and c_brk < mut_bptr->get_end()));
        // move the ones starting at or past the break, except possibly for mut_left_cbptr
        if (mut_bptr->get_start() >= c_brk and mut_bptr != mut_left_cbptr)
        {
            erase(mut_bptr);
            rhs_cont.insert(mut_bptr);
        }
    }
    return rhs_cont;
}

void Mutation_Cont::drop_unused()
{
    auto mut_it = begin();
    while (mut_it != end())
    {
        auto mut_next_it = mut_it;
        ++mut_next_it;
        Mutation_BPtr mut_bptr = &*mut_it;
        if (mut_bptr->chunk_ptr_cont().empty())
        {
            erase(mut_bptr);
            Mutation_Fact::del_elem(mut_bptr);
        }
        mut_it = mut_next_it;
    }
}

} // namespace MAC
