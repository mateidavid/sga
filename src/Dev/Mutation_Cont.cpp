//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "Mutation_Cont.hpp"


namespace MAC
{

Mutation_Cont::Mutation_Cont(const Cigar& cigar, const Seq_Proxy_Type& qr)
{
    Mutation_BPtr accum_mut_bptr;
    for (size_t i = 0; i <= cigar.n_ops(); ++i)
    {
        if (i == cigar.n_ops() or cigar.op_type(i) == '=')
        {
            // if we're keeping track of a mutation, add it to container
            if (accum_mut_bptr)
            {
                insert(accum_mut_bptr);
                accum_mut_bptr = nullptr;
            }
        }
        else
        {
            // disallow ambigous 'M' operation
            ASSERT(cigar.op_type(i) != 'M');
            // accumulator is either empty, or it extends to the start of this mutation
            ASSERT(not accum_mut_bptr or accum_mut_bptr->rf_end() == cigar.op_rf_pos(i));
            if (not accum_mut_bptr)
            {
                accum_mut_bptr = Mutation_Fact::new_elem();
            }
            if (qr.size() > 0)
            {
                // accept either the matched part of the query (length = cigar.qr_len)
                // or entire query (length > cigar.qr_len)
                ASSERT(qr.size() >= cigar.qr_len());
                Size_Type qr_offset = (not cigar.reversed()?
                                       cigar.op_qr_pos(i)
                                       : cigar.op_qr_pos(i) - cigar.op_qr_len(i));
                if (qr.size() == cigar.qr_len())
                {
                    // assume the query sequence before qr_start is missing
                    qr_offset -= cigar.qr_start();
                }
                accum_mut_bptr->extend(
                    cigar.op_rf_pos(i),
                    cigar.op_rf_len(i),
                    /*
                    (not cigar.reversed()?
                     qr.substr(qr_offset, cigar.op_qr_len(i))
                     : reverseComplement(qr.substr(qr_offset, cigar.op_qr_len(i)))));
                    */
                    qr.substr(qr_offset, cigar.op_qr_len(i)).revcomp(cigar.reversed()));
            }
            else
            {
                accum_mut_bptr->extend(
                    cigar.op_rf_pos(i),
                    cigar.op_rf_len(i),
                    cigar.op_qr_len(i));
            }
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
        if (it->rf_start() < c_pos and c_pos < it->rf_end())
        {
            return &*it;
        }
    }
    return nullptr;
}

Mutation_Cont Mutation_Cont::split(Size_Type c_brk, Mutation_CBPtr mut_left_cbptr)
{
    ASSERT(not mut_left_cbptr or (mut_left_cbptr->rf_start() == c_brk and mut_left_cbptr->is_ins()));
    Mutation_Cont rhs_cont;
    auto it = begin();
    while (it != end())
    {
        Mutation_BPtr mut_bptr = &*it;
        ++it;
        // no Mutation may span c_brk
        ASSERT(not (mut_bptr->rf_start() < c_brk and c_brk < mut_bptr->rf_end()));
        // move the ones starting at or past the break, except possibly for mut_left_cbptr
        if (mut_bptr->rf_start() >= c_brk and mut_bptr != mut_left_cbptr)
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
