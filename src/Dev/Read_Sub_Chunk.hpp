#ifndef __READ_SUB_CHUNK_HPP
#define __READ_SUB_CHUNK_HPP

#include "Read_Chunk.hpp"

namespace MAC
{

struct Read_Sub_Chunk
{
    Size_Type c_start;
    Size_Type c_len;
    Size_Type r_start;
    Size_Type r_len;
    Size_Type mut_start;
    Size_Type mut_len;
    Size_Type rc_c_start;
    Size_Type rc_c_len;
    Mutation_Cont::const_iterator mut_cit;
    Mutation_Cont::const_iterator mut_cit_begin;
    Mutation_Cont::const_iterator mut_cit_end;
    Allele_Idx_Cont::const_iterator al_idx_cit;

    Size_Type mut_rf_start;
    Size_Type mut_rf_len;
    Size_Type mut_al_len;
    unsigned al_idx;

    bool rc;

    Size_Type rc_c_end() const { return rc_c_start + rc_c_len; }
    Size_Type c_end() const { return c_start + c_len; }

    void check() const
    {
#ifndef DISABLE_ASSERTS
        ASSERT(rc_c_start <= c_start);
        ASSERT(c_end() <= rc_c_end());
        if (mut_cit == mut_cit_end
            or c_start < mut_cit->rf_start())
        {
            ASSERT(mut_start == 0);
            ASSERT(mut_len == 0);
            ASSERT(c_len == r_len);
        }
        else // mut_cit->rf_start() <= c_start
        {
            ASSERT(al_idx == al_idx_cit->idx());
            ASSERT(mut_rf_start == mut_cit->rf_start());
            ASSERT(mut_rf_len == mut_cit->rf_len());
            ASSERT(mut_al_len == mut_cit->allele_len(al_idx));
            ASSERT(c_end() <= mut_cit->rf_end());
            ASSERT(mut_start <= max(mut_rf_len, mut_al_len));
            ASSERT(c_start == mut_rf_start + min(mut_rf_len, mut_start));
            if (mut_rf_len <= mut_al_len)
            {
                ASSERT(r_len == mut_len);
                ASSERT(c_len == min(mut_start + mut_len, mut_rf_len) - min(mut_start, mut_rf_len));
            }
            else
            {
                ASSERT(c_len == mut_len);
                ASSERT(r_len == min(mut_start + mut_len, mut_al_len) - min(mut_start, mut_al_len));
            }
        }
#endif
    }

    void set_convenience_fields()
    {
        if (mut_cit != mut_cit_end)
        {
            al_idx = al_idx_cit->idx();
            mut_rf_start = mut_cit->rf_start();
            mut_rf_len = mut_cit->rf_end();
            mut_al_len = mut_cit->allele_len(al_idx);
        }
        else
        {
            al_idx = 0;
            mut_rf_start = 0;
            mut_rf_len = 0;
            mut_al_len = 0;
        }
    }

    void next_mutation()
    {
        ASSERT(mut_cit != mut_cit_end);
        ++mut_cit;
        ++al_idx_cit;
        set_convenience_fields();
        check();
    }

    void prev_mutation()
    {
        ASSERT(mut_cit != mut_cit_begin);
        --mut_cit;
        --al_idx_cit;
        set_convenience_fields();
        check();
    }

    void next()
    {
        check();
        ASSERT(c_end() < rc_c_end());
        c_start += c_len;
        r_start = (not rc? r_start + r_len : r_start - r_len);
        if (mut_start > 0 or mut_len > 0)
        {
            // crt chunk is mutation
            mut_start += mut_len;
            leftover_mut_len = max(mut_rf_len, mut_al_len) - mut_start;
            if (leftover_mut_len == 0)
            {
                // next chunk is non-mutation
                ASSERT(c_start == mut_rf_start + mut_rf_len);
                mut_start = 0;
                mut_len = 0;
                next_mutation();
                c_len = (mut_cit != mut_cit_end? mut_rf_start : rc_c_end) - c_start;
                r_len = c_len;
            }
            else
            {
                // next chunk is same mutation
                mut_len = leftover_mut_len;
                c_len = min(mut_start + mut_len, mut_rf_len) - min(mut_start, mut_rf_len);
                r_len = min(mut_start + mut_len, mut_al_len) - min(mut_start, mut_al_len);
            }
        }
        else // mut_start == 0 and mut_len == 0
        {
            // crt chunk is non-mutation
            if (c_start == mut_rf_start)
            {
                // next chunk is mutation
                c_len = mut_rf_len;
                r_len = mut_al_len;
                mut_len = max(c_len, r_len);
            }
            else
            {
                // next chunk is non-mutation
                ASSERT(c_start < mut_rf_start);
                c_len = mut_rf_start - c_start;
                r_len = c_len;
            }
        }
    }

    void prev()
    {
        check();
        ASSERT(rc_c_start < c_start);
        if (mut_start > 0 or mut_len > 0)
        {
            // crt chunk is mutation
            if (mut_start == 0)
            {
                // next chunk is non-mutation
                ASSERT(c_start == mut_rf_start);
                mut_len = 0;
                c_len = c_start - (mut_cit != mut_cit_begin? prev(mut_cit)->rf_end() : rc_c_start);
                r_len = c_len;
            }
            else
            {
                // next chunk is same mutation
                mut_len = mut_start;
                mut_start = 0;
                c_len = min(mut_len, mut_rf_len);
                r_len = min(mut_len, mut_al_len);
            }
        }
        else // mut_start == 0 and mut_len == 0
        {
            // crt chunk is non-mutation
            Size_Type leftover = c_start - (mut_cit != mut_cit_begin? prev(mut_cit)->rf_end() : rc_c_start);
            if (leftover == 0)
            {
                // next chunk is mutation
                prev_mutation();
                mut_start = 0;
                mut_len = max(mut_rf_len, mut_al_len);
                c_len = mut_rf_len;
                r_len = mut_al_len;
            }
            else // leftover > 0
            {
                // next chunk is non-mutation
                c_len = leftover;
                r_len = c_len;
            }
        }
        c_start -= c_len;
        r_start = (not rc? r_start - r_len : r_start + r_len);
    }
}; // class Read_Sub_Chunk

} // nameespace MAC


#endif
