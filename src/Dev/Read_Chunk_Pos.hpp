#ifndef __READ_CHUNK_POS_HPP
#define __READ_CHUNK_POS_HPP

#include "Read_Chunk.hpp"

namespace MAC
{

class Read_Chunk_Pos
{
public:
    Read_Chunk_Pos(Read_Chunk_CBPtr rc_cbptr, bool on_contig, bool start);

    Size_Type c_pos() const { return _c_pos_base; }
    Size_Type r_pos() const { return advance_r_pos(_r_pos_base, _mut_offset, true); }
    Size_Type match_len(bool i) const { return _match_len[i]; }
    bool in_mutation() const { return _mut_offset > 0; }
    Mutation_Cont::const_iterator mut_cit() const { return _mut_cit; }
    Allele_Idx_Cont::const_iterator al_idx_cit() const { return _al_idx_cit; }

private:
    Size_Type _c_pos_base;
    Size_Type _r_pos_base;
    Size_Type _mut_offset;

    Size_Type _match_len[2];
    Size_Type _mut_leftover;

    Size_Type _rc_c_start;
    Size_Type _rc_c_len;
    Mutation_Cont::const_iterator _mut_cit;
    Mutation_Cont::const_iterator _mut_cit_begin;
    Mutation_Cont::const_iterator _mut_cit_end;
    Allele_Idx_Cont::const_iterator _al_idx_cit;

    bool rc;

    static Size_Type advance_r_pos(Size_Type orig_pos, Size_Type len, bool c_forward)
    {
        return not rc == c_forward? orig_pos + len : orig_pos - len;
    }

    void set_convenience_fields()
    {
        if (mut_offset == 0)
        {
            _match_len[0] = c_pos - (_mut_cit != _mut_cit_begin? prev(_mut_cit)->rf_end() : _rc_c_start);
            _match_len[1] = (_mut_cit != _mut_cit_end? _mut_cit->rf_start() : _rc_c_end) - _c_pos;
            _mut_leftover = 0;
        }
        else // mut_offset > 0
        {
            _match_len[0] = 0;
            _match_len[1] = 0;
            _mut_leftover = mut_cit->allele_len(al_idx_cit->idx()) - _mut_offset;
        }
    }

    void increment(Size_Type r_brk)
    {
        check();
        ASSERT(not rc? r_pos() < r_brk : r_brk < r_pos());
        if (in_mutation() or match_len(1) == 0)
        {
            if (not rc? r_pos() + mut_leftover < r_brk : r_pos() - mut_leftover > r_brk)
            {
                // advance past variation
                _c_pos_base = mut_cit->rf_end();
                _r_pos_base = advance_r_pos(_r_pos_base, _mut_offset + mut_leftover, true);
                _mut_offset = 0;
                ++_mut_cit;
                ++_al_idx_cit;
            }
            else // r_brk <= _r_pos_base + mut_al_len (if not rc)
            {
                // stop inside (or at end of) variation
                _mut_offset = not rc? r_brk - _r_pos_base : _r_pos_base - r_brk;
            }
        }
        else // not in_mutation and match_len(1) > 0
        {
            if (not rc? r_pos() + match_len(1) <= r_brk : r_pos() - match_len(1) >= r_brk)
            {
                // advance past match
                _c_pos_base += match_len(1);
                _r_pos_base = advance_r_pos(_r_pos_base, match_len(1), true);
            }
            else // r_brk < r_pos() + match_len(1) (if not rc)
            {
                // stop inside match stretch
                _c_pos_base += (not rc? r_brk - r_pos() : r_pos() - r_brk);
                _r_pos_base = r_brk;
            }
        }
        set_convenience_fields();
        check();
    }

    void decrement(Size_Type r_brk)
    {
        check();
        ASSERT(not rc? r_brk < r_pos() : r_pos() < r_brk);
        if (in_mutation() or match_len(0) == 0)
        {
            if (not rc? r_brk < _r_pos_base : _r_pos_base < r_brk)
            {
                // advance past variation

            }
        }
        check();
    }


    void prev()
    {
        check();
        if (in_mut)
        {
            // crt chunk is mutation
            if (mut_start == 0)
            {
                // next chunk is non-mutation
                in_mut = false;
                ASSERT(c_start == mut_rf_start);
                mut_len = 0;
                c_len = c_start - (mut_cit != mut_cit_begin? std::prev(mut_cit)->rf_end() : rc_c_start);
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
            Size_Type leftover = c_start - (mut_cit != mut_cit_begin? std::prev(mut_cit)->rf_end() : rc_c_start);
            if (mut_cit != mut_cit_begin and leftover == 0)
            {
                // next chunk is mutation
                prev_mutation();
                mut_start = 0;
                mut_len = max(mut_rf_len, mut_al_len);
                in_mut = true;
                c_len = mut_rf_len;
                r_len = mut_al_len;
            }
            else // leftover > 0 or (leftover == 0 and mut_cit == mut_cit_begin)
            {
                // next chunk is non-mutation
                c_len = leftover;
                r_len = c_len;
            }
        }
        c_start -= c_len;
        r_start = (not rc? r_start - r_len : r_start + r_len);
        check();
    }

    void advance(bool forward) { if (forward) next(); else prev(); }

    friend bool operator == (const Read_Sub_Chunk& lhs, const Read_Sub_Chunk& rhs)
    {
        return lhs.c_start == rhs.c_start
            and lhs.c_len == rhs.c_len
            and lhs.r_start == rhs.r_start
            and lhs.r_len == rhs.r_len
            and lhs.mut_start == rhs.mut_start
            and lhs.mut_len == rhs.mut_len;
    }
    friend bool operator != (const Read_Sub_Chunk& lhs, const Read_Sub_Chunk& rhs) { return not(lhs == rhs); }

    void check() const
    {
#ifndef DISABLE_ASSERTS
        ASSERT(in_mut == (mut_start > 0 or mut_len > 0));
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
}; // class Read_Chunk_Pos

} // namespace MAC

#endif
