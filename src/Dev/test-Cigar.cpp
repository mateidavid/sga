#include <iostream>

#ifdef DISABLE_ASSERTS
#undef DISABLE_ASSERTS
#endif

#include "global_assert.hpp"
#include "Cigar.hpp"
#include "RC_Sequence.hpp"


using namespace std;
typedef cigar::Cigar< unsigned > cigar_type;

void check_cigar(const string& cigar_string, const string& rf_seq, const string& qr_seq)
{
    clog << "checking: cigar_string=[" << cigar_string << "]" << endl;
    for (int st = 0; st < 2; ++st)
    {
        cigar_type c(cigar_string, st == 1);
        c.check(rf_seq, qr_seq);
        c.disambiguate(rf_seq, qr_seq);
        c.check(rf_seq, qr_seq);
        // check complement
        {
            auto c_comp = c.complement();
            c_comp.check(qr_seq, rf_seq);
            auto c_comp_comp = c_comp.complement();
            c_comp_comp.check(rf_seq, qr_seq);
            ASSERT(c == c_comp_comp);
        }
        // check substring operations
        {
            auto c_tmp = c.subcigar(0, 0);
            c_tmp.check(rf_seq, qr_seq);
            c_tmp = c.subcigar(0, c.n_ops());
            c_tmp.check(rf_seq, qr_seq);
            ASSERT(c == c_tmp);
        }
        // check cutting operations
        {
            for (size_t i = 0; i < c.n_ops(); ++i)
            {
                for (size_t j = 0; j <= c.op_len(i); ++j)
                {
                    cigar_type c_tmp(cigar_string, st == 1);
                    c_tmp.disambiguate(rf_seq, qr_seq);
                    c_tmp.cut_op(i, j);
                    c_tmp.check(rf_seq, qr_seq);
                }
            }
        }
        // check trim_end operations
        {
            cigar_type c_tmp;
            if (c.rf_len() > 0)
            {
                c_tmp = cigar_type(cigar_string, st == 1);
                c_tmp.disambiguate(rf_seq, qr_seq);
                c_tmp.trim_end(true, false, c.rf_start() + 1);
                c_tmp.check(rf_seq, qr_seq);
                c_tmp = cigar_type(cigar_string, st == 1);
                c_tmp.disambiguate(rf_seq, qr_seq);
                c_tmp.trim_end(true, true, c.rf_end() - 1);
                c_tmp.check(rf_seq, qr_seq);
            }
            if (c.qr_len() > 0)
            {
                c_tmp = cigar_type(cigar_string, st == 1);
                c_tmp.disambiguate(rf_seq, qr_seq);
                c_tmp.trim_end(false, false, c.qr_start() + 1);
                c_tmp.check(rf_seq, qr_seq);
                c_tmp = cigar_type(cigar_string, st == 1);
                c_tmp.disambiguate(rf_seq, qr_seq);
                c_tmp.trim_end(false, true, c.qr_end() - 1);
                c_tmp.check(rf_seq, qr_seq);
            }
        }
    }
}

int main(int, char * argv[])
{
    global_assert::prog_name() = argv[0];
    vector< pair< string, pair< string, string > > > test_v = {
        { "1M", { "A", "A" } },
        { "2M", { "AA", "AA" } },
        { "2M", { "AA", "AT" } },
        { "3M", { "AAA", "AAA" } },
        { "3M", { "AAA", "ATA" } },
        { "3M", { "AAA", "TTT" } },
        { "1M1I", { "A", "AT" } },
        { "1M1D", { "AT", "A" } },
        { "1M1I1M", { "AT", "ACT" } },
        { "1M1D1M", { "ACT", "AT" } }
    };
    for (size_t i = 0; i < test_v.size(); ++i)
    {
        check_cigar(test_v[i].first, test_v[i].second.first, test_v[i].second.second);
    }
}
