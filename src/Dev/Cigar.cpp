#include "Cigar.hpp"

#include <sstream>
#include <cstdlib>
#include <cassert>

using namespace std;


namespace MAC
{
    Cigar::Cigar(const string& cigar, bool reversed)
    : _reversed(reversed)
    {
        istringstream is(cigar);
        char op;
        size_t op_len;
        Size_Type rf_len = 0;
        Size_Type qr_len = 0;
        while (is >> op_len >> op)
        {
            assert(is_match_op(op) or is_deletion_op(op) or is_insertion_op(op));
            _op_vect.push_back(op);
            _op_len_vect.push_back(op_len);
            if (is_match_op(op) or is_deletion_op(op))
                rf_len += (Size_Type)op_len;
            if (is_match_op(op) or is_insertion_op(op))
                qr_len += (Size_Type)op_len;
        }

        Size_Type r1_pos = 0;
        Size_Type r2_pos = (reversed? qr_len : 0);
        for (size_t i = 0; i < _op_vect.size(); ++i)
        {
            Size_Type r1_pos_next = (is_match_op(_op_vect[i]) or is_deletion_op(_op_vect[i])?
                                     r1_pos + _op_len_vect[i] :
                                     r1_pos);
            _rf_offset.push_back(r1_pos);
            Size_Type r2_pos_next = (is_match_op(_op_vect[i]) or is_insertion_op(_op_vect[i])?
                                     (reversed? r2_pos - _op_len_vect[i] : r2_pos + _op_len_vect[i]) :
                                     r2_pos);
            _qr_offset.push_back(min(r2_pos, r2_pos_next));
            r1_pos = r1_pos_next;
            r2_pos = r2_pos_next;
        }
    }

    Cigar Cigar::complement() const
    {
        Cigar res;
        res._reversed = _reversed;
        for (size_t i = 0; i < _op_vect.size(); ++i)
        {
            size_t j = _reversed? _op_vect.size() - 1 - i : i;
            res._op_vect.push_back(complement_op(_op_vect[j]));
            res._op_len_vect.push_back(_op_len_vect[j]);
            res._rf_offset.push_back(_qr_offset[j]);
            res._qr_offset.push_back(_rf_offset[j]);
        }
        return res;
    }

    string Cigar::to_string() const
    {
        ostringstream os;
        for (size_t i = 0; i < _op_vect.size(); ++i)
            os << (size_t)_op_len_vect[i] << _op_vect[i];
        return os.str();
    }

    bool Cigar::operator == (const Cigar& rhs)
    {
        if (_op_vect.size() != rhs._op_vect.size()) return false;
        if (_reversed != rhs._reversed) return false;
        for (size_t i = 0; i < _op_vect.size(); ++i)
        {
            if (_op_vect[i] != rhs._op_vect[i]
                or _op_len_vect[i] != rhs._op_len_vect[i]
                or _rf_offset[i] != rhs._rf_offset[i]
                or _qr_offset[i] != rhs._qr_offset[i])
                return false;
        }
        return true;
    }

    char Cigar::complement_op(char op)
    {
        static const string from_op = "M=XDNPISH";
        static const string to_op = "M=XIIIDDD";
        size_t i = from_op.find_first_of(op);
        assert(i != string::npos);
        return to_op[i];
    }

    string const Cigar::match_ops = "M=X";
    string const Cigar::deletion_ops = "DNP";
    string const Cigar::insertion_ops = "ISH";
}
