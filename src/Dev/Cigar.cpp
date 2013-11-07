#include "Cigar.hpp"

#include <sstream>
#include <cstdlib>
#include <cassert>

#include "indent.hpp"
#include "print_seq.hpp"

using namespace std;


namespace MAC
{
    Cigar::Cigar(const string& cigar, bool reversed, Size_Type rf_start, Size_Type qr_start)
    : _rf_start(rf_start), _qr_start(qr_start), _rf_len(0), _qr_len(0), _reversed(reversed)
    {
        istringstream is(cigar);
        Cigar_Op tmp;
        size_t tmp_len;
        while (is >> tmp_len >> tmp.op)
        {
            assert(Cigar::is_match_op(tmp.op) or Cigar::is_deletion_op(tmp.op) or Cigar::is_insertion_op(tmp.op));
            tmp.len = (Size_Type)tmp_len;
            tmp.rf_offset = _rf_len;
            _op_vect.push_back(tmp);

            if (Cigar::is_match_op(tmp.op) or Cigar::is_deletion_op(tmp.op))
            {
                _rf_len += tmp.len;
            }
            if (Cigar::is_match_op(tmp.op) or Cigar::is_insertion_op(tmp.op))
            {
                _qr_len += tmp.len;
            }
        }
        set_qr_offsets();
    }

    void Cigar::set_qr_offsets()
    {
        Size_Type r2_pos = (_reversed? _qr_len : 0);
        for (size_t i = 0; i < _op_vect.size(); ++i)
        {
            Size_Type r2_pos_next = (is_match(i) or is_insertion(i)?
                                     (not _reversed? r2_pos + _op_vect[i].len : r2_pos - _op_vect[i].len) :
                                     r2_pos);
            _op_vect[i].qr_offset = r2_pos;
            r2_pos = r2_pos_next;
        }
    }

    Cigar Cigar::complement() const
    {
        Cigar res;
        res._reversed = _reversed;
        res._rf_start = _qr_start;
        res._rf_len = _qr_len;
        res._qr_start = _rf_start;
        res._qr_len = _rf_len;
        for (size_t i = 0; i < _op_vect.size(); ++i)
        {
            size_t j = _reversed? _op_vect.size() - 1 - i : i;
            Cigar_Op tmp;
            tmp.op = complement_op(_op_vect[j].op);
            tmp.len = _op_vect[j].len;
            tmp.rf_offset = (not _reversed? _op_vect[j].qr_offset : _op_vect[j].qr_offset - get_qr_op_len(j));
            tmp.qr_offset = (not _reversed? _op_vect[j].rf_offset : _op_vect[j].rf_offset + get_rf_op_len(j));
            res._op_vect.push_back(tmp);
        }
        return res;
    }

    string Cigar::to_string() const
    {
        ostringstream os;
        for (size_t i = 0; i < _op_vect.size(); ++i)
            os << (size_t)_op_vect[i].len << _op_vect[i].op;
        return os.str();
    }

    bool Cigar::operator == (const Cigar& rhs)
    {
        if (_op_vect.size() != rhs._op_vect.size()) return false;
        if (_reversed != rhs._reversed) return false;
        if (_rf_start != rhs._rf_start) return false;
        if (_rf_len != rhs._rf_len) return false;
        if (_qr_start != rhs._qr_start) return false;
        if (_qr_len != rhs._qr_len) return false;
        for (size_t i = 0; i < _op_vect.size(); ++i)
        {
            if (_op_vect[i].op != rhs._op_vect[i].op
                or _op_vect[i].len != rhs._op_vect[i].len
                or _op_vect[i].rf_offset != rhs._op_vect[i].rf_offset
                or _op_vect[i].qr_offset != rhs._op_vect[i].qr_offset)
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

    Cigar Cigar::substring(size_t start, size_t len)
    {
        assert(start + len <= get_n_ops());
        Cigar res;
        res._reversed = _reversed;
        Size_Type rf_start_offset = (start < get_n_ops()? _op_vect[start].rf_offset : 0);
        Size_Type rf_end_offset = (start + len < get_n_ops()? _op_vect[start + len].rf_offset : _rf_len);
        Size_Type qr_start_offset = (start < get_n_ops()? _op_vect[start].qr_offset : (not _reversed? 0 : _qr_len));
        Size_Type qr_end_offset = (start + len < get_n_ops()? _op_vect[start + len].qr_offset : (not _reversed? _qr_len : 0));
        res._rf_start = _rf_start + rf_start_offset;
        res._rf_len = rf_end_offset - rf_start_offset;
        res._qr_start = _qr_start + (not _reversed? qr_start_offset : qr_end_offset);
        res._qr_len = (not _reversed? qr_end_offset - qr_start_offset : qr_start_offset - qr_end_offset);

        Size_Type rf_offset = 0;
        for (size_t i = 0; i < len and start + i < get_n_ops(); i++)
        {
            Cigar_Op tmp = _op_vect[start + i];
            tmp.rf_offset = rf_offset;
            if (is_match(start + i) or is_deletion(start + i))
            {
                rf_offset += tmp.len;
            }
            res._op_vect.push_back(tmp);
        }

        res.set_qr_offsets();
        return res;
    }

    void Cigar::cut_op(size_t idx, Size_Type len)
    {
        assert(idx <= get_n_ops());
        if (idx == get_n_ops())
            return;
        assert(len <= _op_vect[idx].len);
        if (len == 0 or len == _op_vect[idx].len)
            return;
        Cigar_Op& old_op = _op_vect[idx];
        Cigar_Op new_op;
        new_op.op = old_op.op;
        new_op.len = old_op.len - len;
        new_op.rf_offset = old_op.rf_offset + (is_match(idx) or is_deletion(idx)? len : 0);
        new_op.qr_offset = (is_deletion(idx)? old_op.qr_offset : (not _reversed? old_op.qr_offset + len : old_op.qr_offset - len));
        old_op.len = len;
        _op_vect.insert(_op_vect.begin() + idx + 1, new_op);
    }

    string const Cigar::match_ops = "M=X";
    string const Cigar::deletion_ops = "DNP";
    string const Cigar::insertion_ops = "ISH";

    ostream& operator << (ostream& os, const Cigar_Op& rhs)
    {
        os << (size_t)rhs.len << rhs.op << ",rf_offset=" << (size_t)rhs.rf_offset << ",qr_offset=" << (size_t)rhs.qr_offset;
        return os;
    }

    ostream& operator << (ostream& os, const Cigar& rhs)
    {
        os << indent::nl << "(Cigar" << indent::inc
           << indent::nl << "rf_start=" << rhs._rf_start << ",rf_len=" << rhs._rf_len
           << indent::nl << "qr_start=" << rhs._qr_start << ",qr_len=" << rhs._qr_len
           << indent::nl << "reversed=" << rhs._reversed
           << indent::nl << "ops=" << indent::inc;
        print_seq(os, rhs._op_vect, indent::nl, indent::nl);
        os << indent::dec << indent::dec << indent::nl << ")";
    }
}
