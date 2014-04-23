#include "Cigar.hpp"

#include <sstream>
#include <cstdlib>
#include <boost/property_tree/json_parser.hpp>

#include "indent.hpp"
#include "print_seq.hpp"
#include "../Util/Util.h"

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
        ASSERT(Cigar::is_match_op(tmp.op) or Cigar::is_deletion_op(tmp.op) or Cigar::is_insertion_op(tmp.op));
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
                                 (not _reversed? r2_pos + _op_vect[i].len : r2_pos - _op_vect[i].len) : r2_pos);
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
    {
        os << (size_t)_op_vect[i].len << _op_vect[i].op;
    }
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
        {
            return false;
        }
    }
    return true;
}

char Cigar::complement_op(char op)
{
    static const string from_op = "M=XDNPISH";
    static const string to_op = "M=XIIIDDD";
    size_t i = from_op.find_first_of(op);
    ASSERT(i != string::npos);
    return to_op[i];
}

Cigar Cigar::substring(size_t start, size_t end)
{
    ASSERT(start <= end and end <= get_n_ops());
    Cigar res;
    res._reversed = _reversed;
    Size_Type rf_start_offset = (start < get_n_ops()? _op_vect[start].rf_offset : 0);
    Size_Type rf_end_offset = (end < get_n_ops()? _op_vect[end].rf_offset : _rf_len);
    Size_Type qr_start_offset = (start < get_n_ops()? _op_vect[start].qr_offset : (not _reversed? 0 : _qr_len));
    Size_Type qr_end_offset = (end < get_n_ops()? _op_vect[end].qr_offset : (not _reversed? _qr_len : 0));
    res._rf_start = _rf_start + rf_start_offset;
    res._rf_len = rf_end_offset - rf_start_offset;
    res._qr_start = _qr_start + (not _reversed? qr_start_offset : qr_end_offset);
    res._qr_len = (not _reversed? qr_end_offset - qr_start_offset : qr_start_offset - qr_end_offset);

    Size_Type rf_offset = 0;
    for (size_t i = start; i < end; i++)
    {
        Cigar_Op tmp = _op_vect[i];
        tmp.rf_offset = rf_offset;
        if (is_match(i) or is_deletion(i))
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
    ASSERT(idx <= get_n_ops());
    if (idx == get_n_ops())
    {
        return;
    }

    ASSERT(len <= _op_vect[idx].len);
    if (len == 0 or len == _op_vect[idx].len)
    {
        return;
    }

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

void Cigar::disambiguate(const string& rf_seq, const string& qr_seq)
{
    ASSERT(rf_seq.size() == _rf_len);
    ASSERT(qr_seq.size() == _qr_len);

    //cerr << indent::tab << "disambiguating cigar:\n" << indent::inc << *this << indent::dec
    //     << indent::tab << "rf_seq= " << rf_seq
    //     << indent::nl << "qr_seq= " << (not is_reversed()? qr_seq : reverseComplement(qr_seq)) << '\n';

    auto get_qr_pos = [&] (Size_Type pos) -> char {
        ASSERT(_reversed or pos < _qr_len);
        ASSERT(not _reversed or (0 < pos and pos <= _qr_len));
        return (not _reversed? qr_seq[pos] : ::complement(char(qr_seq[pos - 1])));
    };
    for (size_t i = 0; i < get_n_ops(); ++i)
    {
        if (get_op(i) == 'M')
        {
            vector< Cigar_Op > v;
            Cigar_Op op;
            op.rf_offset = _op_vect[i].rf_offset;
            op.qr_offset = _op_vect[i].qr_offset;
            op.op = (rf_seq[op.rf_offset] == get_qr_pos(op.qr_offset)? '=' : 'X');
            op.len = 1;

            while (op.rf_offset + op.len < _op_vect[i].rf_offset + _op_vect[i].len)
            {
                if ((op.op == '=')
                        == (rf_seq[op.rf_offset + op.len]
                            == get_qr_pos((not _reversed? op.qr_offset + op.len : op.qr_offset - op.len))))
                {
                    ++op.len;
                }
                else
                {
                    v.push_back(op);
                    op.op = (op.op == '='? 'X' : '=');
                    op.rf_offset += op.len;
                    op.qr_offset = (not _reversed? op.qr_offset + op.len : op.qr_offset - op.len);
                    op.len = 1;
                }
            }
            ASSERT(op.rf_offset + op.len == _op_vect[i].rf_offset + _op_vect[i].len);
            v.push_back(op);

            _op_vect.erase(_op_vect.begin() + i);
            _op_vect.insert(_op_vect.begin() + i, v.begin(), v.end());
            i += v.size() - 1;
        }
    }
    //cerr << indent::tab << "result:\n" << indent::inc << *this << indent::dec;
}

bool Cigar::check(const string& rf_seq, const string& qr_seq) const
{
    ASSERT(rf_seq.size() == get_rf_len());
    ASSERT(qr_seq.size() == get_qr_len());
    Size_Type check_rf_len = 0;
    Size_Type check_qr_len = 0;
    for (size_t i = 0; i < _op_vect.size(); ++i)
    {
        // check rf_offset
        ASSERT(not (i == 0) or get_rf_offset(i) == get_rf_start());
        ASSERT(i == 0 or get_rf_offset(i) == get_rf_offset(i-1) + get_rf_op_len(i-1));
        // check qr_offset
        ASSERT(not (i == 0) or get_qr_offset(i) == (not _reversed? get_qr_start() : get_qr_start() + get_qr_len()));
        ASSERT(i == 0 or get_qr_offset(i) == (not _reversed? get_qr_offset(i-1) + get_qr_op_len(i-1) : get_qr_offset(i-1) - get_qr_op_len(i-1)));
        // recompute lengths
        check_rf_len += get_rf_op_len(i);
        check_qr_len += get_qr_op_len(i);
        // check '=' ops
        if (get_op(i) == '=')
        {
            ASSERT(rf_seq.substr(get_rf_offset(i) - get_rf_start(), get_rf_op_len(i))
                   == (not _reversed? qr_seq.substr(get_qr_offset(i) - get_qr_start(), get_qr_op_len(i))
                       : reverseComplement(qr_seq.substr(get_qr_offset(i) - get_qr_op_len(i) - get_qr_start(), get_qr_op_len(i)))));
        }
    }
    ASSERT(check_rf_len == get_rf_len());
    ASSERT(check_qr_len == get_qr_len());
    return true;
}

ostream& operator << (ostream& os, const Cigar_Op& rhs)
{
    os << '(' << (size_t)rhs.len << rhs.op << ",rf_offset=" << (size_t)rhs.rf_offset << ",qr_offset=" << (size_t)rhs.qr_offset << ')';
    return os;
}

ostream& operator << (ostream& os, const Cigar& rhs)
{
    /*
    os << indent::tab << "(Cigar" << indent::inc
       << indent::nl << "rf_start=" << rhs._rf_start << ",rf_len=" << rhs._rf_len
       << indent::nl << "qr_start=" << rhs._qr_start << ",qr_len=" << rhs._qr_len
       << indent::nl << "reversed=" << rhs._reversed
       << indent::nl << "ops:" << indent::inc << '\n';
    print_seq(os, rhs._op_vect, indent::nl, indent::tab, '\n');
    os << indent::dec << indent::dec << indent::tab << ")\n";
    */
    boost::property_tree::write_json(os, rhs.to_ptree(), false);
    return os;
}

boost::property_tree::ptree Cigar_Op::to_ptree() const
{
    boost::property_tree::ptree pt;
    ostringstream tmp;
    tmp << len << op;
    pt.put("op", tmp.str());
    pt.put("rf_offset", rf_offset);
    pt.put("qr_offset", qr_offset);
    return pt;
}

boost::property_tree::ptree Cigar::to_ptree() const
{
    boost::property_tree::ptree pt;
    pt.put("rf_start", get_rf_start());
    pt.put("rf_len", get_rf_len());
    pt.put("qr_start", get_qr_start());
    pt.put("qr_len", get_qr_len());
    pt.put("rc", is_reversed());
    pt.put_child("ops", cont_to_ptree(_op_vect));
    return pt;
}

} // namespace MAC
