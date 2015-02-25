//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __CIGAR_HPP
#define __CIGAR_HPP

#include <string>
#include <vector>
#include <boost/range/adaptor/reversed.hpp>
#include "global.hpp"
#include "global_assert.hpp"
#include "Cigar_Op.hpp"
//#include "Util.h"
#include "RC_Sequence.hpp"
#include "logger.hpp"


namespace cigar
{

namespace ba = boost::adaptors;
using std::string;
using std::vector;

/** Cigar string holder. */
template < typename Size_Type = unsigned >
class Cigar
{
public:
    typedef Size_Type size_type;
    typedef Cigar_Op< Size_Type > cigar_op_type;
    typedef rc_sequence::Sequence< std::string > sequence_type;
    typedef sequence_type::proxy_type sequence_proxy_type;

    /** Default constructor */
    Cigar() : _rf_start(0), _qr_start(0), _rf_len(0), _qr_len(0), _reversed(false) {}

    // allow move only
    DELETE_COPY_CTOR(Cigar);
    DEFAULT_MOVE_CTOR(Cigar);
    DELETE_COPY_ASOP(Cigar);
    DEFAULT_MOVE_ASOP(Cigar);

    /** Constructor from string object.
     * @param cigar Cigar in string form.
     * @param rev Bool; if true, pos strand of rf mapped to neg strand of query.
     * @param rf_start Start of match in reference.
     * @param qr_start Start of match in query.
     * @param rf_neg_strand Cigar describes mapping of the neg strand of reference.
     */
    explicit Cigar(const string& cigar_string, bool reversed = false,
                   Size_Type rf_start = 0, Size_Type qr_start = 0, bool rf_neg_strand = false)
    : _rf_start(rf_start), _qr_start(qr_start), _rf_len(0), _qr_len(0), _reversed(reversed)
    {
        std::istringstream is(cigar_string);
        cigar_op_type op;
        // read in series of ops
        while (is >> op)
        {
            _op_vect.push_back(op);
        }
        ASSERT(is.eof());
        if (rf_neg_strand)
        {
            // if cigar describes mapping of neg rf strand, reverse the op vector
            reverse();
        }
        else
        {
            // still need to build offsets
            recompute_offsets();
        }
    }

    /** Reverse mapping coordinates on query */
    void reverse()
    {
        decltype(_op_vect) op_v_rev(_op_vect.rbegin(), _op_vect.rend());
        _op_vect = std::move(op_v_rev);
        recompute_offsets();
    }

    bool reversed() const { return _reversed; }
    bool same_st() const { return not reversed(); }
    bool diff_st() const { return reversed(); }
    size_t n_ops() const { return _op_vect.size(); }
    GETTER(Size_Type, rf_start, _rf_start)
    Size_Type rf_len() const { return _rf_len; }
    Size_Type rf_end() const { return _rf_start + _rf_len; }
    GETTER(Size_Type, qr_start, _qr_start)
    Size_Type qr_len() const { return _qr_len; }
    Size_Type qr_end() const { return _qr_start + _qr_len; }

    /** Get operation. */
    const cigar_op_type& op(size_t i) const
    {
#ifndef NDEBUG
        return _op_vect.at(i);
#else
        return _op_vect[i];
#endif
    }
    cigar_op_type& op(size_t i)
    {
        return const_cast< cigar_op_type& >(const_cast< const Cigar* >(this)->op(i));
    }
    /** Get operation type. */
    char op_type(size_t i) const { return op(i).op_type(); }
    /** Get operation length. */
    Size_Type op_len(size_t i) const { return op(i).len(); }
    /** Check if operation is a match. */
    bool op_is_match(size_t i) const { return op(i).is_match(); }
    /** Check if operation is a deletion. */
    bool op_is_deletion(size_t i) const { return op(i).is_deletion(); }
    /** Check if operation is an insertion. */
    bool op_is_insertion(size_t i) const { return op(i).is_insertion(); }
    /** Get rf start offset of operation (ignoring rf_start). */
    Size_Type op_rf_offset(size_t i) const { return i == n_ops()? _rf_len : op(i).rf_offset(); }
    /** Get qr start offset of operation (ignoring qr_start). */
    Size_Type op_qr_offset(size_t i) const { return i == n_ops()? (same_st()? _qr_len : 0) : op(i).qr_offset(); }
    /** Get qr start offset of operation (ignoring qr_start) relative to its orientation. */
    Size_Type op_qr_roffset(size_t i) const { return same_st()? op_qr_offset(i) : _qr_len - op_qr_offset(i); }
    /** Get rf start position of operation (including rf_start). */
    Size_Type op_rf_pos(size_t i) const { return rf_start() + op_rf_offset(i); }
    /** Get qr start position of operation (including qr_start). */
    Size_Type op_qr_pos(size_t i) const { return qr_start() + op_qr_offset(i); }
    /** Get rf length of sequence of operations. */
    Size_Type op_rf_len(size_t i, size_t cnt = 1) const
    {
        return op_rf_offset(i + cnt) - op_rf_offset(i);
    }
    /** Get qr length of sequence of operations. */
    Size_Type op_qr_len(size_t i, size_t cnt = 1) const
    {
        return same_st()? op_qr_offset(i + cnt) - op_qr_offset(i) : op_qr_offset(i) - op_qr_offset(i + cnt);
    }
    /** Get length of the op before given rf position.
     * @param i Index of the op to consider.
     * @param pos Absolute rf position where to stop.
     */
    Size_Type get_rf_op_prefix_len(size_t i, Size_Type pos) const
    {
        ASSERT(op_rf_pos(i) <= pos and pos <= op_rf_pos(i + 1));
        return pos - op_rf_pos(i);
    }
    /** Get length of the op before given qr position.
     * @param i Index of the op to consider.
     * @param pos Absolute qr position where to stop.
     */
    Size_Type get_qr_op_prefix_len(size_t i, Size_Type pos) const
    {
        if (same_st())
        {
            ASSERT(op_qr_pos(i) <= pos and pos <= op_qr_pos(i + 1));
            return pos - op_qr_pos(i);
        }
        else // reversed
        {
            ASSERT(op_qr_pos(i + 1) <= pos and pos <= op_qr_pos(i));
            return op_qr_pos(i) - pos;
        }
    }

    /** Build complementary cigar (qr->rf). */
    Cigar complement() const
    {
        Cigar res;
        for (size_t i = 0; i < _op_vect.size(); ++i)
        //for (const auto& op : (same_st()? _op_vect : _op_vect | ba::reversed))
        {
            size_t j = (same_st()? i : _op_vect.size() - 1 - i);
            res._op_vect.emplace_back(_op_vect[j].complement());
            //res._op_vect.emplace_back(op.complement());
        }
        res._reversed = _reversed;
        res._rf_start = _qr_start;
        res._qr_start = _rf_start;
        res.recompute_offsets();
        return res;
    }

    /** Compute sub-cigar.
     * @param i Start op.
     * @param cnt Number of ops to include.
     */
    Cigar subcigar(size_t i, size_t cnt) const
    {
        ASSERT(i + cnt <= n_ops());
        Cigar res;
        res._op_vect = decltype(_op_vect)(_op_vect.begin() + i, _op_vect.begin() + i + cnt);
        res._reversed = _reversed;
        res.rf_start() = op_rf_pos(i);
        res.qr_start() = same_st()? op_qr_pos(i) : op_qr_pos(i + cnt);
        res.recompute_offsets();
        return res;
    }

    /** Cut op.
     * @param idx Index of the op to cut.
     * @param len Op length before the cut.
     */
    void cut_op(size_t i, Size_Type len)
    {
        if (i == n_ops())
        {
            return;
        }
        ASSERT(len <= op_len(i));
        if (len == 0 or len == op_len(i))
        {
            return;
        }
        cigar_op_type& old_op = op(i);
        cigar_op_type new_op(old_op.op_type(), old_op.len() - len);
        new_op.rf_offset() = old_op.rf_offset() + (op_is_match(i) or op_is_deletion(i)? len : 0);
        new_op.qr_offset() = (op_is_match(i) or op_is_insertion(i)?
                              (same_st()? old_op.qr_offset() + len : old_op.qr_offset() - len)
                              : old_op.qr_offset());
        old_op.len() = len;
        _op_vect.insert(_op_vect.begin() + i + 1, new_op);
    }

    /** Disambiguate M operations. */
    void disambiguate(const sequence_proxy_type& rf_seq, const sequence_proxy_type& qr_seq)
    {
        ASSERT(rf_seq.size() == rf_len());
        ASSERT(qr_seq.size() == qr_len());

        LOG("cigar", debug) << ptree("disambiguate_start")
            .put("cigar", this->to_ptree())
            .put("rf_seq", rf_seq)
            //.put("qr_seq", (same_st()? qr_seq : reverseComplement(qr_seq)));
            .put("qr_seq", qr_seq.revcomp(diff_st()));

        /*
        auto get_qr_pos = [&] (Size_Type pos) -> char {
            ASSERT(not same_st() or pos < qr_len());
            ASSERT(not diff_st() or (0 < pos and pos <= qr_len()));
            //return (same_st()? qr_seq[pos] : ::complement(char(qr_seq[pos - 1])));
            return same_st()? qr_seq[pos] : dnasequence::complement::base(qr_seq[pos - 1]);
        };
        */
        for (size_t i = 0; i < n_ops(); ++i)
        {
            if (op_type(i) == 'M')
            {
                vector< cigar_op_type > v;
                //cigar_op_type op(rf_seq[op_rf_offset(i)] == get_qr_pos(op_qr_offset(i))? '=' : 'X', 1);
                cigar_op_type op(rf_seq[op_rf_offset(i)] == qr_seq.revcomp(diff_st())[op_qr_roffset(i)]? '=' : 'X', 1);
                op.rf_offset() = op_rf_offset(i);
                op.qr_offset() = op_qr_offset(i);
                size_t qr_roffset = op_qr_roffset(i);

                while (op.rf_offset() + op.len() < op_rf_offset(i + 1))
                {
                    if ((op.op_type() == '=')
                        == (rf_seq[op.rf_offset() + op.len()]
                            //== get_qr_pos((same_st()? op.qr_offset() + op.len() : op.qr_offset() - op.len()))))
                            == qr_seq.revcomp(diff_st())[qr_roffset + op.len()]))
                    {
                        ++op.len();
                    }
                    else
                    {
                        v.push_back(op);
                        op.op_type() = (op.op_type() == '='? 'X' : '=');
                        op.rf_offset() += op.len();
                        op.qr_offset() = (same_st()? op.qr_offset() + op.len() : op.qr_offset() - op.len());
                        qr_roffset += op.len();
                        op.len() = 1;
                    }
                }
                ASSERT(op.rf_offset() + op.len() == op_rf_offset(i + 1));
                v.push_back(op);
                _op_vect.erase(_op_vect.begin() + i);
                _op_vect.insert(_op_vect.begin() + i, v.begin(), v.end());
                i += v.size() - 1;
            }
        }
        LOG("cigar", debug) << ptree("disambiguate_end")
            .put("cigar", this->to_ptree());
    }

    /** Check Cigar. */
    void check(const sequence_proxy_type& rf_seq, const sequence_proxy_type& qr_seq) const
    {
        static_cast< void >(rf_seq);
        static_cast< void >(qr_seq);
#ifndef BOOST_DISABLE_ASSERTS
        ASSERT(rf_seq.size() == rf_len());
        ASSERT(qr_seq.size() == qr_len());
        Size_Type check_rf_len = 0;
        Size_Type check_qr_len = 0;
        for (size_t i = 0; i < _op_vect.size(); ++i)
        {
            // check rf_offset
            ASSERT(not (i == 0) or op_rf_pos(i) == rf_start());
            ASSERT(not (i > 0) or op_rf_pos(i) == op_rf_pos(i - 1) + op_rf_len(i - 1));
            // check qr_offset
            ASSERT(not (i == 0) or op_qr_pos(i) == (same_st()? qr_start() : qr_end()));
            ASSERT(not (i > 0) or op_qr_pos(i) == (same_st()? op_qr_pos(i - 1) + op_qr_len(i - 1) : op_qr_pos(i - 1) - op_qr_len(i - 1)));
            // recompute lengths
            check_rf_len += op_rf_len(i);
            check_qr_len += op_qr_len(i);
            // check '=' ops
            if (op_type(i) == '=')
            {
                /*
                ASSERT(rf_seq.substr(op_rf_offset(i), op_rf_len(i))
                       == (same_st()? qr_seq.substr(op_qr_offset(i), op_qr_len(i))
                           : reverseComplement(qr_seq.substr(op_qr_offset(i) - op_qr_len(i), op_qr_len(i)))));
                */
                ASSERT(rf_seq.substr(op_rf_offset(i), op_rf_len(i))
                       == qr_seq.revcomp(diff_st()).substr(op_qr_roffset(i), op_qr_len(i)));
            }
        }
        ASSERT(check_rf_len == rf_len());
        ASSERT(check_qr_len == qr_len());
#endif
    }

    /** Comparison operators. */
    friend bool operator == (const Cigar& lhs, const Cigar& rhs)
    {
        return (lhs._reversed == rhs._reversed
                and lhs._rf_start == rhs._rf_start
                and lhs._rf_len == rhs._rf_len
                and lhs._qr_start == rhs._qr_start
                and lhs._qr_len == rhs._qr_len
                and lhs._op_vect == rhs._op_vect);
    }
    friend bool operator != (const Cigar& lhs, const Cigar& rhs) { return !(lhs == rhs); }

    /** Formatted i/o operators */
    friend std::ostream& operator << (std::ostream& os, const Cigar& rhs)
    {
        for (const auto& op : rhs._op_vect)
        {
            os << op;
        }
        return os;
    }

    /** Get string representation. */
    string to_string() const
    {
        std::ostringstream os;
        os << *this;
        return os.str();
    }

    /** Get ptree representation. */
    boost::property_tree::ptree to_ptree() const
    {
        return ptree()
            .put("rf_start", rf_start())
            .put("rf_len", rf_len())
            .put("qr_start", qr_start())
            .put("qr_len", qr_len())
            .put("rc", reversed())
            .put("ops", cont_to_ptree(_op_vect));
    }

private:
    /** Recompute _rf_len, _qr_len, and all offsets stored in the op vector */
    void recompute_offsets()
    {
        _rf_len = 0;
        _qr_len = 0;
        for (auto& op : _op_vect)
        {
            op.rf_offset() = _rf_len;
            if (op.is_match() or op.is_deletion())
            {
                _rf_len += op.len();
            }
            if (op.is_match() or op.is_insertion())
            {
                _qr_len += op.len();
            }
        }
        Size_Type qr_pos = (same_st()? 0 : _qr_len);
        for (auto& op : _op_vect)
        {
            op.qr_offset() = qr_pos;
            qr_pos = (op.is_match() or op.is_insertion()?
                      (same_st()? qr_pos + op.len() : qr_pos - op.len())
                      : qr_pos);
        }
    }

    vector< cigar_op_type > _op_vect;
    Size_Type _rf_start;
    Size_Type _qr_start;
    Size_Type _rf_len;
    Size_Type _qr_len;
    bool _reversed;
}; // class Cigar

} // namespace cigar


#endif
