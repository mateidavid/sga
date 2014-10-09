//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __CIGAR_OP_HPP
#define __CIGAR_OP_HPP

#include <string>
#include <iostream>
#include <sstream>
#include "shortcuts.hpp"
#include "ptree.hpp"


namespace cigar
{

using std::string;

template < typename Size_Type >
class Cigar_Op
{
public:
    typedef Size_Type size_type;

    // Allow def ctor, copy&move ctor&asop
    DEFAULT_DEF_CTOR(Cigar_Op)
    DEFAULT_COPY_CTOR(Cigar_Op)
    DEFAULT_MOVE_CTOR(Cigar_Op)
    DEFAULT_COPY_ASOP(Cigar_Op)
    DEFAULT_MOVE_ASOP(Cigar_Op)

    /** Explicit constructor. */
    Cigar_Op(char op_type, Size_Type len, Size_Type rf_offset = 0, Size_Type qr_offset = 0)
        : _len(len), _rf_offset(rf_offset), _qr_offset(qr_offset), _op_type(op_type)
    {
        if (not (is_match() or is_insertion() or is_deletion()))
        {
            _op_type = '?';
            _len = 0;
        }
    }

    /** Get complement operation. */
    Cigar_Op complement() const
    {
        return Cigar_Op(complement_op_type(op_type()), len());
    }

    /** Getters */
    GETTER(char, op_type, _op_type)
    GETTER(Size_Type, len, _len)
    GETTER(Size_Type, rf_offset, _rf_offset)
    GETTER(Size_Type, qr_offset, _qr_offset)

    /** Bool conversion */
    BOOL_CONVERSION(Cigar_Op, public, (_len > 0))

    /** Check operation type. */
    bool is_match() const { return match_ops().find_first_of(op_type()) != string::npos; }
    bool is_deletion() const { return deletion_ops().find_first_of(op_type()) != string::npos; }
    bool is_insertion() const { return insertion_ops().find_first_of(op_type()) != string::npos; }

    /** Comparison operators */
    friend bool operator == (const Cigar_Op& lhs, const Cigar_Op& rhs)
    {
        return lhs.op_type() == rhs.op_type()
            and lhs.len() == rhs.len()
            and lhs.rf_offset() == rhs.rf_offset()
            and lhs.qr_offset() == rhs.qr_offset();
    }
    friend bool operator != (const Cigar_Op& lhs, const Cigar_Op& rhs) { return !(lhs == rhs); }

    /** Formatted i/o operators */
    friend std::istream& operator >> (std::istream& is, Cigar_Op& op)
    {
        size_t len;
        char op_type;
        if (is >> len >> op_type)
        {
            op = Cigar_Op(op_type, len);
            // validate input
            if (not op)
            {
                is.setstate(std::ios_base::failbit);
            }
        }
        return is;
    }
    friend std::ostream& operator << (std::ostream& os, const Cigar_Op& op)
    {
        os << static_cast< size_t >(op.len()) << op.op_type();
        return os;
    }

    /** ptree output */
    boost::property_tree::ptree to_ptree() const
    {
        std::ostringstream tmp;
        tmp << *this;
        return ptree().put("op", tmp.str())
            .put("rf_offset", rf_offset())
            .put("qr_offset", qr_offset());
    }

private:
    static char complement_op_type(char op_type)
    {
        static const string from_op = "M=XDNPISH";
        static const string to_op = "M=XIIIDDD";
        auto i = from_op.find_first_of(op_type);
        return to_op.at(i);
    }

    static const string& match_ops()
    {
        static const string _match_ops = "M=X";
        return _match_ops;
    }
    static const string& deletion_ops()
    {
        static const string _deletion_ops = "DNP";
        return _deletion_ops;
    }
    static const string& insertion_ops()
    {
        static const string _insertion_ops = "ISH";
        return _insertion_ops;
    }

    Size_Type _len;
    Size_Type _rf_offset;
    Size_Type _qr_offset;
    char _op_type;
}; // class Cigar_Op

} // namespace cigar


#endif
