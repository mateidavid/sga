//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __CIGAR_HPP
#define __CIGAR_HPP

#include <string>
#include <vector>

#include "MAC_forward.hpp"


namespace MAC
{

struct Cigar_Op
{
    Size_Type len;
    Size_Type rf_offset;
    Size_Type qr_offset;
    char op;
    boost::property_tree::ptree to_ptree() const;
};

/** Cigar string holder. */
class Cigar
{
public:
    // allow move only
    DELETE_COPY_CTOR(Cigar)
    DEFAULT_MOVE_CTOR(Cigar)
    DELETE_COPY_ASOP(Cigar)
    DEFAULT_MOVE_ASOP(Cigar)

    /** Constructor from string object.
     * @param cigar Cigar in string form.
     * @param rev Flag indicating coordinates are reversed on the query string.
     * @param rf_start Start of match in reference.
     * @param qr_start Start of match in query.
     */
    Cigar(const std::string& cigar = "", bool reversed = false, Size_Type rf_start = 0, Size_Type qr_start = 0);

    bool is_reversed() const { return _reversed; }
    size_t get_n_ops() const { return _op_vect.size(); }
    Size_Type get_rf_start() const { return _rf_start; }
    Size_Type get_rf_len() const { return _rf_len; }
    Size_Type get_qr_start() const { return _qr_start; }
    Size_Type get_qr_len() const { return _qr_len; }

    /** Get operation. */
    char get_op(size_t i) const
    {
        ASSERT(i < get_n_ops());
        return _op_vect[i].op;
    }

    /** Get operation length, by index. */
    Size_Type get_op_len(size_t i) const
    {
        ASSERT(i < get_n_ops());
        return _op_vect[i].len;
    }

    /** Get reference length of operation, by index. */
    Size_Type get_rf_op_len(size_t i) const
    {
        ASSERT(i <= get_n_ops());
        return (i < get_n_ops() and (is_match(i) or is_deletion(i)) ? _op_vect[i].len : 0);
    }

    /** Get query length of operation, by index. */
    Size_Type get_qr_op_len(size_t i) const
    {
        ASSERT(i <= get_n_ops());
        return (i < get_n_ops() and (is_match(i) or is_insertion(i)) ? _op_vect[i].len : 0);
    }

    /** Get rf length spanned by cigar ops [start...end-1]. */
    Size_Type get_rf_sub_len(size_t start, size_t end) const
    {
        ASSERT(start <= end and end <= get_n_ops());
        return get_rf_offset(end) - get_rf_offset(start);
    }

    /** Get qr length spanned by cigar ops [start...end-1]. */
    Size_Type get_qr_sub_len(size_t start, size_t end) const
    {
        ASSERT(start <= end and end <= get_n_ops());
        return (not _reversed? get_qr_offset(end) - get_qr_offset(start) : get_qr_offset(start) - get_qr_offset(end));
    }

    /** Get reference offset of operation, by index. */
    Size_Type get_rf_offset(size_t i) const
    {
        ASSERT(i <= get_n_ops());
        return _rf_start + (i < get_n_ops() ? _op_vect[i].rf_offset : _rf_len);
    }

    /** Get query offset of operation, by index. */
    Size_Type get_qr_offset(size_t i) const
    {
        ASSERT(i <= get_n_ops());
        return _qr_start + (i < get_n_ops() ? _op_vect[i].qr_offset : (not _reversed? _qr_len : 0));
    }

    /** Get length of the op before given rf position.
     * @param i Index of the op to consider.
     * @param pos Absolute rf position where to stop.
     */
    Size_Type get_rf_op_prefix_len(size_t i, Size_Type pos) const
    {
        ASSERT(i < get_n_ops());
        ASSERT(get_rf_offset(i) <= pos and pos <= get_rf_offset(i + 1));
        return pos - get_rf_offset(i);
    }

    /** Get length of the op before given qr position.
     * @param i Index of the op to consider.
     * @param pos Absolute qr position where to stop.
     */
    Size_Type get_qr_op_prefix_len(size_t i, Size_Type pos) const
    {
        ASSERT(i < get_n_ops());
        if (not _reversed)
        {
            ASSERT(get_qr_offset(i) <= pos and pos <= get_qr_offset(i + 1));
            return pos - get_qr_offset(i);
        }
        else // _reversed
        {
            ASSERT(get_qr_offset(i + 1) <= pos and pos <= get_qr_offset(i));
            return get_qr_offset(i) - pos;
        }
    }

    void set_rf_start(Size_Type rf_start)
    {
        _rf_start = rf_start;
    }
    void set_qr_start(Size_Type qr_start)
    {
        _qr_start = qr_start;
    }

    bool is_match(size_t i) const
    {
        ASSERT(i < get_n_ops());
        return Cigar::is_match_op(_op_vect[i].op);
    }
    bool is_deletion(size_t i) const
    {
        ASSERT(i < get_n_ops());
        return Cigar::is_deletion_op(_op_vect[i].op);
    }
    bool is_insertion(size_t i) const
    {
        ASSERT(i < get_n_ops());
        return Cigar::is_insertion_op(_op_vect[i].op);
    }

    /** Get complementary cigar (qr->rf). */
    Cigar complement() const;

    /** Compute sub-cigar.
     * @param start Start op.
     * @param end First op to leave out.
     */
    Cigar substring(size_t start, size_t end);

    /** Cut op.
     * @param idx Index of the op to cut.
     * @param len Op length before the cut.
     */
    void cut_op(size_t idx, Size_Type len);

    /** Disambiguate M operations. */
    void disambiguate(const std::string& rf_seq, const std::string& qr_seq);

    /** Check Cigar. */
    bool check(const std::string& rf_seq, const std::string& qr_seq) const;

    /** Get string representation. */
    std::string to_string() const;

    /** Equality comparison. */
    bool operator == (const Cigar& rhs);

    friend std::ostream& operator << (std::ostream& os, const Cigar& rhs);
    boost::property_tree::ptree to_ptree() const;

private:
    std::vector< Cigar_Op > _op_vect;
    Size_Type _rf_start;
    Size_Type _qr_start;
    Size_Type _rf_len;
    Size_Type _qr_len;
    bool _reversed;

    void set_qr_offsets();

    static char complement_op(char op);

    /** Is it a match operation. */
    static bool is_match_op(char op)
    {
        return match_ops.find_first_of(op) != std::string::npos;
    }

    /** Is it a deletion operation. */
    static bool is_deletion_op(char op)
    {
        return deletion_ops.find_first_of(op) != std::string::npos;
    }

    /** Is it an insertion operation. */
    static bool is_insertion_op(char op)
    {
        return insertion_ops.find_first_of(op) != std::string::npos;
    }

    static std::string const match_ops;
    static std::string const deletion_ops;
    static std::string const insertion_ops;
}; // class Cigar

} // namespace MAC


#endif
