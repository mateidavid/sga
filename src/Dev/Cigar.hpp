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
    };

    /** Cigar string holder. */
    class Cigar
    {
    public:
        /** Constructor from string object.
         * @param cigar Cigar in string form.
         * @param rev Flag indicating coordinates are reversed on the query string.
         */
        Cigar(const std::string& cigar, bool reversed = false, Size_Type rf_start = 0, Size_Type qr_start = 0);

        bool is_reversed() const { return _reversed; }
        size_t get_n_ops() const { return _op_vect.size(); }
        Size_Type get_rf_len() const { return _rf_len; }
        Size_Type get_qr_len() const { return _qr_len; }
        Size_Type get_rf_start() const { return _rf_start; }
        Size_Type get_qr_start() const { return _qr_start; }

        /** Get operation struct, by index. */
        const Cigar_Op& get_op_struct(size_t i) const { return _op_vect[i]; }

        /** Get operation. */
        char get_op(size_t i) const { return _op_vect[i].op; }

        /** Get operation length, by index. */
        Size_Type get_op_len(size_t i) const { return _op_vect[i].len; }

        /** Get reference length of operation, by index. */
        Size_Type get_rf_op_len(size_t i) const { return is_match(i) or is_deletion(i)? _op_vect[i].len : 0; }

        /** Get query length of operation, by index. */
        Size_Type get_qr_op_len(size_t i) const { return is_match(i) or is_insertion(i)? _op_vect[i].len : 0; }

        /** Get reference offset of operation, by index. */
        Size_Type get_rf_offset(size_t i) const { return _rf_start + _op_vect[i].rf_offset; }

        /** Get query offset of operation, by index. */
        Size_Type get_qr_offset(size_t i) const { return _qr_start + _op_vect[i].qr_offset; }

        void set_rf_start(Size_Type rf_start) { _rf_start = rf_start; }
        void set_qr_start(Size_Type qr_start) { _qr_start = qr_start; }

        bool is_match(size_t i) const { return Cigar::is_match_op(_op_vect[i].op); }
        bool is_deletion(size_t i) const { return Cigar::is_deletion_op(_op_vect[i].op); }
        bool is_insertion(size_t i) const { return Cigar::is_insertion_op(_op_vect[i].op); }

        /** Get complementary cigar (qr->rf). */
        Cigar complement() const;

        /** Compute sub-cigar.
         * @param start Start op.
         * @param len Number of ops.
         */
        Cigar substring(size_t start, size_t len);

        /** Cut op.
         * @param idx Index of the op to cut.
         * @param len Op length before the cut.
         */
        void cut_op(size_t idx, Size_Type len);

        /** Get string representation. */
        std::string to_string() const;

        /** Equality comparison. */
        bool operator == (const Cigar& rhs);

        friend std::ostream& operator << (std::ostream& os, const Cigar& rhs);

    private:
        std::vector< Cigar_Op > _op_vect;
        Size_Type _rf_start;
        Size_Type _qr_start;
        Size_Type _rf_len;
        Size_Type _qr_len;
        bool _reversed;

        Cigar() : _reversed(false) {}
        void set_qr_offsets();

        static char complement_op(char op);

        /** Is it a match operation. */
        static bool is_match_op(char op) { return match_ops.find_first_of(op) != std::string::npos; }

        /** Is it a deletion operation. */
        static bool is_deletion_op(char op) { return deletion_ops.find_first_of(op) != std::string::npos; }

        /** Is it an insertion operation. */
        static bool is_insertion_op(char op) { return insertion_ops.find_first_of(op) != std::string::npos; }

        static std::string const match_ops;
        static std::string const deletion_ops;
        static std::string const insertion_ops;
    };
}


#endif
