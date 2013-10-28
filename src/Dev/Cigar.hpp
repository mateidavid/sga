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
    /** Cigar string holder. */
    class Cigar
    {
    public:
        /** Constructor from string object.
         * @param cigar Cigar in string form.
         * @param rev Flag indicating coordinates are reversed on the query string.
         */
        Cigar(const std::string& cigar, bool reversed = false);

        /** Direction of the mapping. */
        bool is_reversed() const { return _reversed; }

        /** Get the number of operations. */
        size_t get_n_ops() const { return _op_vect.size(); }

        /** Get operation, by index. */
        char get_op(size_t i) const { return _op_vect[i]; }

        /** Get operation length, by index. */
        Size_Type get_op_len(size_t i) const { return _op_len_vect[i]; }

        /** Get reference length of operation, by index. */
        Size_Type get_rf_op_len(size_t i) const { return is_match_op(_op_vect[i]) or is_deletion_op(_op_vect[i])? _op_len_vect[i] : 0; }

        /** Get query length of operation, by index. */
        Size_Type get_qr_op_len(size_t i) const { return is_match_op(_op_vect[i]) or is_insertion_op(_op_vect[i])? _op_len_vect[i] : 0; }

        /** Get reference offset of operation, by index. */
        Size_Type get_rf_offset(size_t i) const { return _rf_offset[i]; }

        /** Get query offset of operation, by index. */
        Size_Type get_qr_offset(size_t i) const { return _qr_offset[i]; }

        /** Get total reference length. */
        Size_Type get_rf_len() const
        {
            return _op_vect.size() == 0?
                   0
                   : _rf_offset.back() + get_rf_op_len(_op_vect.size() - 1);
        }

        /** Get total reference length. */
        Size_Type get_qr_len() const
        {
            return _op_vect.size() == 0?
                   0
                   : _reversed ?
                     _qr_offset.front() + get_qr_op_len(0)
                     : _qr_offset.back() + get_qr_op_len(_op_vect.size() - 1);
        }

        /** Get complementary cigar (qr->rf). */
        Cigar complement() const;

        /** Get string representation. */
        std::string to_string() const;

        /** Is it a match operation. */
        static bool is_match_op(char op) { return match_ops.find_first_of(op) != std::string::npos; }

        /** Is it a deletion operation. */
        static bool is_deletion_op(char op) { return deletion_ops.find_first_of(op) != std::string::npos; }

        /** Is it an insertion operation. */
        static bool is_insertion_op(char op) { return insertion_ops.find_first_of(op) != std::string::npos; }

        /** Equality comparison. */
        bool operator == (const Cigar& rhs);

    private:
        std::vector<char> _op_vect;
        std::vector<Size_Type> _op_len_vect;
        std::vector<Size_Type> _rf_offset;
        std::vector<Size_Type> _qr_offset;
        bool _reversed;

        Cigar() : _reversed(false) {}

        static char complement_op(char op);

        static std::string const match_ops;
        static std::string const deletion_ops;
        static std::string const insertion_ops;
    };
}


#endif
