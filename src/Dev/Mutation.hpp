//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __MUTATION_HPP
#define __MUTATION_HPP

#include <string>
#include <iostream>

/**
 * @brief Holds information about a mutation from a base sequence.
 *
 * @tparam Seq_Type Base sequence type.
 * @tparam Size_Type Integral type used for offsets and lengths.
 */
template <class Seq_Type = std::string, typename Size_Type = size_t>
class Mutation
{
public:
    typedef typename Seq_Type::value_type value_type;

    Mutation()
    : start(0), len(0) {}

    Mutation(Size_Type _start, Size_Type _len, const Seq_Type& _seq = Seq_Type())
    : start(_start), len(_len), seq(_seq) {}

    Mutation(const Mutation& rhs)
    : start(rhs.start), len(rhs.len), seq(rhs.seq) {}

    static Mutation ins(Size_Type start, value_type symbol)
    {
        return Mutation(start, 0, Seq_Type(1, symbol));
    }

    static Mutation snp(Size_Type start, value_type symbol)
    {
        return Mutation(start, 1, Seq_Type(1, symbol));
    }

    static Mutation del(Size_Type start)
    {
        return Mutation(start, 1);
    }

    Size_Type get_start() const { return start; }
    Size_Type get_len() const { return len; }
    const Seq_Type& get_seq() const { return seq; }

    bool is_ins() const { return len == 0 and seq.size() == 1; }
    bool is_snp() const { return len == 1 and seq.size() == 1; }
    bool is_del() const { return len == 1 and seq.size() == 0; }
    bool is_empty() const { return len == 0 and seq.size() == 0; }

    friend std::ostream& operator <<(std::ostream& os, const Mutation& rhs)
    {
        os << "(start=" << (size_t)rhs.start << ",len=" << (size_t)rhs.len << ",seq=" << rhs.seq << ")";
        return os;
    }

private:
    Seq_Type seq;
    Size_Type start;
    Size_Type len;
};


#endif
