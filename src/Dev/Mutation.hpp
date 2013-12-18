//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __MUTATION_HPP
#define __MUTATION_HPP

#include <iostream>
#include <cassert>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>
#include "MAC_forward.hpp"
#include "Cigar.hpp"
#include "../Util/Util.h"


namespace MAC
{
    /** Holds information about a mutation from a base sequence.
     *
     * The class holds the start and span of the base sequence region affected by a mutation,
     * as well as the alternate sequence.
     */
    class Mutation
    {
    public:
        /** Type of key used to store Mutation objects. */
        typedef Size_Type key_type;

        typedef std::function<const Mutation*(const Mutation&)> add_mut_mod_type;

        /** Default constructor. */
        Mutation()
        : _start(0), _len(0), _seq_len(0) {}

        /** Constructor.
         * @param start Start of the mutation, i.e., length of base sequence prior to mutation.
         * @param len Length of the base sequence affected by the mutation.
         * @param seq Alternate sequence.
         */
        Mutation(Size_Type start, Size_Type len, const Seq_Type& seq = Seq_Type())
        : _seq(seq), _start(start), _len(len), _seq_len(seq.size()) {}

        /** Constructor.
         * @param start Start of the mutation, i.e., length of base sequence prior to mutation.
         * @param len Length of the base sequence affected by the mutation.
         * @param seq_len Length of alternate sequence.
         */
        Mutation(Size_Type start, Size_Type len, Size_Type seq_len)
        : _start(start), _len(len), _seq_len(seq_len) {}

        /** Copy constructor. */
        Mutation(const Mutation& rhs)
        : _seq(rhs._seq), _start(rhs._start), _len(rhs._len), _seq_len(rhs._seq_len) {}

        /** Static insertion constructor.
         * @param start Start of the mutation, i.e., length of base sequence prior to mutation.
         * @param seq Inserted sequence.
         */
        static Mutation make_ins(Size_Type start, const Seq_Type& seq)
        {
            assert(seq.size() > 0);
            return Mutation(start, 0, seq);
        }

        /** Static insertion constructor.
         * @param start Start of the mutation, i.e., length of base sequence prior to mutation.
         * @param seq_len Length of the inserted sequence.
         */
        static Mutation make_ins(Size_Type start, Size_Type seq_len)
        {
            assert(seq_len > 0);
            return Mutation(start, 0, seq_len);
        }

        /** Static mnp (multi-nucleotide polymorphism) constructor.
         * @param start Start of the mutation, i.e., length of base sequence prior to mutation.
         * @param seq Alternate sequence.
         */
        static Mutation make_mnp(Size_Type start, const Seq_Type& seq)
        {
            assert(seq.size() > 0);
            return Mutation(start, seq.size(), seq);
        }

        /** Static mnp constructor.
         * @param start Start of the mutation, i.e., length of base sequence prior to mutation.
         * @param seq_len Length of the alternate sequence.
         */
        static Mutation make_mnp(Size_Type start, Size_Type seq_len)
        {
            assert(seq_len > 0);
            return Mutation(start, seq_len, seq_len);
        }

        /** Static deletion constructor.
         * @param start Start of the mutation, i.e., length of base sequence prior to mutation.
         * @param len Length of the deletion.
         */
        static Mutation make_del(Size_Type start, Size_Type len)
        {
            return Mutation(start, len);
        }

        /** @name Getters */
        /**@{*/
        Size_Type get_start() const { return _start; }
        Size_Type get_len() const { return _len; }
        Size_Type get_end() const { return _start + _len; }
        Size_Type get_seq_len() const { return _seq_len; }
        bool have_seq() const { return _seq.size() == _seq_len; }
        const Seq_Type& get_seq() const { return _seq; }
        key_type get_key() const { return _start; }
        /**@}*/

        /** @name Basic queries */
        /**@{*/
        bool is_ins() const { return _len == 0 and _seq_len > 0; }
        bool is_snp() const { return _len == 1 and _len == _seq_len; }
        bool is_del() const { return _len > 0 and _seq_len == 0; }
        bool is_empty() const { return _len == 0 and _seq_len == 0; }
        /**@}*/

        /** Merge with given Mutation.
         * Pre: Mutations must be adjacent on rf.
         * @param rhs Next Mutation.
         */
        void merge(const Mutation& rhs)
        {
            assert(get_end() == rhs.get_start());
            assert(have_seq() == rhs.have_seq());
            _len += rhs._len;
            _seq_len += rhs._seq_len;
            _seq += rhs._seq;
        }

         /** Cut mutation at given offsets.
         * @param base_offset Base offset, 0-based.
         * @param alt_offset Alternate sequence offset, 0-based.
         * @return The part of the original mutation that was cut from this object.
         */
        Mutation cut(Size_Type base_offset, Size_Type alt_offset);

        /** Simplify Mutation by dropping the ends of rf and qr if they match.
         * @param rf Reference sequence spanned by the mutation.
         */
        void simplify(const Seq_Type& rf);

        /** Reverse the mutation.
         * @param c_len The contig length.
         */
        void reverse(Size_Type c_len)
        {
            _start = c_len - (_start + _len);
            if (have_seq())
                _seq = reverseComplement(_seq);
        }

        /** Change Mutation base by adding a prefix of the given length. */
        void add_base_prefix(Size_Type len) { _start += len; }

        bool operator == (const Mutation&) const;

        friend std::ostream& operator << (std::ostream&, const Mutation&);

    private:
        Seq_Type _seq;
        Size_Type _start;
        Size_Type _len;
        Size_Type _seq_len;
    };

    typedef const Mutation* Mutation_CPtr;

    namespace detail
    {
        /** Key extractor struct for boost::multi_index_container. */
        struct Mutation_Key
        {
            typedef Mutation::key_type result_type;
            result_type operator() (const Mutation& m) const { return m.get_key(); }
            result_type operator() (const Mutation_CPtr& m_cptr) const { return m_cptr->get_key(); }
        };
    }

    /** Container for Mutation objects. */
    typedef boost::multi_index_container<
      Mutation,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_non_unique<
          detail::Mutation_Key
        >
      >
    > Mutation_Cont;

    /** Create a set of mutations to a reference string based on a cigar object.
     * Pre: Cigar contains no 'M' operations (use disambiguate() first).
     * Post: Adjacent non-match operations are merged.
     * @param cigar Cigar object describing the match.
     * @param qr Query string; optional: if not given, Mutation objects contain alternate string lengths only.
     * @return Container of Mutation objects.
     */
    std::shared_ptr< Mutation_Cont > make_mutations_from_cigar(const Cigar& cigar, const std::string& qr = std::string());

    /** Add Mutation to container, use existing Mutation if it already exists.
     * @param mut_cont Mutation container.
     * @param mut Mutation to add.
     * @return Pointer to Mutation inside container.
     */
    Mutation_CPtr add_mut_to_cont(Mutation_Cont& mut_cont, const Mutation& mut);

    /** Mutation Translation object. */
    struct Mutation_Trans
    {
        Mutation_CPtr old_mut_cptr;
        Mutation_CPtr new_mut_cptr;
        std::vector< std::tuple< Mutation_CPtr, Size_Type, Size_Type > > new_mut_rev_list;
    };

    /** Container for Mutation Translation objects. */
    typedef boost::multi_index_container<
      Mutation_Trans,
      boost::multi_index::indexed_by<
        boost::multi_index::hashed_unique<
          boost::multi_index::member< Mutation_Trans, Mutation_CPtr, &Mutation_Trans::old_mut_cptr >
        >
      >
    > Mutation_Trans_Cont;
}


#endif
