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
#include <boost/tuple/tuple.hpp>

#include "MAC_forward.hpp"


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
        bool is_mnp() const { return _len > 0 and _len == _seq_len; }
        bool is_del() const { return _len > 0 and _seq_len == 0; }
        bool is_empty() const { return _len == 0 and _seq_len == 0; }
        /**@}*/

        /** Cut mutation at given offsets.
         * @param base_offset Base offset, 0-based.
         * @param alt_offset Alternate sequence offset, 0-based.
         * @return The part of the original mutation that was cut from this object.
         */
        Mutation cut(Size_Type base_offset, Size_Type alt_offset);

        friend std::ostream& operator << (std::ostream&, const Mutation&);

    private:
        Seq_Type _seq;
        Size_Type _start;
        Size_Type _len;
        Size_Type _seq_len;
    };

    struct Mutation_Extra
    {
        const Mutation* mut_cptr;
        Size_Type alt_start;

        friend std::ostream& operator << (std::ostream&, const Mutation_Extra&);
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

        struct Mutation_Extra_Key
        {
            typedef Mutation::key_type result_type;
            result_type operator() (const Mutation_Extra& me) const { return me.mut_cptr->get_key(); }
        };
        struct Mutation_Extra_Alt_Key
        {
            typedef Mutation::key_type result_type;
            result_type operator() (const Mutation_Extra& me) const { return me.alt_start; }
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

    /** Container for Mutation_With_Alt_Pos objects. */
    typedef boost::multi_index_container<
      Mutation_Extra,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_non_unique<
          detail::Mutation_Extra_Key
        >,
        boost::multi_index::ordered_non_unique<
          detail::Mutation_Extra_Alt_Key
        >
      >
    > Mutation_Extra_Cont;

}

#endif
