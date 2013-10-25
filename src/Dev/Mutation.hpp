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
        : _start(0), _len(0) {}

        /** Constructor.
         * @param start Start of the mutation, i.e., length of base sequence prior to mutation.
         * @param len Length of the base sequence affected by the mutation.
         * @param seq Alternate sequence.
         */
        Mutation(Size_Type start, Size_Type len, const Seq_Type& seq = Seq_Type())
        : _seq(seq), _start(start), _len(len) {}

        /** Copy constructor. */
        Mutation(const Mutation& rhs)
        : _seq(rhs._seq), _start(rhs._start), _len(rhs._len) {}

        /** Static insertion constructor.
         * @param start Start of the mutation, i.e., length of base sequence prior to mutation.
         * @param seq Inserted sequence.
         */
        static Mutation make_ins(Size_Type start, const Seq_Type& seq)
        {
            assert(seq.size() > 0);
            return Mutation(start, 0, seq);
        }

        /** Static snp constructor.
         * @param start Start of the mutation, i.e., length of base sequence prior to mutation.
         * @param symbol Alternate symbol.
         */
        static Mutation make_snp(Size_Type start, const Seq_Type& seq)
        {
            assert(seq.size() == 1);
            return Mutation(start, 1, seq);
        }

        /** Static deletion constructor.
         * @param start Start of the mutation, i.e., length of base sequence prior to mutation.
         */
        static Mutation make_del(Size_Type start)
        {
            return Mutation(start, 1);
        }

        /** @name Getters */
        /**@{*/
        Size_Type get_start() const { return _start; }
        Size_Type get_len() const { return _len; }
        Size_Type get_end() const { return _start+_len; }
        const Seq_Type& get_seq() const { return _seq; }
        key_type get_key() const { return _start; }
        /**@}*/

        /** @name Basic queries */
        /**@{*/
        bool is_ins() const { return _len == 0 and _seq.size() == 1; }
        bool is_snp() const { return _len == 1 and _seq.size() == 1; }
        bool is_del() const { return _len == 1 and _seq.size() == 0; }
        bool is_empty() const { return _len == 0 and _seq.size() == 0; }
        /**@}*/

        friend std::ostream& operator << (std::ostream&, const Mutation&);

    private:
        Seq_Type _seq;
        Size_Type _start;
        Size_Type _len;
    };

    namespace detail
    {
        /** Key extractor struct for boost::multi_index_container. */
        struct Mutation_Key
        {
            typedef Mutation::key_type result_type;
            result_type operator() (const Mutation& m) const { return m.get_key(); }
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
}

#endif
