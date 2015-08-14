//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __MUTATION_HPP
#define __MUTATION_HPP

#include "MAC_forward.hpp"
#include "Allele_Cont.hpp"


namespace MAC
{

/**
 * Holds information about a mutation from a base sequence.
 * The class holds the start and span of the base sequence region affected by a mutation,
 * as well as the alternate sequence.
 */
class Mutation
{
private:
    // Can only be created by Factory object
    friend class bounded::Factory< Mutation >;

    /// Default constructor.
    Mutation()
        : _start(0), _len(0) {}

    // disallow copy or move
    DELETE_COPY_CTOR(Mutation);
    DELETE_MOVE_CTOR(Mutation);
    DELETE_COPY_ASOP(Mutation);
    DELETE_MOVE_ASOP(Mutation);

    /**
     * Constructor from sequence.
     * @param start Start of the mutation, i.e., length of base sequence prior to mutation.
     * @param len Length of the base sequence affected by the mutation.
     * @param ce_cbptr Contig_Entry for this Mutation.
     * @param seq Alternate allele.
     */
    Mutation(Size_Type start, Size_Type len, Contig_Entry_CBPtr ce_cbptr, Seq_Type&& seq)
        : _start(start), _len(len), _ce_cbptr(ce_cbptr)
    {
        _alt_cont.push_back(move(seq));
    }
    Mutation(Size_Type start, Size_Type len, Contig_Entry_CBPtr ce_cbptr)
        : _start(start), _len(len), _ce_cbptr(ce_cbptr)
    {}

    ~Mutation()
    {
        ASSERT(is_unlinked());
        _alt_cont.clear_and_dispose();
    }

public:
    /// @name Getters
    /// @{
    Size_Type rf_start() const { return _start; }
    Size_Type rf_len() const { return _len; }
    Size_Type rf_end() const { return _start + _len; }
    GETTER(Contig_Entry_CBPtr, ce_cbptr, _ce_cbptr)

    unsigned n_alleles() const { return 1 + _alt_cont.size(); }

    /** Find equivalent alt, or add this one. */
    unsigned find_or_add_alt(const Seq_Proxy_Type& seq);
    /**
     * Remove an alternate allele.
     * After removal, the last alternate allele is moved to index i.
     * @return The index before this operation of the alt allele which is now moved at index i.
     */
    unsigned remove_alt(unsigned i);
    Size_Type allele_len(unsigned i) const
    {
        return (i == 0
                ? rf_len()
                : _alt_cont.at(i - 1).seq().size());
    }
    Seq_Proxy_Type allele_seq(unsigned i) const;

    /// @}

    /// @name Basic queries
    /// @{
    //bool is_ins() const { return _len == 0 and _seq_len > 0; }
    //bool is_snp() const { return _len == 1 and _len == _seq_len; }
    //bool is_del() const { return _len > 0 and _seq_len == 0; }
    //bool is_empty() const { return _len == 0 and _seq_len == 0; }
    /// @}

    /**
     * Merge with adjacent Mutation.
     * Pre: Both Mutations are unlinked.
     * @param rhs_mut_bptr Mutation to the right.
     * @param alt_map Map; keys are pairs of alt alleles to keep;
     * values are overwritten with their new alt index.
     */
    void merge_right(Mutation_BPtr rhs_mut_bptr, map< pair< unsigned, unsigned >, unsigned >& allele_map);

    /**
     * Cut mutation at given offsets, allocate new Mutation to keep leftover.
     * @param offset Offset, 0-based.
     * @return The second part of the Mutation that was cut;
     * and map where key = old alt index, value = new alt indexes.
     */
    pair< Mutation_BPtr, map< unsigned, pair< unsigned, unsigned > > >
    cut(Size_Type offset);

    /**
     * Simplify Mutation by dropping the ends of rf and alt alleles if they match.
     */
    void simplify();

    /**
     * Reverse the mutation.
     */
    void reverse();

    /**
     * Shift Mutation.
     * @param delta Signed integer value to add to start point.
     */
    template < typename delta_type >
    void shift(delta_type delta)
    {
        ASSERT(delta_type(_start) + delta >= 0);
        _start = Size_Type(delta_type(_start) + delta);
    }

    /*
    friend bool operator == (const Mutation& lhs, const Mutation& rhs)
    {
        return (lhs._start == rhs._start
                and lhs._len == rhs._len
                and lhs._seq_len == rhs._seq_len
                and lhs.have_seq() == rhs.have_seq()
                and (not lhs.have_seq() or lhs._seq == rhs._seq));
    }
    */
    friend bool operator < (const Mutation& lhs, const Mutation& rhs)
    {
        return lhs._start < rhs._start;
    }

    /*
    boost::property_tree::ptree to_ptree() const
    {
        return ptree().put("addr", static_cast< const void* >(this))
            .put("rf_start", rf_start())
            .put("rf_len", rf_len())
            .put("seq_len", seq_len())
            .put("seq", seq());
    }
    static void to_stream(ostream& os, Mutation_CBPtr mut_cbptr, Contig_Entry_CBPtr ce_cbptr);
    */

    void save_strings(ostream& os, size_t& n_strings, size_t& n_bytes) const;
    void init_strings();
    void load_strings(istream& is, size_t& n_strings, size_t& n_bytes);

private:
    // Hooks for storage in intrusive interval trees inside Contig_Entry objects.
    friend struct detail::Mutation_Set_Node_Traits;
    bool is_unlinked() const { return not(_parent or _l_child or _r_child); }

    Size_Type _start;
    Size_Type _len;
    Allele_Cont _alt_cont;
    Contig_Entry_CBPtr _ce_cbptr;
    Mutation_BPtr _parent;
    Mutation_BPtr _l_child;
    Mutation_BPtr _r_child;
    bool _color;
}; // class Mutation

struct Mutation_Comparator
{
    bool operator () (const Mutation& lhs, const Mutation& rhs) const { return lhs < rhs; }
    bool operator () (const Mutation& lhs, Size_Type rhs_rf_start) const { return lhs.rf_start() < rhs_rf_start; }
    bool operator () (Size_Type lhs_rf_start, const Mutation& rhs) const { return lhs_rf_start < rhs.rf_start(); }
    bool operator () (Size_Type lhs_rf_start, Size_Type rhs_rf_start) const { return lhs_rf_start < rhs_rf_start; }
}; // struct Mutation_Comparator

} // namespace MAC


#endif
