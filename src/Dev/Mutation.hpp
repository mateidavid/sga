//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __MUTATION_HPP
#define __MUTATION_HPP

#include "MAC_forward.hpp"
#include "Read_Chunk_Ptr_Cont.hpp"


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
        : _start(0), _len(0), _seq_len(0) {}

    /**
     * Constructor from sequence.
     * @param start Start of the mutation, i.e., length of base sequence prior to mutation.
     * @param len Length of the base sequence affected by the mutation.
     * @param seq Alternate sequence.
     */
    Mutation(Size_Type start, Size_Type len, Seq_Type&& seq)
        : _seq(move(seq)), _start(start), _len(len), _seq_len(_seq.size()) {}

    /**
     * Constructor from sequence length.
     * @param start Start of the mutation, i.e., length of base sequence prior to mutation.
     * @param len Length of the base sequence affected by the mutation.
     * @param seq_len Length of alternate sequence.
     */
    Mutation(Size_Type start, Size_Type len, Size_Type seq_len = 0)
        : _start(start), _len(len), _seq_len(seq_len) {}

    DELETE_COPY_CTOR(Mutation);
    Mutation(Mutation&& rhs) : Mutation() { *this = move(rhs); }
    DELETE_COPY_ASOP(Mutation);

    ~Mutation()
    {
        ASSERT(chunk_ptr_cont().empty());
        ASSERT(is_unlinked());
    }

public:
    Mutation& operator = (Mutation&& rhs)
    {
        if (this != &rhs)
        {
            // allow move only when unlinked
            ASSERT(is_unlinked() and rhs.is_unlinked());
            _seq = move(rhs._seq);
            _start = move(rhs._start);
            _len = move(rhs._len);
            _seq_len = move(rhs._seq_len);
            _chunk_ptr_cont = move(rhs._chunk_ptr_cont);
        }
        return *this;
    }

    /// @name Getters
    /// @{
    Size_Type rf_start() const { return _start; }
    Size_Type rf_len() const { return _len; }
    Size_Type rf_end() const { return _start + _len; }
    Size_Type seq_len() const { return _seq_len; }
    bool have_seq() const { return _seq.size() == _seq_len; }
    GETTER(Seq_Type, seq, _seq)
    /// @}

    /// @name Basic queries
    /// @{
    bool is_ins() const { return _len == 0 and _seq_len > 0; }
    bool is_snp() const { return _len == 1 and _len == _seq_len; }
    bool is_del() const { return _len > 0 and _seq_len == 0; }
    bool is_empty() const { return _len == 0 and _seq_len == 0; }
    /// @}

    /**
     * Extend Mutation given sequence.
     * Pre: This Mutation is unlinked.
     * Pre: Mutation contains its alternate sequence.
     * @param start Reference position.
     * @param extra_len Extra reference length.
     * @param extra_seq Extra alternate sequence.
     */
    void extend(Size_Type start, Size_Type extra_len, const Seq_Proxy_Type& extra_seq)
    {
        ASSERT(is_unlinked());
        ASSERT(have_seq());
        if (is_empty())
        {
            _start = start;
        }
        ASSERT(rf_end() == start);
        _len += extra_len;
        _seq += extra_seq;
        _seq_len += extra_seq.size();
    }

    /**
     * Extend Mutation given sequence size.
     * Pre: This Mutation is unlinked.
     * Pre: Mutation contains its alternate sequence.
     * @param start Reference position.
     * @param extra_len Extra reference length.
     * @param extra_seq Extra alternate sequence.
     */
    void extend(Size_Type start, Size_Type extra_len, Size_Type extra_seq_size)
    {
        ASSERT(is_unlinked());
        ASSERT(_seq.empty());
        if (is_empty())
        {
            _start = start;
        }
        ASSERT(rf_end() == start);
        _len += extra_len;
        _seq_len += extra_seq_size;
    }

    /**
     * Extend Mutation; if empty, copy the given Mutation.
     * Pre: This Mutation is unlinked.
     * Pre: Both Mutations contain alternate sequences.
     * @param extra_mut_cbptr Extra mutation.
     */
    void extend(Mutation_CBPtr extra_mut_cbptr)
    {
        ASSERT(is_unlinked());
        ASSERT(have_seq() and extra_mut_cbptr->have_seq());
        if (is_empty())
        {
            _start = extra_mut_cbptr->rf_start();
        }
        ASSERT(rf_end() == extra_mut_cbptr->rf_start());
        _len += extra_mut_cbptr->rf_len();
        _seq += extra_mut_cbptr->seq();
        _seq_len += extra_mut_cbptr->seq_len();
    }

    /**
     * Cut mutation at given offsets, allocate new Mutation to keep leftover.
     * @param base_offset Base offset, 0-based.
     * @param alt_offset Alternate sequence offset, 0-based.
     * @return The second part of the Mutation that was cut.
     */
    Mutation_CBPtr cut(Size_Type base_offset, Size_Type alt_offset);

    /**
     *Simplify Mutation by dropping the ends of rf and qr if they match.
     * @param rf Reference sequence spanned by the mutation.
     */
    void simplify(const Seq_Proxy_Type& rf);

    /**
     * Reverse the mutation.
     * @param c_len The contig length.
     */
    void reverse(Size_Type c_len)
    {
        _start = c_len - (_start + _len);
        if (have_seq())
        {
            _seq = _seq.revcomp();
        }
    }

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

    const Read_Chunk_Ptr_Cont& chunk_ptr_cont() const { return _chunk_ptr_cont; }
    Read_Chunk_Ptr_Cont& chunk_ptr_cont() { return _chunk_ptr_cont; }

    friend bool operator == (const Mutation& lhs, const Mutation& rhs)
    {
        return (lhs._start == rhs._start
                and lhs._len == rhs._len
                and lhs._seq_len == rhs._seq_len
                and lhs.have_seq() == rhs.have_seq()
                and (not lhs.have_seq() or lhs._seq == rhs._seq));
    }

    boost::property_tree::ptree to_ptree() const
    {
        return ptree().put("addr", static_cast< const void* >(this))
            .put("rf_start", rf_start())
            .put("rf_len", rf_len())
            .put("seq_len", seq_len())
            .put("seq", seq());
    }
    static void to_stream(ostream& os, Mutation_CBPtr mut_cbptr, Contig_Entry_CBPtr ce_cbptr);

    void save_strings(ostream& os, size_t& n_strings, size_t& n_bytes) const;
    void init_strings();
    void load_strings(istream& is, size_t& n_strings, size_t& n_bytes);

private:
    // Hooks for storage in intrusive interval trees inside Contig_Entry objects.
    friend struct detail::Mutation_ITree_Node_Traits;
    bool is_unlinked() const { return not(_parent or _l_child or _r_child); }

    Seq_Type _seq;
    Size_Type _start;
    Size_Type _col_n_max_end;
    Read_Chunk_Ptr_Cont _chunk_ptr_cont;
    Mutation_BPtr _parent;
    Mutation_BPtr _l_child;
    Mutation_BPtr _r_child;
    uint32_t _len;
    uint32_t _seq_len;
}; // class Mutation

} // namespace MAC


#endif
