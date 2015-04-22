//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __READ_ENTRY_HPP
#define __READ_ENTRY_HPP

#include "MAC_forward.hpp"
#include "Read_Chunk.hpp"
#include "Read_Chunk_Cont.hpp"


namespace MAC
{

/** Holds information about a read.
 *
 * The read sequence gets assigned to various contigs. The Read_Entry object holds the read name,
 * as well as the (order) sequence of read chunks which are mapped to contigs.
 */
class Read_Entry
{
private:
    friend class bounded::Factory< Read_Entry >;

    /** Constructor.
     * @param name String containing read name; (take ownership).
     * @param len Length of the read.
     */
    Read_Entry(string&& name, Size_Type len) : _name(move(name)), _orig_len(len), _len(len), _start(0), _delta(0) {}

    // allow move only
    DEFAULT_DEF_CTOR(Read_Entry);
    DELETE_COPY_CTOR(Read_Entry);
    Read_Entry(Read_Entry&& rhs) { *this = move(rhs); }

    ~Read_Entry()
    {
        ASSERT(chunk_cont().empty());
        ASSERT(is_unlinked());
    }

public:
    DELETE_COPY_ASOP(Read_Entry);
    Read_Entry& operator = (Read_Entry&& rhs)
    {
        if (this != &rhs)
        {
            ASSERT(is_unlinked() and rhs.is_unlinked());
            _name = move(rhs._name);
            _chunk_cont = move(rhs._chunk_cont);
            _len = move(rhs._len);
            _start = move(rhs._start);
        }
        return *this;
    }

    /** @name Getters */
    /**@{*/
    GETTER(string, name, _name)
    GETTER(Read_Chunk_RE_Cont, chunk_cont, _chunk_cont)
    private:
    GETTER(Size_Type, len, _len)
    public:
    GETTER(Size_Type, start, _start)
    Size_Type end() const { return start() + len(); }
    //Size_Type len() const { return _len; }
    Seq_Type get_seq(bool trimmed) const;
    Size_Type get_len(bool trimmed) const
    {
        return trimmed? len() : len() + _start_seq.size() + _end_seq.size();
    }
    /**@}*/

    /** Check if this read ends the last contig where it is mapped.
     * @param check_start True to check read start, false to check read end.
     * @return True if there are no more bases in the contig past the read end.
     */
    bool is_terminal(bool check_start) const;

    /**
     * Trim read entry.
     * @param r_end Bool; if true, trim from the end; if false, trim from the start.
     * @param s Sequence that was trimmed.
     */
    void trim_end(bool r_end, const Seq_Proxy_Type& s)
    {
        Size_Type trim_len = s.size();
        ASSERT(len() >= trim_len);
        if (r_end)
        {
            Seq_Type tmp(s);
            tmp += _end_seq;
            swap(tmp, _end_seq);
            _len -= trim_len;
        }
        else
        {
            _start_seq += s;
            _len -= trim_len;
            _start += trim_len;
        }
    }

    void add_edit(Size_Type start, Size_Type len, const Seq_Proxy_Type& seq)
    {
        (void)start;
        (void)len;
        (void)seq;
        ptrdiff_t delta = static_cast< ptrdiff_t >(len) - seq.size();
        _len = delta + _len;
        _delta += delta;
    }

    /** Integrity check. */
    void check() const;

    //friend std::ostream& operator << (std::ostream& os, const Read_Entry& rhs);
    boost::property_tree::ptree to_ptree() const;

private:
    friend class Graph; // for saving&loading _start_seq and _end_seq
    string _name;
    Seq_Type _start_seq;
    Seq_Type _end_seq;
    Read_Chunk_RE_Cont _chunk_cont;
    Size_Type _orig_len;
    Size_Type _len;
    Size_Type _start;
    ptrdiff_t _delta;

    /** Hooks for storage in intrusive set inside Graph object. */
    friend struct detail::Read_Entry_Set_Node_Traits;
    Read_Entry_BPtr _parent;
    Read_Entry_BPtr _l_child;
    Read_Entry_BPtr _r_child;
    bool _col;
    bool is_unlinked() const { return not(_parent or _l_child or _r_child); }
}; // class Read_Entry

} // namespace MAC


#endif
