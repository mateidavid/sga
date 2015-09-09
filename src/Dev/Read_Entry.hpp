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

    /// Trivial default constructor
    DEFAULT_DEF_CTOR(Read_Entry);

    // disallow copy or move
    DELETE_COPY_CTOR(Read_Entry);
    DELETE_MOVE_CTOR(Read_Entry);
    DELETE_COPY_ASOP(Read_Entry);
    DELETE_MOVE_ASOP(Read_Entry);

    /**
     * Constructor.
     * @param name String containing read name; (take ownership).
     * @param len Length of the read.
     */
    Read_Entry(string&& name, Size_Type len) : _name(move(name)), _orig_len(len), _len(len), _start(0) {}

    ~Read_Entry()
    {
        ASSERT(chunk_cont().empty());
        ASSERT(is_unlinked());
    }

public:
    /** @name Getters */
    /**@{*/
    GETTER(string, name, _name)
    GETTER(Read_Chunk_RE_Cont, chunk_cont, _chunk_cont)
    GETTER(Size_Type, len, _len)
    GETTER(Size_Type, start, _start)
    Size_Type end() const { return start() + len(); }
    Seq_Type get_seq() const;
    /**@}*/

    /**
     * Check if this read ends the last contig where it is mapped.
     * @param check_start True to check read start, false to check read end.
     * @return True if there are no more bases in the contig past the read end.
     */
    bool is_terminal(bool check_start) const;

    /**
     * Edit read sequence.
     * This operation irreversibly changes the read sequence of this entry
     * by altering its implicit representation.
     * @param rc_bptr Chunk where the edit occurs.
     * @param start Absolute read position of the edit (it must be inside rc_bptr).
     * @param old_seq Old read sequence which is being deleted.
     * @param new_seq New read sequence which is being inserted.
     */
    void edit(Read_Chunk_BPtr rc_bptr, Size_Type start, const Seq_Proxy_Type& old_seq, const Seq_Proxy_Type& new_seq)
    {
        LOG("Read_Entry", info) << ptree("begin")
            .put("re_name", name())
            .put("start", start)
            .put("old_seq", old_seq)
            .put("new_seq", new_seq);
        Size_Type old_seq_len = old_seq.size();
        Size_Type new_seq_len = new_seq.size();
        ptrdiff_t delta = static_cast< ptrdiff_t >(new_seq_len) - static_cast< ptrdiff_t >(old_seq_len);
        ASSERT(rc_bptr->get_r_start() <= start);
        ASSERT(start + old_seq_len <= rc_bptr->get_r_end());
        chunk_cont().implement_edit(rc_bptr, delta);
        _len = static_cast< ptrdiff_t >(_len) + delta;
    }

    /** Integrity check. */
    void check() const;

    //friend std::ostream& operator << (std::ostream& os, const Read_Entry& rhs);
    boost::property_tree::ptree to_ptree() const;

    void save_strings(ostream& os, size_t& n_strings, size_t& n_bytes) const;
    void init_strings();
    void load_strings(istream& is, size_t& n_strings, size_t& n_bytes);

private:
    string _name;
    Read_Chunk_RE_Cont _chunk_cont;
    Size_Type _orig_len;
    Size_Type _len;
    Size_Type _start;

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
