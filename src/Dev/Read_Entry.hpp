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

namespace detail
{

struct Read_Entry_Set_Node_Traits;

}

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
    Read_Entry(string&& name, Size_Type len) : _name(move(name)), _len(len) {}

    // allow move only
    DEFAULT_DEF_CTOR(Read_Entry)
    DELETE_COPY_CTOR(Read_Entry)
    Read_Entry(Read_Entry&& rhs) { *this = move(rhs); }

    ~Read_Entry()
    {
        ASSERT(chunk_cont().empty());
        ASSERT(is_unlinked());
    }

public:
    DELETE_COPY_ASOP(Read_Entry)
    Read_Entry& operator = (Read_Entry&& rhs)
    {
        if (this != &rhs)
        {
            ASSERT(is_unlinked() and rhs.is_unlinked());
            _name = move(rhs._name);
            _chunk_cont = move(rhs._chunk_cont);
            _len = move(rhs._len);
        }
        return *this;
    }

    /** @name Getters */
    /**@{*/
    GETTER(string, name, _name)
    GETTER(Read_Chunk_RE_Cont, chunk_cont, _chunk_cont)
    Size_Type len() const { return _len; }
    Seq_Type get_seq() const;
    /**@}*/

    /** Check if this read ends the last contig where it is mapped.
     * @param check_start True to check read start, false to check read end.
     * @return True if there are no more bases in the contig past the read end.
     */
    bool is_terminal(bool check_start) const;

    /** Integrity check. */
    void check() const;

    //friend std::ostream& operator << (std::ostream& os, const Read_Entry& rhs);
    boost::property_tree::ptree to_ptree() const;

private:
    string _name;
    Read_Chunk_RE_Cont _chunk_cont;
    Size_Type _len;

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
