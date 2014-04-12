//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __READ_ENTRY_HPP
#define __READ_ENTRY_HPP

#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <functional>

#include "MAC_forward.hpp"
#include "Read_Chunk.hpp"
#include "Read_Chunk_Cont.hpp"


namespace MAC
{

struct Read_Entry_Set_Node_Traits;

/** Holds information about a read.
 *
 * The read sequence gets assigned to various contigs. The Read_Entry object holds the read name,
 * as well as the (order) sequence of read chunks which are mapped to contigs.
 */
class Read_Entry
{
private:
    friend class Factory< Read_Entry >;

    /** Constructor.
     * @param name String containing read name; (take ownership).
     * @param len Length of the read.
     */
    Read_Entry(std::string&& name, Size_Type len) : _name(std::move(name)), _len(len)
    {
        Read_Chunk_BPtr rc_bptr = Read_Chunk_Fact::new_elem(nullptr, len);
        _chunk_cont.insert(rc_bptr);
    }

    // allow move only
    DEFAULT_DEF_CTOR(Read_Entry)
    DELETE_COPY_CTOR(Read_Entry)
    Read_Entry(Read_Entry&& rhs) { *this = std::move(rhs); }
public:
    DELETE_COPY_ASOP(Read_Entry)
    Read_Entry& operator = (Read_Entry&& rhs)
    {
        if (this != &rhs)
        {
            ASSERT(is_unlinked() and rhs.is_unlinked());
            _name = std::move(rhs._name);
            _chunk_cont = std::move(rhs._chunk_cont);
            _len = std::move(rhs._len);
        }
        return *this;
    }

    /** Find bounded pointer to this object.
     * Pre: Must be linked.
     */
    Read_Entry_CBPtr bptr_to() const
    {
        ASSERT(not is_unlinked());
        if (_parent->_l_child.raw() == this)
        {
            return _parent->_l_child;
        }
        if (_parent->_r_child.raw() == this)
        {
            return _parent->_r_child;
        }
        if (_parent->_parent.raw() == this)
        {
            return _parent->_parent;
        }
        ASSERT(false);
        return nullptr;
    }
    Read_Entry_BPtr bptr_to()
    {
        return boost::intrusive::pointer_traits< Read_Entry_BPtr >::const_cast_from(
            const_cast< const Read_Entry* >(this)->bptr_to());
    }

    /** @name Getters */
    /**@{*/
    const std::string& get_name() const { return _name; }
    Size_Type get_len() const { return _len; }
    const Read_Chunk_RE_Cont& chunk_cont() const { return _chunk_cont; }
    Read_Chunk_RE_Cont& chunk_cont() { return _chunk_cont; }
    Seq_Type get_seq() const;
    /**@}*/

    /** Check if this read ends the last contig where it is mapped.
     * @param check_start True to check read start, false to check read end.
     * @return True if there are no more bases in the contig past the read end.
     */
    bool is_terminal(bool check_start) const;

    /** Integrity check. */
    bool check() const;

    friend std::ostream& operator << (std::ostream& os, const Read_Entry& rhs);

private:
    std::string _name;
    Read_Chunk_RE_Cont _chunk_cont;
    Size_Type _len;

    /** Hooks for storage in intrusive set inside Graph object. */
    friend struct Read_Entry_Set_Node_Traits;
    Read_Entry_BPtr _parent;
    Read_Entry_BPtr _l_child;
    Read_Entry_BPtr _r_child;
    bool _col;
    bool is_unlinked() const { return not(_parent or _l_child or _r_child); }
}; // class Read_Entry

} // namespace MAC


#endif
