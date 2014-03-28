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
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>

#include "MAC_forward.hpp"
#include "Read_Chunk.hpp"


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
     * @param name_ptr Pointer to string containing read name.
     * @param len Length of the read.
     */
    Read_Entry(const std::string* name_ptr, Size_Type len) : _name_ptr(name_ptr), _len(len)
    {
        ASSERT(name_ptr != NULL);
        Read_Chunk_BPtr rc_bptr = Read_Chunk_Fact::new_elem(nullptr, len);
        _chunk_cont.insert(rc_bptr);
    }

    // allow move only
    DELETE_COPY_CTOR(Read_Entry)
    Read_Entry(Read_Entry&& rhs) { *this = std::move(rhs); }
public:
    DELETE_COPY_ASOP(Read_Entry)
    Read_Entry& operator = (Read_Entry&& rhs)
    {
        if (this != &rhs)
        {
            ASSERT(is_unlinked() and rhs.is_unlinked());
            _name_ptr = std::move(rhs._name_ptr);
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
        return static_cast< Read_Entry_BPtr >(const_cast< const Read_Entry* >(this)->bptr_to());
    }

    /** @name Getters */
    /**@{*/
    const std::string& get_name() const { return *_name_ptr; }
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
    std::shared_ptr< const std::string > _name_ptr;
    Read_Chunk_RE_Cont _chunk_cont;
    Size_Type _len;

    /** Hooks for storage in intrusive set inside Graph object. */
    friend struct Read_Entry_Set_Node_Traits;
    Read_Entry_BPtr _parent;
    Read_Entry_BPtr _l_child;
    Read_Entry_BPtr _r_child;
    bool _col;
    bool is_unlinked() const { return not(_parent or _l_child or _r_child); }
};

/** Comparator for storage in tree. */
struct Read_Entry_Compare
{
    bool operator () (const Read_Entry& lhs, const Read_Entry& rhs) const
    {
        return lhs.get_name() < rhs.get_name();
    }
};

struct Read_Entry_Set_Node_Traits
{
    typedef Holder< Read_Entry > node;
    typedef Read_Entry_Fact fact_type;
    typedef fact_type::ptr_type node_ptr;
    typedef fact_type::const_ptr_type const_node_ptr;
    typedef int color;

    static node_ptr get_parent(const_node_ptr n) { return n->_parent; }
    static void set_parent(node_ptr n, node_ptr ptr) { n->_parent = ptr; }
    static node_ptr get_left(const_node_ptr n) { return n->_l_child; }
    static void set_left(node_ptr n, node_ptr ptr) { n->_l_child = ptr; }
    static node_ptr get_right(const_node_ptr n) { return n->_r_child; }
    static void set_right(node_ptr n, node_ptr ptr) { n->_r_child = ptr; }
    static color get_color(const_node_ptr n) { return n->_col; }
    static void set_color(node_ptr n, color c) { n->_col = c ; }
    static color black() { return 0; }
    static color red() { return 1; }
};

struct Read_Entry_Set_Value_Traits
{
    typedef Read_Entry value_type;
    typedef Read_Entry_Set_Node_Traits node_traits;
    typedef node_traits::node_ptr node_ptr;
    typedef node_traits::const_node_ptr const_node_ptr;
    typedef node_ptr pointer;
    typedef const_node_ptr const_pointer;
    typedef node_traits::fact_type::ref_type reference;
    typedef node_traits::fact_type::const_ref_type const_reference;

    static const boost::intrusive::link_mode_type link_mode = boost::intrusive::normal_link;

    static node_ptr to_node_ptr (reference value) { return &value; }
    static const_node_ptr to_node_ptr (const_reference value) { return &value; }
    static pointer to_value_ptr(node_ptr n) { return n; }
    static const_pointer to_value_ptr(const_node_ptr n) { return n; }
};

typedef boost::intrusive::set<
            Read_Entry,
            boost::intrusive::compare< Read_Entry_Compare >,
            boost::intrusive::value_traits< Read_Entry_Set_Value_Traits >
            > Read_Entry_Cont;
} // namespace MAC


#endif
