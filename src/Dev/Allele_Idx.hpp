#ifndef __ALLELE_IDX_HPP
#define __ALLELE_IDX_HPP

#include "MAC_forward.hpp"

namespace MAC
{

class Allele_Idx
{
private:
    // only constructed by Factory object
    friend class bounded::Factory< Allele_Idx >;
    /// Trivial default constructor
    DEFAULT_DEF_CTOR(Allele_Idx);
    // disallow copy or move
    DELETE_COPY_CTOR(Allele_Idx);
    DELETE_MOVE_CTOR(Allele_Idx);
    DELETE_COPY_ASOP(Allele_Idx);
    DELETE_MOVE_ASOP(Allele_Idx);
    /// Constructor from seq.
    explicit Allele_Idx(unsigned idx) : _idx(idx) {}
    /// Destructor: check if unlinked
    ~Allele_Idx() { ASSERT(is_unlinked()); }
public:
    GETTER(unsigned, idx, _idx)
    boost::property_tree::ptree to_ptree() const { return ptree().put("idx", idx()); }
private:
    unsigned _idx;
    /// Hooks for storage in Allele_Idx_Cont
    friend struct detail::Allele_Idx_List_Node_Traits;
    Allele_Idx_BPtr _previous;
    Allele_Idx_BPtr _next;
    bool is_unlinked() const { return not(_previous or _next); }
}; // class Allele_Idx

} // namespace MAC

#endif
