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
    ~Allele() { ASSERT(is_unlinked()); }
public:
    GETTER(Seq_Type, idx, _idx)
private:
    Seq_Type _idx;
    /// Hooks for storage in Allele_Idx_Cont
    friend struct detail::Allele_Idx_List_Node_Traits;
    Allele_BPtr _previous;
    Allele_BPtr _next;
    bool is_unlinked() const { return not(_previous or _next); }
}; // class Allele_Idx

} // namespace MAC

#endif
