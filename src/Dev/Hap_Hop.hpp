#ifndef __HAP_HOP_HPP
#define __HAP_HOP_HPP

#include "MAC_forward.hpp"
#include "Allele_Anchor.hpp"


namespace MAC
{

namespace detail
{

struct Hap_Hop_List_Node_Traits;
struct Hap_Hop_Set_Node_Traits;

}

class Hap_Hop
{
private:
    // only constructed by Factory object
    friend class bounded::Factory< Hap_Hop >;

    // disallow copy or move
    DEFAULT_DEF_CTOR(Hap_Hop)
    DELETE_COPY_CTOR(Hap_Hop)
    DELETE_MOVE_CTOR(Hap_Hop)
    DELETE_COPY_ASOP(Hap_Hop)
    DELETE_MOVE_ASOP(Hap_Hop)

    Hap_Hop(const Allele_Anchor& allele_anchor, const Allele_Specifier& allele_specifier, bool c_direction)
    : _allele_anchor(allele_anchor), _allele_specifier(allele_specifier), _c_direction(c_direction) {}

    ~Hap_Hop() { ASSERT(is_unlinked()); }

public:
    GETTER(Allele_Anchor, allele_anchor, _allele_anchor)
    GETTER(Allele_Specifier, allele_specifier, _allele_specifier)
    GETTER(bool, c_direction, _c_direction)
    Contig_Entry_CBPtr ce_cbptr() const { return _allele_anchor.ce_cbptr(); }

private:
    Hap_Entry_BPtr _hap_bptr;

    /** Hooks for storage in intrusive list. */
    friend struct detail::Hap_Hop_List_Node_Traits;
    Hap_Hop_BPtr _previous;
    Hap_Hop_BPtr _next;

    /** Hooks for storage in intrusive multiset. */
    friend struct detail::Hap_Hop_Set_Node_Traits;
    Hap_Hop_BPtr _parent;
    Hap_Hop_BPtr _l_child;
    Hap_Hop_BPtr _r_child;

    Allele_Anchor _allele_anchor;
    Allele_Specifier _allele_specifier;
    bool _c_direction;
    bool _col;

    bool is_unlinked() const { return not(_previous or _next or _parent or _l_child or _r_child); }
}; // class Hap_Hop

} // namespace MAC


#endif
