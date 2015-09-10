#ifndef __HAP_HOP_HPP
#define __HAP_HOP_HPP

#include "MAC_forward.hpp"
#include "Allele_Anchor.hpp"


namespace MAC
{

class Hap_Hop
{
private:
    // only constructed by Factory object
    friend class bounded::Factory< Hap_Hop >;

    /// Trivial default constructor
    DEFAULT_DEF_CTOR(Hap_Hop);

    // disallow copy or move
    DELETE_COPY_CTOR(Hap_Hop);
    DELETE_MOVE_CTOR(Hap_Hop);
    DELETE_COPY_ASOP(Hap_Hop);
    DELETE_MOVE_ASOP(Hap_Hop);

    Hap_Hop(Hap_Entry_CBPtr he_cbptr, const Allele_Anchor& allele_anchor,
            const Allele_Specifier& allele_specifier, bool c_direction)
        : _he_cbptr(he_cbptr), _allele_anchor(allele_anchor),
          _allele_specifier(allele_specifier), _c_direction(c_direction) {}

    ~Hap_Hop() { ASSERT(is_unlinked()); }

public:
    GETTER(Hap_Entry_CBPtr, he_cbptr, _he_cbptr)
    GETTER(Allele_Anchor, allele_anchor, _allele_anchor)
    GETTER(Allele_Specifier, allele_specifier, _allele_specifier)
    GETTER(bool, c_direction, _c_direction)
    Contig_Entry_CBPtr ce_cbptr() const { return _allele_anchor.ce_cbptr(); }

    ptree to_ptree() const
    {
        return ptree().put("he_cbptr", he_cbptr().to_int())
            .put("anchor", allele_anchor().to_ptree())
            .put("allele", allele_specifier().to_ptree())
            .put("c_dir", c_direction());
    }
    static void to_stream(ostream& os, Hap_Hop_CBPtr hh_cbptr)
    {
        os << "hop " << setw(5) << left << hh_cbptr.to_int()
           << "he " << setw(5) << left << hh_cbptr->he_cbptr().to_int()
           << "ce " << setw(5) << left << hh_cbptr->ce_cbptr().to_int();
        if (hh_cbptr->allele_anchor().is_endpoint())
        {
            os << "end " << hh_cbptr->allele_anchor().c_right()
               << " [ ce_next " << setw(5) << left << hh_cbptr->allele_specifier().ce_next_cbptr().to_int()
               << " " << hh_cbptr->allele_specifier().same_orientation() << " ]";
        }
        else
        {
            os << "mut " << hh_cbptr->allele_anchor().mut_cbptr().to_int()
               << " [ " << hh_cbptr->allele_specifier().strand() << " ]";
        }
    }

private:
    Hap_Entry_CBPtr _he_cbptr;

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
