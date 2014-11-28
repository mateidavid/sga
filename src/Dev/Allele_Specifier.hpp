#ifndef __ALLELE_SPECIFIER_HPP
#define __ALLELE_SPECIFIER_HPP

#include "MAC_forward.hpp"

namespace MAC
{

/** Allele specifier type.
 * When the anchor is a mutation, an allele is specified by a bool:
 * rf=false; qr=true.
 * When the anchor is an endpoint, an allele is specified by the
 * destination of the edge: (ce_next, same_orientation)
 */
class Allele_Specifier
{
public:
    DEFAULT_DEF_CTOR(Allele_Specifierv)
    DEFAULT_COPY_CTOR(Allele_Specifier)
    DEFAULT_MOVE_CTOR(Allele_Specifier)
    DEFAULT_COPY_ASOP(Allele_Specifier)
    DEFAULT_MOVE_ASOP(Allele_Specifier)

    Allele_Specifier(Contig_Entry_CBPtr ce_cbptr, bool same_orientation)
        : _ce_cbptr(ce_cbptr), _same_orientation(same_orientation) {}

    explicit Allele_Specifier(const pair< Contig_Entry_CBPtr, bool >& p)
        : _ce_cbptr(p.first), _same_orientation(p.second) {}

    explicit Allele_Specifier(bool strand)
        : _strand(strand) {}

    bool is_mutation() const { return not is_endpoint(); }
    bool is_endpoint() const { return _ce_cbptr; }

    GETTER(Contig_Entry_CBPtr, ce_cbptr, _ce_cbptr)
    GETTER(bool, same_orientation, _same_orientation)
    GETTER(bool, strand, _strand)

private:
    Contig_Entry_CBPtr _ce_cbptr;
    bool _same_orientation;
    bool _strand;
}; // class Allele_Specifier

} // namespace MAC


#endif
