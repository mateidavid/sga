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
 * Special case: the special allele (nullptr, true) is used to denote that
 * the contig has out degree 0.
 */
class Allele_Specifier
{
public:
    DEFAULT_DEF_CTOR(Allele_Specifier);

    Allele_Specifier(Contig_Entry_CBPtr ce_next_cbptr, bool same_orientation)
        : _ce_next_cbptr(ce_next_cbptr), _same_orientation(same_orientation), _strand(false) {}

    explicit Allele_Specifier(const pair< Contig_Entry_CBPtr, bool >& p)
        : _ce_next_cbptr(p.first), _same_orientation(p.second), _strand(false) {}

    explicit Allele_Specifier(bool strand)
        : _ce_next_cbptr(), _same_orientation(false), _strand(strand) {}

    bool is_mutation() const { return not is_endpoint(); }
    bool is_endpoint() const { return _ce_next_cbptr or _same_orientation; }

    GETTER(Contig_Entry_CBPtr, ce_next_cbptr, _ce_next_cbptr)
    GETTER(bool, same_orientation, _same_orientation)
    GETTER(bool, strand, _strand)

    friend bool operator == (const Allele_Specifier& lhs, const Allele_Specifier& rhs)
    {
        return (lhs._ce_next_cbptr == rhs._ce_next_cbptr)
            and (lhs.is_endpoint()? lhs._same_orientation == rhs._same_orientation : lhs._strand == rhs._strand);
    }
    friend bool operator != (const Allele_Specifier& lhs, const Allele_Specifier& rhs) { return !(lhs == rhs); }
    friend bool operator <  (const Allele_Specifier& lhs, const Allele_Specifier& rhs)
    {
        if (lhs.is_endpoint())
        {
            if (rhs.is_endpoint())
            {
                return make_pair(lhs.ce_next_cbptr(), lhs.same_orientation())
                    < make_pair(rhs.ce_next_cbptr(), rhs.same_orientation());
            }
            else
            {
                return true;
            }
        }
        else
        {
            if (rhs.is_endpoint())
            {
                return false;
            }
            else
            {
                return lhs.strand() < rhs.strand();
            }
        }
    }
    friend bool operator <= (const Allele_Specifier& lhs, const Allele_Specifier& rhs) { return lhs == rhs or lhs < rhs; }
    friend bool operator >  (const Allele_Specifier& lhs, const Allele_Specifier& rhs) { return !(lhs <= rhs); }
    friend bool operator >= (const Allele_Specifier& lhs, const Allele_Specifier& rhs) { return !(lhs < rhs); }

    boost::property_tree::ptree to_ptree() const
    {
        if (is_endpoint())
        {
            return ptree().put("ce_next_cbptr", ce_next_cbptr().to_int())
                .put("same_orientation", same_orientation());
        }
        else
        {
            return ptree().put("strand", _strand);
        }
    }

private:
    Contig_Entry_CBPtr _ce_next_cbptr;
    bool _same_orientation;
    bool _strand;
}; // class Allele_Specifier

} // namespace MAC


#endif
