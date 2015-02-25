#ifndef __HAP_ENTRY_HPP
#define __HAP_ENTRY_HPP

#include "MAC_forward.hpp"
#include "Hap_Hop_List.hpp"


namespace MAC
{

namespace detail
{

struct Hap_Entry_List_Node_Traits;

}

class Hap_Entry
{
private:
    // only constructed by Factory object
    friend class bounded::Factory< Hap_Entry >;

    // disallow copy or move
    DEFAULT_DEF_CTOR(Hap_Entry);
    DELETE_COPY_CTOR(Hap_Entry);
    DELETE_MOVE_CTOR(Hap_Entry);
    DELETE_COPY_ASOP(Hap_Entry);
    DELETE_MOVE_ASOP(Hap_Entry);

    ~Hap_Entry() { ASSERT(is_unlinked()); }

public:
    GETTER(Hap_Hop_List, hh_cont, _hh_cont)

private:
    Hap_Hop_List _hh_cont;

    /** Hooks for storage in intrusive list. */
    friend struct detail::Hap_Entry_List_Node_Traits;
    Hap_Entry_BPtr _previous;
    Hap_Entry_BPtr _next;

    bool is_unlinked() const { return not(_previous or _next); }
}; // class Hap_Entry

} // namespace MAC


#endif
