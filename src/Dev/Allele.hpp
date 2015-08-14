#ifndef __ALLELE_HPP
#define __ALLELE_HPP

#include "MAC_forward.hpp"

namespace MAC
{

class Allele
{
private:
    // only constructed by Factory object
    friend class bounded::Factory< Allele >;
    /// Trivial default constructor
    DEFAULT_DEF_CTOR(Allele);
    // disallow copy or move
    DELETE_COPY_CTOR(Allele);
    DELETE_MOVE_CTOR(Allele);
    DELETE_COPY_ASOP(Allele);
    DELETE_MOVE_ASOP(Allele);
    /// Constructor from seq.
    explicit Allele(Seq_Type&& seq) : _seq(move(seq)) {}
    /// Destructor: check if unlinked
    ~Allele() { ASSERT(is_unlinked()); }
public:
    GETTER(Seq_Type, seq, _seq)

    void save_strings(ostream& os, size_t& n_strings, size_t& n_bytes) const
    {
        os.write(_seq.c_str(), _seq.size() + 1);
        ++n_strings;
        n_bytes += _seq.size() + 1;
    }
    void init_strings()
    {
        new (&_seq) string();
    }
    void load_strings(istream& is, size_t& n_strings, size_t& n_bytes)
    {
        getline(is, _seq, '\0');
        ++n_strings;
        n_bytes += _seq.size() + 1;
    }

private:
    Seq_Type _seq;
    /// Hooks for storage in Allele_Cont
    friend struct detail::Allele_List_Node_Traits;
    Allele_BPtr _previous;
    Allele_BPtr _next;
    bool is_unlinked() const { return not(_previous or _next); }
}; // class Allele

} // namespace MAC

#endif
