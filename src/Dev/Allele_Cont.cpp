#include "Allele_Cont.hpp"

namespace MAC
{

void Allele_Cont::save_strings(ostream& os, size_t& n_strings, size_t& n_bytes) const
{
    for (auto al_cbptr : *this | referenced)
    {
        al_cbptr->save_strings(os, n_strings, n_bytes);
    }
}

void Allele_Cont::init_strings()
{
    Base::get_root_node()->init_strings();
    for (auto al_bptr : *this | referenced)
    {
        al_bptr->init_strings();
    }
}

void Allele_Cont::load_strings(istream& is, size_t& n_strings, size_t& n_bytes)
{
    for (auto al_bptr : *this | referenced)
    {
        al_bptr->load_strings(is, n_strings, n_bytes);
    }
}

} // namespace MAC
