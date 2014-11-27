#include "Hap_Map.hpp"


namespace MAC
{

Hap_Map::Hap_Map(const Graph& g)
{
    build(g);
}

void Hap_Map::build(const Graph& g)
{
    for (auto ce_cbptr : g.ce_cont())
    {
        if (not ce_cptr->is_normal()) continue;

        //TODO
    }
}

} // namespace MAC
