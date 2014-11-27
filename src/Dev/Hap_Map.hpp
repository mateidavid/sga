#ifndef __HAP_MAP_HPP
#define __HAP_MAP_HPP

#include "Graph.hpp"
#include "Hap_Entry_Cont.hpp"
#include "Hap_Hop_Set.hpp"


namespace MAC
{

class Hap_Map
{
public:
    Hap_Map(const Graph& g);

private:
    Hap_Entry_Fact _he_fact;
    Hap_Hop_Fact _hh_fact;

    Hap_Entry_Cont _he_cont;
    Hap_Hop_Set _hh_set;

}; // class Hap_Map

} // namespace MAC


#endif
