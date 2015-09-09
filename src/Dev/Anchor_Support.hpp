#ifndef __ANCHOR_SUPPORT_HPP
#define __ANCHOR_SUPPORT_HPP

#include "Allele_Anchor.hpp"

namespace MAC
{

Anchor_Chunk_Support get_anchor_chunk_support(const Allele_Anchor& anchor);
Anchor_Read_Support get_anchor_read_support(const Anchor_Chunk_Support& anchor_chunk_support, bool c_direction);
Anchor_Read_Support get_anchor_read_support(const Allele_Anchor& anchor, bool c_direction);
Anchor_Read_Support get_hop_read_support(const Hap_Hop_CBPtr& hh_cbptr, bool h_direction);
Anchor_Read_Support get_hop_read_support(const pair< Hap_Hop_CBPtr, bool >& p);

Allele_Chunk_Support get_allele_chunk_support(const Anchor_Chunk_Support& anchor_chunk_support,
                                              const Allele_Specifier& allele);
Allele_Chunk_Support get_allele_chunk_support(Anchor_Chunk_Support&& anchor_chunk_support,
                                              const Allele_Specifier& allele);
Allele_Read_Support get_allele_read_support(const Anchor_Read_Support& anchor_read_support,
                                            const Allele_Specifier& allele);
Allele_Read_Support get_allele_read_support(Anchor_Read_Support&& anchor_read_support,
                                            const Allele_Specifier& allele);
Allele_Read_Support get_allele_read_support(const Allele_Anchor& anchor,
                                            const Allele_Specifier& allele,
                                            bool c_direction);

Allele_Read_Support collapsed_support(const Anchor_Read_Support& anchor_read_support);
void subset_support(Anchor_Read_Support& anchor_read_support, const Allele_Read_Support& common_support, bool direction);

Allele_Read_Support common_allele_support(const Allele_Read_Support& al1_read_support,
                                          const Allele_Read_Support& al2_read_support,
                                          bool direction);

} // namespace MAC

#endif
