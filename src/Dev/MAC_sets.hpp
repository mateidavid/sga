#ifndef __MAC_SETS_HPP
#define __MAC_SETS_HPP

#include "Allele_Specifier.hpp"
#include "RC_DSet.hpp"
#include "RE_DSet.hpp"

namespace MAC
{

typedef RC_DSet Allele_Chunk_Support;
typedef map< Allele_Specifier, Allele_Chunk_Support > Anchor_Chunk_Support;
typedef RE_DSet Allele_Read_Support;
typedef map< Allele_Specifier, Allele_Read_Support > Anchor_Read_Support;

} // namespace MAC

#endif
