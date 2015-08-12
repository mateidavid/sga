//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __MAC_FORWARD_HPP
#define __MAC_FORWARD_HPP

#include <tuple>
#include <string>
#include <list>
#include <vector>
#include <map>
#include <set>
#include <memory>
#include <functional>
#include <utility>
#include <iostream>
#include <boost/intrusive/list.hpp>
#include <boost/intrusive/set.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/range/join.hpp>
#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include "shortcuts.hpp"
#include "global_assert.hpp"
#include "factory.hpp"
#include "ptree.hpp"
#include "logger.hpp"
#include "Cigar.hpp"
#include "bits.hpp"
#include "ref_range.hpp"
#include "RC_Sequence.hpp"
#include "Range.hpp"
#include "Range_Cont.hpp"


#include <boost/intrusive/itree.hpp>


/** Multi-Allelic Contig namespace. */
namespace MAC
{

using namespace std;
namespace bi = boost::intrusive;
namespace ba = boost::adaptors;
using boost::range::join;
using boost::optional;


/** Type for absolute and relative offsets in read and contig sequences. */
typedef size_t Size_Type;
/** Type holding a read sequence. */
typedef rc_sequence::Sequence< std::string > Seq_Type;
typedef Seq_Type::proxy_type Seq_Proxy_Type;
/** Type holding a single base pair. */
typedef Seq_Type::value_type Symb_Type;

/** Cigar object */
typedef cigar::Cigar< Size_Type > Cigar;

/** Range and Range_Cont */
typedef range::Range< Size_Type > Range_Type;
typedef range::Range_Cont< Size_Type > Range_Cont;

class Allele;
class Allele_Idx;
class Mutation;
class Mutation_Chunk_Adapter;
class Read_Chunk;
class Read_Entry;
class Contig_Entry;
class Hap_Hop;
class Hap_Entry;
class Graph;

namespace detail
{
    struct Allele_List_Node_Traits;
    struct Allele_Idx_List_Node_Traits;
    struct Mutation_ITree_Node_Traits;
    struct Mutation_Ptr_List_Node_Traits;
    struct Read_Chunk_Ptr_List_Node_Traits;
    struct Read_Chunk_ITree_Node_Traits;
    struct Read_Chunk_Set_Node_Traits;
    struct Read_Entry_Set_Node_Traits;
    struct Contig_Entry_List_Node_Traits;
    struct Hap_Entry_List_Node_Traits;
    struct Hap_Hop_List_Node_Traits;
    struct Hap_Hop_Set_Node_Traits;
} // namespace detail

typedef bounded::Factory< Allele > Allele_Fact;
typedef Allele_Fact::ptr_type Allele_BPtr;
typedef Allele_Fact::const_ptr_type Allele_CBPtr;
typedef Allele_Fact::ref_type Allele_BRef;
typedef Allele_Fact::const_ref_type Allele_CBRef;

typedef bounded::Factory< Allele_Idx > Allele_Idx_Fact;
typedef Allele_Idx_Fact::ptr_type Allele_Idx_BPtr;
typedef Allele_Idx_Fact::const_ptr_type Allele_Idx_CBPtr;
typedef Allele_Idx_Fact::ref_type Allele_Idx_BRef;
typedef Allele_Idx_Fact::const_ref_type Allele_Idx_CBRef;

typedef bounded::Factory< Mutation > Mutation_Fact;
typedef Mutation_Fact::ptr_type Mutation_BPtr;
typedef Mutation_Fact::const_ptr_type Mutation_CBPtr;
typedef Mutation_Fact::ref_type Mutation_BRef;
typedef Mutation_Fact::const_ref_type Mutation_CBRef;

typedef bounded::Factory< Mutation_Chunk_Adapter > Mutation_Chunk_Adapter_Fact;
typedef Mutation_Chunk_Adapter_Fact::ptr_type Mutation_Chunk_Adapter_BPtr;
typedef Mutation_Chunk_Adapter_Fact::const_ptr_type Mutation_Chunk_Adapter_CBPtr;
typedef Mutation_Chunk_Adapter_Fact::ref_type Mutation_Chunk_Adapter_BRef;
typedef Mutation_Chunk_Adapter_Fact::const_ref_type Mutation_Chunk_Adapter_CBRef;

typedef bounded::Factory< Read_Chunk > Read_Chunk_Fact;
typedef Read_Chunk_Fact::ptr_type Read_Chunk_BPtr;
typedef Read_Chunk_Fact::const_ptr_type Read_Chunk_CBPtr;
typedef Read_Chunk_Fact::ref_type Read_Chunk_BRef;
typedef Read_Chunk_Fact::const_ref_type Read_Chunk_CBRef;

typedef bounded::Factory< Read_Entry > Read_Entry_Fact;
typedef Read_Entry_Fact::ptr_type Read_Entry_BPtr;
typedef Read_Entry_Fact::const_ptr_type Read_Entry_CBPtr;
typedef Read_Entry_Fact::ref_type Read_Entry_BRef;
typedef Read_Entry_Fact::const_ref_type Read_Entry_CBRef;

typedef bounded::Factory< Contig_Entry > Contig_Entry_Fact;
typedef Contig_Entry_Fact::ptr_type Contig_Entry_BPtr;
typedef Contig_Entry_Fact::const_ptr_type Contig_Entry_CBPtr;
typedef Contig_Entry_Fact::ref_type Contig_Entry_BRef;
typedef Contig_Entry_Fact::const_ref_type Contig_Entry_CBRef;

typedef bounded::Factory< Hap_Hop > Hap_Hop_Fact;
typedef Hap_Hop_Fact::ptr_type Hap_Hop_BPtr;
typedef Hap_Hop_Fact::const_ptr_type Hap_Hop_CBPtr;
typedef Hap_Hop_Fact::ref_type Hap_Hop_BRef;
typedef Hap_Hop_Fact::const_ref_type Hap_Hop_CBRef;

typedef bounded::Factory< Hap_Entry > Hap_Entry_Fact;
typedef Hap_Entry_Fact::ptr_type Hap_Entry_BPtr;
typedef Hap_Entry_Fact::const_ptr_type Hap_Entry_CBPtr;
typedef Hap_Entry_Fact::ref_type Hap_Entry_BRef;
typedef Hap_Entry_Fact::const_ref_type Hap_Entry_CBRef;

// algorithms
class Unmap_Mut_Clusters;
class Hap_Map;

} // namespace MAC


#endif
