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

#include "shortcuts.hpp"
#include "global.hpp"
#include "global_assert.hpp"
#include "factory.hpp"
#include "ptree.hpp"
#include "logger.hpp"
#include "Cigar.hpp"
//#include "Util.h"
#include "bits.hpp"
#include "ref_range.hpp"
#include "DNA_Sequence.hpp"
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
typedef dnasequence::Sequence< std::string, size_t > Seq_Type;
typedef Seq_Type::proxy_type Seq_Proxy_Type;
/** Type holding a single base pair. */
typedef Seq_Type::value_type Symb_Type;

/** Cigar object */
typedef cigar::Cigar< Size_Type > Cigar;

/** Range and Range_Cont */
typedef range::Range< Size_Type > Range_Type;
typedef range::Range_Cont< Size_Type > Range_Cont;

class Mutation;
class Mutation_Chunk_Adapter;
class Read_Chunk;
class Read_Entry;
class Contig_Entry;
class Graph;

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

// algorithms
class Unmap_Mut_Clusters;

} // namespace MAC


#endif
