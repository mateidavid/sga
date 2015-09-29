#include <boost/preprocessor.hpp>

#include "Mutation.hpp"
#include "Mutation_Chunk_Adapter.hpp"
#include "Read_Chunk.hpp"
#include "Read_Entry.hpp"
#include "Contig_Entry.hpp"
//#include "Hap_Hop.hpp"
//#include "Hap_Entry.hpp"

#define BOUNDED_TYPE(I) MAC::BOUNDED_TYPE_##I
#define BOUNDED_TYPE_0 Mutation
#define BOUNDED_TYPE_1 Mutation_Chunk_Adapter
#define BOUNDED_TYPE_2 Read_Chunk
#define BOUNDED_TYPE_3 Read_Entry
#define BOUNDED_TYPE_4 Contig_Entry
//#define BOUNDED_TYPE_5 Hap_Hop
//#define BOUNDED_TYPE_6 Hap_Entry
#define BOUNDED_TYPE_NUM 5

#define INSTANTIATE_BOUNDED_TYPE(Z, I, _) \
    template union bounded::detail::Value_Wrapper< BOUNDED_TYPE(I) >; \
    template class bounded::detail::Storage< BOUNDED_TYPE(I) >; \
    template struct bounded::detail::Identifier< BOUNDED_TYPE(I) >; \
    template class bounded::Pointer< BOUNDED_TYPE(I) >; \
    template class bounded::Pointer< const BOUNDED_TYPE(I) >; \
    template class bounded::Reference< BOUNDED_TYPE(I) >; \
    template class bounded::Reference< const BOUNDED_TYPE(I) >; \
    template class bounded::Factory< BOUNDED_TYPE(I) >;

namespace bounded
{

BOOST_PP_REPEAT(BOUNDED_TYPE_NUM, INSTANTIATE_BOUNDED_TYPE, _)

} // namespace bounded
