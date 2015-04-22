#include "Read_Chunk_Ptr_Cont.hpp"
#include "Read_Chunk.hpp"

namespace MAC
{

void Read_Chunk_Ptr_Cont::clear_and_dispose()
{
    Base::clear_and_dispose([] (Mutation_Chunk_Adapter_CBPtr mca_cbptr) {
            Read_Chunk_BPtr rc_bptr = mca_cbptr->chunk_cbptr().unconst();
            rc_bptr->mut_ptr_cont().erase(mca_cbptr);
            Mutation_Chunk_Adapter_Fact::del_elem(mca_cbptr);
        });
} // clear_and_dispose

} // namespace MAC
