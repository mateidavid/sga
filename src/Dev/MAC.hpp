//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __MAC_HPP
#define __MAC_HPP

#include "MAC_forward.hpp"
#include "Mutation.hpp"
#include "Read_Chunk.hpp"
#include "Read_Entry.hpp"
#include "Contig_Entry.hpp"


namespace MAC
{
    /** Add a read.
     *
     * Create basic Read_Entry, Read_Chunk, and Contig_Entry objects,
     * initialize them, and place them in their respective containers.
     *
     * @param name_ptr Pointer to string with read name. (Read container takes ownership.)
     * @param seq_ptr Pointer to container with read sequence. (Contig container takes ownership.)
     * @param re_cont Container for read entry objects.
     * @param ce_cont Container for contig entry objects.
     */
    void add_read(Read_Entry_Cont& re_cont, Contig_Entry_Cont& ce_cont,
                  const std::string* name_ptr, const Seq_Type* seq_ptr);

    std::vector<std::pair< Read_Chunk, Mutation_Cont >> make_chunks_from_cigar(
        Size_Type r1_start, Size_Type r1_len,
        Size_Type r2_start, Size_Type r2_len, bool r2_rc,
        const std::string& cigar);

    /** Add an overlap between 2 reads.
     *
     * The contigs holding overlapping chunks of each read are collapsed into one.
     * In the process, reads and contigs are first fragmented into matching chunks.
     *
     * @param r1_name Name of read 1.
     * @param r2_name Name of read 2.
     * @param r1_start Overlap start in r1, 0-based. Equivalently, length of r1 before the overlap.
     * @param r1_len Overlap length in r1.
     * @param r2_start Overlap start in r2, 0-based. Equivalently, length of r2 before the overlap.
     * @param r2_len Overlap length in r2.
     * @param r2_rc True iff r2 is reverse complemented.
     * @param cigar Cigar string (r1:reference, r2:query).
     */
    void add_overlap(Read_Entry_Cont& re_cont, Contig_Entry_Cont& ce_cont,
                     const std::string& r1_name, const std::string& r2_name,
                     Size_Type r1_start, Size_Type r1_len,
                     Size_Type r2_start, Size_Type r2_len, bool r2_rc,
                     const std::string& cigar);
}


#endif
