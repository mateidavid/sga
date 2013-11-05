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
    class Graph
    {
    public:
        const Read_Entry_Cont& re_cont() const { return _re_cont; }
        const Contig_Entry_Cont& ce_cont() const { return _ce_cont; }

        /** Retrieve Read_Entry corresponding to a given read.
         * @param name Read name.
         * @return CPtr to Read_Entry object, or NULL if not found.
         */
        const Read_Entry* get_read_entry(const std::string& name) const
        {
            Read_Entry_Cont::iterator it = _re_cont.find(name);
            return (it != _re_cont.end()? &(*it) : NULL);
        }

        /** Add a read.
         *
         * Create basic Read_Entry, Read_Chunk, and Contig_Entry objects,
         * initialize them, and place them in their respective containers.
         *
         * @param name_ptr Pointer to string with read name. (Read container takes ownership.)
         * @param seq_ptr Pointer to container with read sequence. (Contig container takes ownership.)
         */
        void add_read(const std::string* name_ptr, Seq_Type* seq_ptr);

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
        void add_overlap(const std::string& r1_name, const std::string& r2_name,
                         Size_Type r1_start, Size_Type r1_len,
                         Size_Type r2_start, Size_Type r2_len, bool r2_rc,
                         const std::string& cigar);

        /** Integrity check. */
        void check() const;

        friend std::ostream& operator << (std::ostream&, const Graph&);

    private:
        Read_Entry_Cont _re_cont;
        Contig_Entry_Cont _ce_cont;

        void modify_read_entry(const Read_Entry* re_cptr, Read_Entry::modifier_type modifier)
        {
            modify_element<Read_Entry_Cont>(_re_cont, re_cptr, modifier);
        }

        void modify_contig_entry(const Contig_Entry* ce_cptr, Contig_Entry::modifier_type modifier)
        {
            modify_element<Contig_Entry_Cont>(_ce_cont, ce_cptr, modifier);
        }

        void modify_read_chunk(Read_Chunk_CPtr rc_cptr, Read_Chunk::modifier_type modifier)
        {
            modify_read_entry(rc_cptr->get_re_ptr(), [&] (Read_Entry& re) { re.modify_read_chunk(rc_cptr, modifier); });
        }

        const Read_Entry* insert_read_entry(const Read_Entry& re);
        const Contig_Entry* insert_contig_entry(const Contig_Entry& ce);
        void cut_read_entry(const Read_Entry* re_cptr, Size_Type r_brk, bool force = false);
        void cut_read_chunk(Read_Chunk_CPtr rc_cptr, Size_Type r_brk);
        void cut_mutation(const Contig_Entry* ce_cptr, const Mutation* mut_cptr, Size_Type c_offset, Size_Type r_offset);
        void cut_contig_entry(const Contig_Entry* ce_cptr, Size_Type c_brk, const Mutation* mut_left_cptr);
    };

    std::vector<std::pair< Read_Chunk, Mutation_Cont >> make_chunks_from_cigar(
        Size_Type r1_start, Size_Type r1_len,
        Size_Type r2_start, Size_Type r2_len, bool r2_rc,
        const std::string& cigar);

}


#endif
