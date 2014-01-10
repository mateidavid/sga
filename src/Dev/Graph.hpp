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

        /** Merge all contigs. */
        void merge_all_read_contigs();

        /** Integrity checks. */
        bool check_all() const;
        bool check(const std::set< const Read_Entry* >& re_set, const std::set< const Contig_Entry* >& ce_set = std::set< const Contig_Entry* >()) const;
        bool check(const std::set< const Contig_Entry* >& ce_set) const { return check(std::set< const Read_Entry* >(), ce_set); }

        /** Stats. */
        void dump_detailed_counts(std::ostream& os) const;

        friend std::ostream& operator << (std::ostream&, const Graph&);

    private:
        Read_Entry_Cont _re_cont;
        Contig_Entry_Cont _ce_cont;

        void modify_read_entry(const Read_Entry* re_cptr, Read_Entry::mod_type re_mod)
        {
            modify_element<Read_Entry_Cont>(_re_cont, re_cptr, re_mod);
        }

        void modify_contig_entry(const Contig_Entry* ce_cptr, Contig_Entry::mod_type ce_mod)
        {
            modify_element<Contig_Entry_Cont>(_ce_cont, ce_cptr, ce_mod);
        }

        void modify_read_chunk(Read_Chunk_CPtr rc_cptr, Read_Chunk::mod_type rc_mod)
        {
            modify_read_entry(rc_cptr->get_re_ptr(), [&] (Read_Entry& re) { re.modify_read_chunk(rc_cptr, rc_mod); });
        }

        const Read_Entry* insert_read_entry(const Read_Entry& re);
        const Contig_Entry* insert_contig_entry(const Contig_Entry& ce);
        void erase_contig_entry(const Contig_Entry* ce_cptr);
        bool cut_read_entry(const Read_Entry* re_cptr, Size_Type r_brk, bool force = false);
        bool cut_read_chunk(Read_Chunk_CPtr rc_cptr, Size_Type r_brk, bool force = false);
        void cut_mutation(const Contig_Entry* ce_cptr, const Mutation* mut_cptr, Size_Type c_offset, Size_Type r_offset);
        bool cut_contig_entry(const Contig_Entry* ce_cptr, Size_Type c_brk, const Mutation* mut_left_cptr);

        void remap_chunks(std::map< Read_Chunk_CPtr, std::shared_ptr< Read_Chunk > >& rc_map, Mutation_Cont& extra_mut_cont);
        void merge_read_chunks(Read_Chunk_CPtr c1rc1_chunk_cptr, Read_Chunk_CPtr c2rc2_chunk_cptr, Cigar& rc1rc2_cigar);
        std::shared_ptr< std::vector< std::tuple< Read_Chunk_CPtr, Read_Chunk_CPtr, Cigar > > > chunker(
            const Read_Entry* re1_cptr, const Read_Entry* re2_cptr, Cigar& cigar);
        void reverse_contig(const Contig_Entry* ce_cptr);
        bool try_merge_contig(const Contig_Entry* ce_cptr, bool forward);
        void merge_read_contigs(const Read_Entry* re_cptr);
        void unmap_chunk(Read_Chunk_CPtr rc_cptr);
        void extend_unmapped_chunk(const Read_Entry* re_cptr, Size_Type rc_start, Size_Type rc_end);
        void extend_unmapped_chunk_dir(const Read_Entry* re_cptr, Size_Type pos, bool dir);
        void merge_unmappable_chunks(Read_Chunk_CPtr rc1_cptr, Read_Chunk_CPtr rc2_cptr);
        void scan_read_for_unmappable_chunks(const Read_Entry* re_cptr, Size_Type rc_start, Size_Type rc_end);
        void scan_contig_for_unmappable_chunks(const Contig_Entry* ce_cptr);
    };
}


#endif
