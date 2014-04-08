//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __MAC_HPP
#define __MAC_HPP

#include "MAC_forward.hpp"
#include "Mutation.hpp"
#include "Mutation_Cont.hpp"
#include "Mutation_Ptr_Cont.hpp"
#include "Read_Chunk.hpp"
#include "Read_Chunk_Cont.hpp"
#include "Read_Chunk_Ptr_Cont.hpp"
#include "Read_Entry.hpp"
#include "Read_Entry_Cont.hpp"
#include "Contig_Entry.hpp"
#include "Contig_Entry_Cont.hpp"


namespace MAC
{

class Graph
{
public:
    /** Default constructor. */
    Graph() {}

    // disallow copy or move
    DELETE_COPY_CTOR(Graph)
    DELETE_MOVE_CTOR(Graph)
    DELETE_COPY_ASOP(Graph)
    DELETE_MOVE_ASOP(Graph)

    const Read_Entry_Cont& re_cont() const { return _re_cont; }
    const Contig_Entry_Cont& ce_cont() const { return _ce_cont; }

    /** Add a read.
     * Create basic Read_Entry, Read_Chunk, and Contig_Entry objects,
     * initialize them, and place them in their respective containers.
     * @param name_ptr Pointer to string with read name. (Read container takes ownership.)
     * @param seq_ptr Pointer to container with read sequence. (Contig container takes ownership.)
     */
    void add_read(const std::string* name_ptr, Seq_Type* seq_ptr);
    
    /** Add an overlap between 2 reads.
     * The contigs holding overlapping chunks of each read are collapsed into one.
     * In the process, reads and contigs are first fragmented into matching chunks.
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
    //void merge_all_read_contigs();
    
    /** Clear contig colors. */
    //void clear_contig_colours();
    //void clear_contig_visit();
    
    /** Mark contig endpoints and branches. */
    /*
    std::tuple< Size_Type, bool > visit_contig(const Contig_Entry* ce_cptr, bool dir);
    void dfs_scontig(const Contig_Entry* ce_cptr, bool ce_endpoint, bool dir,
                     std::list< std::tuple< const Contig_Entry*, size_t, size_t, bool > >& l, bool& cycle);
    void print_supercontig_lengths(std::ostream& os);
    void print_supercontig_lengths_2(std::ostream& os);
    void set_contig_ids();
    void unmap_single_chunks();
    */
    
    /** Integrity checks. */
    bool check_all() const;
    bool check(const std::set< Read_Entry_CBPtr >& re_set, const std::set< Contig_Entry_CBPtr >& ce_set = std::set< Contig_Entry_CBPtr >()) const;
    bool check(const std::set< Contig_Entry_CBPtr >& ce_set) const { return check(std::set< Read_Entry_CBPtr >(), ce_set); }
    //bool check_colours() const;
    
    /** Stats. */
    //void dump_detailed_counts(std::ostream& os) const;
    
    /*
    void print_separated_het_mutations(
        std::ostream& os, size_t min_support_report, Size_Type min_separation) const;
        void print_unmappable_contigs(std::ostream& os) const;
    */
        
    friend std::ostream& operator << (std::ostream&, const Graph&);
        
private:
    Mutation_Fact _mut_fact;
    Mutation_Chunk_Adapter_Fact _mca_fact;
    Read_Chunk_Fact _rc_fact;
    Read_Entry_Fact _re_fact;
    Contig_Entry_Fact _ce_fact;

    Read_Entry_Cont _re_cont;
    Contig_Entry_Cont _ce_cont;

    /** Cut Read_Entry at given read coordinate.
     * NOTE: A cut must be forced iff it is at the edge of a read.
     * @param re_bptr Read_Entry to cut.
     * @param r_brk Cut coordinate in read.
     * @param force False: enough to cut Read_Chunk only; True: always cut underlying Contrig_Entry.
     * @return True iff a cut was made.
     */
    bool cut_read_entry(Read_Entry_BPtr re_bptr, Size_Type r_brk, bool force);

    /** Cut Read_Chunk at given read coordinate.
     * NOTE: A cut must be forced iff it is at the edge of a read chunk.
     * @param rc_bptr Read_Chunk to cut.
     * @param r_brk Cut coordinate in read.
     * @param force True: always cut underlying Contrig_Entry.
     * @return True iff a cut was made.
     */
    bool cut_read_chunk(Read_Chunk_BPtr rc_bptr, Size_Type r_brk, bool force);

    /** Cut Contig_Entry object.
     * @param ce_bptr Contig_Entry to cut.
     * @param c_brk Contig position of the cut.
     * @param mut_left_cbptr If not NULL, single insertion at breakpoint
     * to be kept on the LHS (others are moved to RHS).
     */
    bool cut_contig_entry(Contig_Entry_BPtr ce_bptr, Size_Type c_brk, Mutation_CBPtr mut_left_cbptr);

    /*
    void erase_contig_entry(const Contig_Entry* ce_cptr);
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
    */
}; // class Graph

} // namespace MAC


#endif
