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
#include "Range_Cont.hpp"


namespace MAC
{

class Graph
{
public:
    /** Default constructor. */
    Graph() : _unmap_trigger_len(9) {}

    // disallow copy or move
    DELETE_COPY_CTOR(Graph)
    DELETE_MOVE_CTOR(Graph)
    DELETE_COPY_ASOP(Graph)
    DELETE_MOVE_ASOP(Graph)

    ~Graph()
    {
        ASSERT(ce_cont().empty());
        ASSERT(re_cont().empty());
    }

    const Read_Entry_Cont& re_cont() const { return _re_cont; }
    Read_Entry_Cont& re_cont() { return _re_cont; }
    const Contig_Entry_Cont& ce_cont() const { return _ce_cont; }
    Contig_Entry_Cont& ce_cont() { return _ce_cont; }
    const Size_Type& unmap_trigger_len() const { return _unmap_trigger_len; }
    Size_Type& unmap_trigger_len() { return _unmap_trigger_len; }

    /** Add a read.
     * Create basic Read_Entry, Read_Chunk, and Contig_Entry objects,
     * initialize them, and place them in their respective containers.
     * @param name_ptr String with read name. (Read container takes ownership.)
     * @param seq_ptr Read sequence. (Contig container takes ownership.)
     */
    void add_read(string&& name, Seq_Type&& seq);

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
     * @param cat_at_step Catenate contigs after adding the overlap.
     */
    void add_overlap(const string& r1_name, const string& r2_name,
                     Size_Type r1_start, Size_Type r1_len,
                     Size_Type r2_start, Size_Type r2_len, bool r2_rc,
                     const string& cigar, bool cat_at_step);

    /** Merge all contigs. */
    void cat_all_read_contigs();

    /** Get reads beyond which coverage drops (that is, supercontigs end with out-degree zero).
     * For each such read, we print the strand with which the coverage ends.
     * Thus, to extend the current supercontig, one should find mappings beyond the 3' end
     * of the sequences included in the output. A read's positive and negative strand may both
     * appear in the output.
     * @param os Output stream.
     */
    void get_terminal_reads(ostream& os) const;

    /** Clear contig colors. */
    //void clear_contig_colours();
    //void clear_contig_visit();

    /** Mark contig endpoints and branches. */
    /*
    tuple< Size_Type, bool > visit_contig(const Contig_Entry* ce_cptr, bool dir);
    void dfs_scontig(const Contig_Entry* ce_cptr, bool ce_endpoint, bool dir,
                     list< tuple< const Contig_Entry*, size_t, size_t, bool > >& l, bool& cycle);
    void print_supercontig_lengths(ostream& os);
    void print_supercontig_lengths_2(ostream& os);
    */

    /** Unmap read chunks not mapped to lone contigs. */
    void unmap_single_chunks();

    /** Unmap read ends not mapped to anything else. */
    void unmap_read_ends();

    /** Clear CE and RE containers and deallocate all objects. */
    void clear_and_dispose();

    /** Integrity checks. */
    void check_all() const;
    void check(const set< Read_Entry_CBPtr >& re_set,
               const set< Contig_Entry_CBPtr >& ce_set = set< Contig_Entry_CBPtr >()) const;
    void check(const set< Contig_Entry_CBPtr >& ce_set) const
    {
        check(set< Read_Entry_CBPtr >(), ce_set);
    }
    void check_leaks() const;

    /** Stats. */
    void dump_detailed_counts(ostream& os) const;

    /*
    void print_separated_het_mutations(
        ostream& os, size_t min_support_report, Size_Type min_separation) const;
    */

    /** Print unmappable contigs between each pair of mappable ones. */
    void print_unmappable_contigs(ostream& os) const;

    /** Resolve unmappable contigs. */
    void resolve_unmappable_regions();
    void resolve_unmappable_inner_region(Contig_Entry_CBPtr ce_cbptr, bool c_right,
                                         Contig_Entry_CBPtr ce_next_cbptr, bool same_orientation);
    void resolve_unmappable_terminal_region(Contig_Entry_CBPtr ce_cbptr, bool c_right);

    /** Assign unmappable read sequences to a set of base sequences.
     * Base sequences are picked from the set of read sequences, the rest of read sequences
     * are mapped to these base sequences.
     * @param seq_cnt_map Read sequences, with count of appearances.
     * @param bseq_v Vector of base sequences.
     * @param seq_bseq_map Mapping (read sequence) -> (index of bseq, cigar_string).
     */
    void resolve_unmappable_fully_mapped(
        const map< string, size_t >& seq_cnt_map,
        vector< string >& bseq_v,
        map< string, tuple< size_t, string > >& seq_bseq_map);
    /** Remap set of unmappable seq ends to a set of base sequences.
     * @param seq_set Set of pairs (sequence, side):
     * the mapping of sequence to abse sequences should be anchored on side.
     * @param bseq_v Vector of base sequences to choose from
     * @param seq_bseq_map Mapping (seq, side) -> (index of bseq, cigar_string).
     */
    void resolve_unmappable_partially_mapped(
        const set< tuple< string, bool > >& seq_set,
        const vector< string >& bseq_v,
        map< tuple< string, bool >, tuple< size_t, string > >& seq_bseq_map);

    /** Process interactive commands. */
    void interactive_commands(istream& is, ostream& os);

    //friend ostream& operator << (ostream&, const Graph&);
    boost::property_tree::ptree to_ptree() const;
    boost::property_tree::ptree factory_stats() const;

    void save(ostream&) const;
    void load(istream&);

private:
    Mutation_Fact _mut_fact;
    Mutation_Chunk_Adapter_Fact _mca_fact;
    Read_Chunk_Fact _rc_fact;
    Read_Entry_Fact _re_fact;
    Contig_Entry_Fact _ce_fact;

    Read_Entry_Cont _re_cont;
    Contig_Entry_Cont _ce_cont;

    Size_Type _unmap_trigger_len;

    /** Cut Read_Entry at given read coordinate.
     * NOTE: A cut must be forced iff it is at the edge of a read.
     * @param re_bptr Read_Entry to cut.
     * @param r_brk Cut coordinate in read.
     * @return True iff a cut was made.
     */
    bool cut_read_entry(Read_Entry_BPtr re_bptr, Size_Type r_brk);

    /** Cut Read_Chunk at given read coordinate.
     * NOTE: A cut must be forced iff it is at the edge of a read chunk.
     * @param rc_bptr Read_Chunk to cut.
     * @param r_brk Cut coordinate in read.
     * @return True iff a cut was made.
     */
    bool cut_read_chunk(Read_Chunk_BPtr rc_bptr, Size_Type r_brk);

    /** Cut Contig_Entry object.
     * @param ce_bptr Contig_Entry to cut.
     * @param c_brk Contig position of the cut.
     * @param mut_left_cbptr If not NULL, single insertion at breakpoint
     * to be kept on the LHS (others are moved to RHS).
     */
    bool cut_contig_entry(Contig_Entry_BPtr ce_bptr, Size_Type c_brk, Mutation_CBPtr mut_left_cbptr);

    /** Repeatedly cut read entries and match as described in the cigar string.
     * @param re1_bptr Read_Entry reference.
     * @param re2_bptr Read_Entry query.
     * @param cigar Cigar string describing the mapping.
     * @return Vector of tuples of the form (chunk r_start, chunk r_start, cigar).
     */
    vector< tuple< Size_Type, Size_Type, Cigar > >
    chunker(Read_Entry_BPtr re1_bptr, Read_Entry_BPtr re2_bptr, Cigar& cigar);

    /** Make chunk unmappable.
     * The contig containing the chunk is first trimmed to the extent of the chunk.
     * After potential contig trimming, all other chunks in the contig now fully spanned
     * by this chunk are also made unmappable.
     */
    void unmap_chunk(Read_Chunk_BPtr rc_bptr);

    /** Unmap read regions.
     * The read is cut at region endpoints, and all chunks spanning the given regions
     * are made unmappable.
     */
    void unmap_regions(Read_Entry_BPtr re_bptr, const Range_Cont< Size_Type >& region_cont);

    /** Try to extend an unmappable read region in both directions.
     * @param re_bptr Read_Entry object.
     * @param rc_start Extend left of this position on read.
     * @param rc_end Extend right of this position on read.
     */
    void extend_unmapped_chunk(Read_Entry_BPtr re_bptr, Size_Type rc_start, Size_Type rc_end);

    /** Extend an unmappable read region in the given direction.
     * @param re_bptr Read_Entry object.
     * @param pos Read position where to extend from.
     * @param r_right Bool; true: extend right on read; false: extend left on read.
     */
    void extend_unmapped_chunk_dir(Read_Entry_BPtr re_bptr, Size_Type pos, bool r_right);

    /** Attempt to catenate contig with neighbour in given direction.
     * @param ce_bptr Contig_Entry to consider.
     * @param c_right Bool; true: merge past contig end; false: merge past contig start.
     * @return True iff the merge was successful.
     */
    bool cat_contigs(Contig_Entry_BPtr ce_bptr, bool c_right);

    /** Merge Contig_Entry objects of 2 read chunks according to the cigar mapping of the chunks.
     * NOTE: If chunks are already mapped to the same contig, nothing is changed.
     * Pre: Chunks must be mapped to the full length of their contigs.
     * Post: All chunks from c2 (including rc2) are remapped to c1;
     * new mutations are added to c1 as needed;
     * rc2 is mapped as follows: c1<-rc1<-rc2;
     * all other chunks rcX are mapped indirectly through c2 as follows: c1<-rc1<-rc2<-c2<-rcX;
     * c2 and all its chunks, mutations, and MCAs, are deallocated.
     * @param c1rc1_chunk_bptr Mapping of rc1 to c1.
     * @param c2rc2_chunk_bptr Mapping of rc2 to c2.
     * @param rc1rc2_cigar Mapping of rc2 to rc1.
     */
    void merge_chunk_contigs(Read_Chunk_BPtr c1rc1_chunk_bptr, Read_Chunk_BPtr c2rc2_chunk_bptr, Cigar& rc1rc2_cigar);

    /** Find unmappable regions in a read entry.
     * @param re_cbptr Read to scan.
     * @param r_start Read position to scan from.
     * @param r_end Read position to scan to.
     * @return A set of disjoint ranges which should be made unmappable in this read.
     */
    Range_Cont< Size_Type >
    find_unmappable_regions(Read_Entry_CBPtr re_cbptr, Size_Type r_start, Size_Type r_end) const;

    /** Find unmappable regions in a read chunk.
     * @param rc_cbptr Read chunk to scan.
     * @param region_cont Set of ranges of the chunk's read which should be made unmappable.
     */
    void find_unmappable_regions(Read_Chunk_CBPtr rc_cbptr, Range_Cont< Size_Type >& region_cont) const;

    void cat_read_contigs(Read_Entry_BPtr re_bptr);

    /** Unmap terminal read chunk if it is mapped to contig end with no other supporting chunk.
     * @param rc_bptr Chunk to check.
     * @param r_start Bool; true: unmap in the direction of read start; false: unmap in the direction of read end.
     */
    void unmap_single_terminal_chunk(Read_Chunk_BPtr rc_bptr, bool r_start);
}; // class Graph

} // namespace MAC


#endif
