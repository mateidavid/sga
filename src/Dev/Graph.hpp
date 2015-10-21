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
#include "Allele_Anchor.hpp"
#include "Unmapper.hpp"

#include "BWT.h"
#include "BWTAlgorithms.h"

namespace MAC
{

class Graph
{
public:
    /** Default constructor. */
    Graph()
        : _u{this},
        _unmap_trigger_len(9),
        _aux_coverage(-1),
        _cat_at_step(false),
        _trim_tuc_step(false) {}

    // disallow copy or move
    DELETE_COPY_CTOR(Graph);
    DELETE_MOVE_CTOR(Graph);
    DELETE_COPY_ASOP(Graph);
    DELETE_MOVE_ASOP(Graph);

    ~Graph()
    {
        ASSERT(ce_cont().empty());
        ASSERT(re_cont().empty());
        if (_index_set.pBWT) delete _index_set.pBWT;
        if (_index_set.pSSA) delete _index_set.pSSA;
        if (_aux_index_set.pBWT) delete _aux_index_set.pBWT;
    }

    const Read_Entry_Cont& re_cont() const { return _re_cont; }
    Read_Entry_Cont& re_cont() { return _re_cont; }
    const Contig_Entry_Cont& ce_cont() const { return _ce_cont; }
    Contig_Entry_Cont& ce_cont() { return _ce_cont; }
    const Unmapper& unmapper() const { return _u; }
    const Size_Type& unmap_trigger_len() const { return _unmap_trigger_len; }
    Size_Type& unmap_trigger_len() { return _unmap_trigger_len; }
    GETTER(bool, cat_at_step, _cat_at_step)
    GETTER(bool, trim_tuc_step, _trim_tuc_step)
    GETTER(BWTIndexSet, index_set, _index_set)
    GETTER(BWTIndexSet, aux_index_set, _aux_index_set)
    int get_aux_coverage() const { return _aux_coverage; }

    /** Add a read.
     * Create basic Read_Entry, Read_Chunk, and Contig_Entry objects,
     * initialize them, and place them in their respective containers.
     * @param name_ptr String with read name. (Read container takes ownership.)
     * @param seq_ptr Read sequence. (Contig container takes ownership.)
     */
    void add_read(string&& name, Seq_Type&& seq);

    /** Remove a read */
    void remove_read(Read_Entry_BPtr re_bptr);

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
                     const string& cigar);

    /** Merge all contigs. */
    void cat_all_read_contigs();

    /** Trim terminal unmappable chunks. */
    void trim_tucs();

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

    typedef list< pair< Contig_Entry_CBPtr, bool > > supercontig_type;
    typedef list< supercontig_type > supercontig_list;
    /// Return a list of supercontigs in the graph
    /// Each supercontig is a list of (contig, orientation) pairs
    /// Also sets Contig_Entry::supercontig_endpoint_mask bits in the contig tags
    supercontig_list get_supercontigs(int unmappable_policy, size_t ignore_threshold = 1) const;
    Size_Type skip_supercontig_bulges(supercontig_list & l) const;

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
    void print_basic_stats(ostream& os) const;
    void print_detailed_counts(ostream& os) const;
    void print_supercontig_stats(ostream& os) const;
    void print_mutations(ostream& os, size_t min_support = 1, Size_Type flank_len = 20) const;

    /** Print unmappable contigs between each pair of mappable ones. */
    void print_unmappable_contigs(ostream& os) const;

    void test_mutation_allele_swapping();

    /// Load BWT, SSA and read id list.
    void load_bwt(const string& bwt_prefix);

    /// Load BWT of Illumina data.
    void load_aux_bwt(const string& aux_bwt_file);

    typedef array< vector< pair< string, Size_Type > >, 2 > find_reads_with_seq_type;
    /// Find reads containing a given sequence
    find_reads_with_seq_type find_reads_with_seq(const Seq_Proxy_Type& seq, unsigned max_count = 0) const;

    typedef pair< array< vector< pair< Read_Entry_CBPtr, unsigned > >, 2 >, bool > find_read_entries_with_seq_type;
    /**
     * Find REs containing a given sequence.
     * @return A pair (vector< re, pos >[2], flag). flag is true iff
     * there are reads containing seq in the BWT but not in the graph.
     */
    find_read_entries_with_seq_type
    find_read_entries_with_seq(const Seq_Proxy_Type& seq, unsigned max_count = 0) const;

    /** Split Read_Entry after given Read_Chunk in the given read direction. */
    pair< Read_Entry_CBPtr, Read_Entry_CBPtr >
    split_read(Read_Chunk_CBPtr rc_cbptr);
    /** Split RE into 2 overlapping REs, duplicating the sequence between l_anchor and its right neighbour. */
    pair< Read_Entry_CBPtr, Read_Entry_CBPtr >
    split_read(Read_Chunk_CBPtr rc_cbptr, const Allele_Anchor& l_anchor);

    void compute_unique_flank_len(Contig_Entry_BPtr ce_bptr, Mutation_BPtr mut_bptr, bool allele,
                                  Size_Type min_flank_len, Size_Type max_flank_len);
    void compute_unique_flank_lens(Size_Type min_flank_len, Size_Type max_flank_len);
    void compute_mutation_copy_num();
    void compute_aux_coverage();

    void print_contig(ostream& os, Contig_Entry_CBPtr ce_cbptr, bool c_direction,
                      map< Read_Chunk_CBPtr, unsigned >& rc_grid_pos) const;

    ptree to_ptree() const;
    ptree factory_stats() const;

    void save(ostream&) const;
    void load(istream&);
    void export_gfa(ostream&, bool) const;

private:
    friend class Unmap_Mut_Clusters;
    friend class Unmapper;
    friend class Read_Merger;

    /**
     * Cut Read_Entry at given read coordinate.
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

    /**
     * Trim contig entry to the extent of the given chunk.
     * Return the chunk, now spanning the entire contig.
     */
    //Read_Chunk_BPtr trim_contig_to_chunk(Read_Chunk_BPtr rc_bptr);

    /**
     * Trim terminal unmappable chunk.
     */
    void trim_tuc(Read_Chunk_BPtr rc_bptr);

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

    /** Find overlapping mappings in a read entry.
     * @param re_cbptr RE to scan.
     * @param r_start Read position start.
     * @param r_end Read position end.
     * @return A set of disjoint ranges which contain overlapping mappings (of this or other reads).
     */
    Range_Cont
    find_overlapping_regions(Read_Entry_CBPtr re_cbptr, Size_Type r_start, Size_Type r_end) const;

    /** Find overlapping mappings in a read chunk.
     * @param ce_cbptr RC to scan.
     * @param region_cont Set of ranges of the chunk's read which contain overlapping mappings.
     */
    void find_overlapping_regions(Read_Chunk_CBPtr rc_cbptr, Range_Cont& region_cont) const;

    void cat_read_contigs(Read_Entry_BPtr re_bptr);

    /** Data members */
    Unmapper _u;
    Mutation_Fact _mut_fact;
    Mutation_Chunk_Adapter_Fact _mca_fact;
    Read_Chunk_Fact _rc_fact;
    Read_Entry_Fact _re_fact;
    Contig_Entry_Fact _ce_fact;

    Read_Entry_Cont _re_cont;
    Contig_Entry_Cont _ce_cont;

    Size_Type _unmap_trigger_len;
    int _aux_coverage;
    bool _cat_at_step;
    bool _trim_tuc_step;

    BWTIndexSet _index_set;
    BWTIndexSet _aux_index_set;
    map< unsigned, string > _iid_to_sid_m;
    map< string, unsigned > _sid_to_iid_m;

}; // class Graph

} // namespace MAC


#endif
