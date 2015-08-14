//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include <boost/range/adaptor/map.hpp>
#include "Mutation.hpp"
#include "Contig_Entry.hpp"


namespace MAC
{

unsigned Mutation::find_or_add_alt(const Seq_Proxy_Type& seq)
{
    return 1 + _alt_cont.find_or_add(seq);
}

unsigned Mutation::remove_alt(unsigned i)
{
    ASSERT(i > 0);
    --i;
    if (i < _alt_cont.size() - 1)
    {
        _alt_cont.swap_alleles(i, _alt_cont.size() - 1);
    }
    _alt_cont.pop_back();
    return n_alleles();
}

Seq_Proxy_Type Mutation::allele_seq(unsigned i) const
{
    if (i == 0)
    {
        return ce_cbptr()->substr(rf_start(), rf_len());
    }
    else
    {
        return _alt_cont.at(i - 1).seq();
    }
}

void Mutation::merge_right(Mutation_BPtr rhs_mut_bptr, map< pair< unsigned, unsigned >, unsigned >& allele_map)
{
    ASSERT(is_unlinked());
    ASSERT(rhs_mut_bptr->is_unlinked());
    ASSERT(ce_cbptr() == rhs_mut_bptr->ce_cbptr());
    if (rf_end() != rhs_mut_bptr->rf_start())
    {
        abort();
    }
    Allele_Cont new_cont;
    for (const auto& p : allele_map | ba::map_keys)
    {
        if (p.first == 0 and p.second == 0)
        {
            allele_map[p] = 0;
        }
        else
        {
            new_cont.push_back(Seq_Type(allele_seq(p.first)) + Seq_Type(rhs_mut_bptr->allele_seq(p.second)));
            allele_map[p] = new_cont.size();
        }
    }
    _len += rhs_mut_bptr->rf_len();
    _alt_cont = move(new_cont);
}

pair< Mutation_BPtr, map< unsigned, pair< unsigned, unsigned > > >
Mutation::cut(Size_Type offset)
{
    ASSERT(is_unlinked());
    Size_Type rf_offset = min(offset, rf_len());
    Mutation_BPtr rhs_mut_bptr = Mutation_Fact::new_elem(rf_start() + rf_offset, rf_len() - rf_offset, ce_cbptr());
    Allele_Cont new_cont;
    map< Seq_Type, unsigned > lhs_alleles;
    map< Seq_Type, unsigned > rhs_alleles;
    map< unsigned, pair< unsigned, unsigned > > allele_map;
    lhs_alleles[ce_cbptr()->substr(rf_start(), rf_offset)] = 0;
    rhs_alleles[ce_cbptr()->substr(rf_start() + rf_offset, rf_len() - rf_offset)] = 0;
    allele_map[0] = make_pair(0, 0);
    unsigned i = 1;
    for (const auto a_cbptr : _alt_cont | referenced)
    {
        Size_Type al_offset = min(offset, a_cbptr->seq().size());
        Seq_Type lhs_seq = a_cbptr->seq().substr(0, al_offset);
        Seq_Type rhs_seq = a_cbptr->seq().substr(al_offset);
        if (lhs_alleles.count(lhs_seq) == 0)
        {
            new_cont.push_back(Seq_Type(lhs_seq));
            lhs_alleles[lhs_seq] = new_cont.size();
        }
        if (rhs_alleles.count(rhs_seq) == 0)
        {
            rhs_mut_bptr->_alt_cont.push_back(Seq_Type(rhs_seq));
            rhs_alleles[rhs_seq] = rhs_mut_bptr->_alt_cont.size();
        }
        allele_map[i] = make_pair(lhs_alleles[lhs_seq], rhs_alleles[rhs_seq]);
        ++i;
    }
    _len = rf_offset;
    _alt_cont = move(new_cont);
    return make_pair(rhs_mut_bptr, allele_map);
}

void Mutation::simplify()
{
    ASSERT(is_unlinked());
    // compute common suffix len
    {
        Seq_Proxy_Type rf_seq = allele_seq(0);
        Size_Type common_suffix_len = rf_seq.size();
        for (auto al_bptr : _alt_cont | referenced)
        {
            if (common_suffix_len == 0) break;
            unsigned i = 0;
            while (i < min(common_suffix_len, al_bptr->seq().size())
                   and al_bptr->seq()[al_bptr->seq().size() - 1 - i] == rf_seq[rf_seq.size() - 1 - i]) ++i;
            if (i < common_suffix_len) common_suffix_len = i;
        }
        if (common_suffix_len > 0)
        {
            _len -= common_suffix_len;
            for (auto al_bptr : _alt_cont | referenced)
            {
                al_bptr->seq().resize(al_bptr->seq().size() - common_suffix_len);
            }
        }
    }
    // compute common prefix len
    {
        Seq_Proxy_Type rf_seq = allele_seq(0);
        Size_Type common_prefix_len = rf_seq.size();
        for (auto al_bptr : _alt_cont | referenced)
        {
            if (common_prefix_len == 0) break;
            unsigned i = 0;
            while (i < min(common_prefix_len, al_bptr->seq().size())
                   and al_bptr->seq()[i] == rf_seq[i]) ++i;
            if (i < common_prefix_len) common_prefix_len = i;
        }
        if (common_prefix_len > 0)
        {
            _start += common_prefix_len;
            _len -= common_prefix_len;
            for (auto al_bptr : _alt_cont | referenced)
            {
                al_bptr->seq() = al_bptr->seq().substr(common_prefix_len);
            }
        }
    }
}

void Mutation::reverse()
{
    _start = ce_cbptr()->len() - rf_end();
    for (auto al_bptr : _alt_cont | referenced)
    {
        al_bptr->seq() = al_bptr->seq().revcomp();
    }
}

/*
void Mutation::to_stream(ostream& os, Mutation_CBPtr mut_cbptr, Contig_Entry_CBPtr ce_cbptr)
{
    ASSERT(mut_cbptr);
    ASSERT(ce_cbptr);
    set< Read_Chunk_CBPtr > qr_set;
    set< Read_Chunk_CBPtr > full_rf_set;
    tie(qr_set, full_rf_set, ignore) = ce_cbptr->mut_support(mut_cbptr);
    os << "mut " << setw(5) << left << mut_cbptr.to_int()
       << setw(5) << right << mut_cbptr->rf_start()
       << setw(5) << right << mut_cbptr->rf_len()
       << setw(5) << right << mut_cbptr->seq_len()
       << setw(5) << right << full_rf_set.size()
       << setw(5) << right << qr_set.size();
}
*/

void Mutation::save_strings(ostream& os, size_t& n_strings, size_t& n_bytes) const
{
    _alt_cont.save_strings(os, n_strings, n_bytes);
}

void Mutation::init_strings()
{
    _alt_cont.init_strings();
}

void Mutation::load_strings(istream& is, size_t& n_strings, size_t& n_bytes)
{
    _alt_cont.load_strings(is, n_strings, n_bytes);
}

} // namespace MAC
