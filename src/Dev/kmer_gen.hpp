//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

// Generate kmers out of a given string, skipping bad bases and bases with low qv.

#ifndef __KMER_GEN_HPP
#define __KMER_GEN_HPP

#include <string>


class Kmer_gen
{
public:
    /**
    * @brief Constructor.
    *
    * @param _sq_p String to generate kmers from
    * @param _k k-mer size
    * @param _allowed_chars_p String of valid chars (optional: if NULL/empty, no check)
    * @param _qv_p Quality value string (optional: if NULL/empty, no check)
    * @param _min_qv Minimum qv (bases with smaller qv are skipped)
    * @param _phred_offset Value to subtract from qvs in _qv_p
    */
    Kmer_gen(const std::string* _sq_p, size_t _k,
             const std::string* _allowed_chars_p = NULL,
             const std::string* _qv_p = NULL, int _min_qv = 0, int _phred_offset = 33)
        : sq_p(_sq_p),
          allowed_chars_p(_allowed_chars_p),
          qv_p(_qv_p),
          k(_k),
          start(0),
          end(0),
          min_qv(_min_qv),
          phred_offset(_phred_offset) {}

    /**
    * @brief Restart kmer generation.
    */
    void reset() { start = 0; end = 0; }

    /**
    * @brief Get next kmer.
    *
    * @return Start index (0-based) of next kmer; string::npos if no kmers left.
    */
    size_t get_next();

private:
    const std::string* sq_p;
    const std::string* allowed_chars_p;
    const std::string* qv_p;
    size_t k;
    size_t start;
    size_t end;
    int min_qv;
    int phred_offset;
};


#endif
