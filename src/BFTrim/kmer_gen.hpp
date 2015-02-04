//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

// Generate fix-sized kmers out of valid positions in a string.
//
// Example use:
//   string s = "abracadara";
//   size_t k = 3;
//   auto kg = kmer_gen(s.size(), k, [&] (size_t i) { return s[i] != 'r'; });
//   for (auto it = kg.begin(); it != kg.end(); ++it)
//   {
//       cout << *it << ": " << s.substr(*it, k) << endl;
//   }
// Output:
//   3: aca
//   4: cad
//   5: ada

#ifndef __KMER_GEN_HPP
#define __KMER_GEN_HPP

#include <string>
#include <functional>
#include <cassert>
#include <boost/iterator/iterator_facade.hpp>


class kmer_gen;

class kmer_gen_iterator
    : public boost::iterator_facade< kmer_gen_iterator,
                                     size_t,
                                     boost::forward_traversal_tag,
                                     size_t
                                     >
{
private:
    typedef boost::iterator_facade< kmer_gen_iterator,
                                    size_t,
                                    boost::forward_traversal_tag,
                                    size_t
                                    > Base;
public:
    kmer_gen_iterator()
        : _validator(), _len(0), _k(0), _start(0), _end(0) {}

    bool is_valid() const
    {
        assert(_len > 0 and _k > 0);
        return _end == _start + _k;
    }

private:
    friend class boost::iterator_core_access;
    friend class kmer_gen;

    kmer_gen_iterator(size_t len, size_t k, std::function< bool(size_t) > validator,
                      size_t start, size_t end)
        : _validator(validator), _len(len), _k(k), _start(start), _end(end) {}

    size_t dereference() const
    {
        assert(is_valid());
        return _start;
    }

    bool equal(kmer_gen_iterator other) const
    {
        return _start == other._start and _end == other._end;
    }

    /** Get start of next valid kmer.
     * When calling this function, the iterator must be in one of 3 possible states:
     * 1. non-initialized, with _end == _start == 0
     * 2. valid, with _end == _start + _k and _start != string::npos
     * 3. past end, with _end == _start == string::npos
     */
    void increment()
    {
        assert(_len > 0 and _k > 0);
        assert((_end == _start and _start == 0)
               or (_end == _start + _k and _start != std::string::npos)
               or (_end == _start and _start == std::string::npos));

        if (_start == std::string::npos)
        {
            // already past end
            return;
        }

        if (_end == _start + _k)
        {
            // discard previously generated kmer
            ++_start;
        }

        // include next base
        ++_end;
        while (_end <= _len)
        {
            // consider string [_start, ..., _end-1]
            // check if last position is valid
            if (_validator and not _validator(_end - 1))
            {
                // invalid base: skip all kmers containing this pos
                _start = _end;
                _end = _start + 1;
                continue;
            }
            assert(_start < _end);
            assert(_end <= _start + _k);
            if (_end == _start + _k)
            {
                return;
            }
            ++_end;
        }
        _start = std::string::npos;
        _end = std::string::npos;
    }

    std::function< bool(size_t) > _validator;
    size_t _len;
    size_t _k;
    size_t _start;
    size_t _end;
}; // class kmer_gen_iterator


class kmer_gen
{
public:
    kmer_gen(size_t len, size_t k, std::function< bool(size_t) > validator = nullptr)
        : _validator(validator), _len(len), _k(k) {}

    kmer_gen_iterator begin() const { return ++kmer_gen_iterator(_len, _k, _validator, 0, 0); }
    kmer_gen_iterator end() const { return kmer_gen_iterator(_len, _k, _validator, std::string::npos, std::string::npos); }

private:
    std::function< bool(size_t) > _validator;
    size_t _len;
    size_t _k;
}; // class kmer_gen


#endif
