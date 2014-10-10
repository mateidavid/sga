#ifndef __DNA_SEQUENCE_HPP
#define __DNA_SEQUENCE_HPP

#include <stdexcept>
#include <mutex>
#include <chrono>
#include <boost/iterator/transform_iterator.hpp>
#include "shortcuts.hpp"


namespace dna_sequence
{

template < typename String_Type, typename Size_Type >
class DNA_Sequence;

/** Static class holding complement table. */
class complement
{
public:
    /** Function that complements one base. */
    static char base(char c)
    {
        return table()[static_cast< unsigned >(c)];
    }
    /** Iterator modifier that complements the output of the iterator. */
    template < typename Iterator >
    static boost::transform_iterator< decltype(&complement::base), Iterator > iterator(const Iterator& it)
    {
        return boost::make_transform_iterator(it, &complement::base);
    }
    /** Function that complements a sequence. */
    template < typename Sequence_Type >
    static Sequence_Type sequence(const Sequence_Type& s)
    {
        return Sequence_Type(iterator(s.cbegin()), iterator(s.cend()));
    }
private:
    /** Complement table holder (thread-safe initialization). */
    static const char* table()
    {
        static char _complement_table[128];
        static bool _inited = false;
        if (not _inited)
        {
            static std::mutex _mutex;
            {
                std::lock_guard< std::mutex > _lock(_mutex);
                if (not _inited)
                {
                    _complement_table['A'] = 'T';
                    _complement_table['a'] = 't';
                    _complement_table['C'] = 'G';
                    _complement_table['c'] = 'g';
                    _complement_table['G'] = 'C';
                    _complement_table['g'] = 'c';
                    _complement_table['T'] = 'A';
                    _complement_table['N'] = 'N';
                    _complement_table['n'] = 'n';
                    _complement_table['X'] = 'X';
                    _complement_table['x'] = 'x';
                    _complement_table['.'] = '.';
                }
            }
            _inited = true;
        }
        return _complement_table;
    }
}; // class complement

namespace detail
{

template < typename String_Type, typename Size_Type >
class Proxy;

/** Iterator for proxy objects */
class Proxy_Iterator
    : public boost::iterator_adaptor< Proxy_Iterator, const char*, char, boost::random_access_traversal_tag, char >
{
private:
    typedef boost::iterator_adaptor< Proxy_Iterator, const char*, char, boost::random_access_traversal_tag, char > Adaptor_Base;
public:
    Proxy_Iterator()
        : Proxy_Iterator::iterator_adaptor_(0), _reversed(false), _complemented(false) {}
    explicit Proxy_Iterator(const char* other, bool reversed = false, bool complemented = false)
        : Proxy_Iterator::iterator_adaptor_(other), _reversed(reversed), _complemented(complemented) {}

    using typename Adaptor_Base::value_type;
    using typename Adaptor_Base::reference;
    using typename Adaptor_Base::difference_type;

private:
    friend class boost::iterator_core_access;

    char dereference() const
    {
        char base = not _reversed? *(this->base()) : *(this->base() - 1);
        return not _complemented? base : complement::base(base);
    }
    void increment() { increment_dir(true); }
    void decrement() { increment_dir(false); }
    void increment_dir(bool dir)
    {
        if (dir != _reversed)
        {
            ++this->base_reference();
        }
        else
        {
            --this->base_reference();
        }
    }
    void advance(difference_type n)
    {
        this->base_reference() += not _reversed? n : -n;
    }
    difference_type distance_to(const Proxy_Iterator& other) const
    {
        return not _reversed? other.base() - this->base() : this->base() - other.base();
    }

    bool _reversed;
    bool _complemented;
}; // class Proxy_Iterator

/** Proxy class.
 * Its role is to accumulate rev(), comp(), and substr() operations
 * applied on a DNA_Sequence object. One can use a proxy object in 3 ways:
 * - retrieve a given position using at() and operator[]()
 * - print object
 * - construct a new sequence object
 */
template < typename String_Type, typename Size_Type >
class Proxy
{
private:
    typedef DNA_Sequence< String_Type, Size_Type > dna_sequence_type;
    typedef Proxy_Iterator iterator;

    friend class DNA_Sequence< String_Type, Size_Type >;

    DELETE_DEF_CTOR(Proxy)
    DELETE_COPY_CTOR(Proxy)
    DEFAULT_MOVE_CTOR(Proxy)
    DELETE_COPY_ASOP(Proxy)
    DELETE_MOVE_ASOP(Proxy)

    Proxy(const dna_sequence_type* seq_p, Size_Type start, Size_Type len, bool reversed, bool complemented)
    : _seq_p(seq_p), _start(start), _len(len), _reversed(reversed), _complemented(complemented) {}

public:
    size_t size() const { return _len; }
    bool empty() const { return size() > 0; }

    /** Further rev/comp/substr applied on this proxy */
    Proxy& rev(bool reversed = true) { _reversed = (_reversed != reversed); return *this; }
    Proxy& comp(bool complemented = true) { _complemented = (_complemented != complemented); return *this; }
    Proxy& revcomp(bool revcomped = true) { return (*this).rev(revcomped).comp(revcomped); return *this; }
    Proxy& substr(size_t pos = 0, size_t len = String_Type::npos)
    {
        if (pos > size())
        {
            throw std::out_of_range("pos > size");
        }
        if (len == String_Type::npos)
        {
            len = size() - pos;
        }
        else
        {
            if (pos + len > size())
            {
                throw std::out_of_range("pos + len > size");
            }
        }
        if (not _reversed)
        {
            _start += pos;
        }
        else
        {
            _start = (_start + _len - (pos + len));
        }
        _len = len;
        return *this;
    }

    /** Element access */
    char operator [] (Size_Type pos) const
    {
        return *(begin() + pos);
    }
    char at(Size_Type pos) const
    {
        if (pos >= _len)
        {
            throw std::out_of_range("pos >= size");
        }
        return (*this)[pos];
    }

    /** Iterators */
    iterator begin() const { return iterator(&(*_seq_p)[not _reversed? _start : _start + _len], _reversed, _complemented); }
    iterator end() const { return iterator(&(*_seq_p)[not _reversed? _start + _len : _start], _reversed, _complemented); }

    /** Convert back to string */
    operator String_Type () const
    {
        if (not _reversed and not _complemented)
        {
            return String_Type(&(*_seq_p)[_start], &(*_seq_p)[_start + _len]);
        }
        else
        {
            return String_Type(begin(), end());
        }
    }

    /** Formatted output operator */
    friend std::ostream& operator << (std::ostream& os, const Proxy& rhs)
    {
        if (not rhs._reversed and not rhs._complemented)
        {
            std::copy(&(*rhs._seq_p)[rhs._start], &(*rhs._seq_p)[rhs._start + rhs._len],
                      std::ostream_iterator< char >(os));
        }
        else
        {
            std::copy(rhs.begin(), rhs.end(),
                      std::ostream_iterator< char >(os));
        }
        return os;
    }

private:
    const dna_sequence_type* const _seq_p;
    Size_Type _start;
    Size_Type _len;
    bool _reversed;
    bool _complemented;
}; // class Proxy

} // namespace detail

/** DNA_Sequence class.
 * Implemented as a non-orthodox derivation from a regular string class.
 * The derived class adds no new data members, but provides new functions
 * rev(), comp(), revcomp(), and overloads substr(). All of these produce
 * proxy objects that accumulate the corresponding restrictions.
 */
template < typename String_Type, typename Size_Type = size_t >
class DNA_Sequence
    : public String_Type
{
private:
    typedef detail::Proxy< String_Type, Size_Type > proxy_type;
public:
    DEFAULT_DEF_CTOR(DNA_Sequence)
    DEFAULT_COPY_CTOR(DNA_Sequence)
    DEFAULT_MOVE_CTOR(DNA_Sequence)
    DEFAULT_COPY_ASOP(DNA_Sequence)
    DEFAULT_MOVE_ASOP(DNA_Sequence)

    /** Construction from base */
    DNA_Sequence(const String_Type& seq) : String_Type(seq) {}
    DNA_Sequence(String_Type&& seq) : String_Type(std::move(seq)) {}
    DNA_Sequence& operator = (const String_Type& seq)
    {
        if (&seq != static_cast< String_Type* >(this))
        {
            *static_cast< String_Type* >(this) = seq;
        }
        return *this;
    }
    DNA_Sequence& operator = (String_Type&& seq)
    {
        if (&seq != static_cast< String_Type* >(this))
        {
            *static_cast< String_Type* >(this) = std::move(seq);
        }
        return *this;
    }

    /** Functions producing proxies */
    proxy_type rev(bool reversed = true) const
    {
        return proxy_type(this, 0, this->size(), reversed, false);
    }
    proxy_type comp(bool complemented = true) const
    {
        return proxy_type(this, 0, this->size(), false, complemented);
    }
    proxy_type revcomp(bool revcomped = true) const
    {
        return proxy_type(this, 0, this->size(), revcomped, revcomped);
    }
    proxy_type substr(size_t pos = 0, size_t len = String_Type::npos) const
    {
        if (pos > this->size())
        {
            throw std::out_of_range("pos > size");
        }
        if (len == String_Type::npos)
        {
            len = this->size() - pos;
        }
        else
        {
            if (pos + len > this->size())
            {
                throw std::out_of_range("pos + len > size");
            }
        }
        return proxy_type(this, pos, len, false, false);
    }
}; // class DNA_Sequence

} // namespace dna_sequence


#endif
