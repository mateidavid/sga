#ifndef __FACTORY_HPP
#define __FACTORY_HPP

#include "shortcuts.hpp" // needed before other includes for is_convertible fix

#include <cstddef>
#include <iostream>
#include <limits>
#include <deque>
#include <vector>
#include <set>
#include <functional>
#include <type_traits>
#include <boost/intrusive/pointer_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/tti/tti.hpp>

#include "global_assert.hpp"
#include "nonconst_methods.hpp"
#include "ptree.hpp"
#include "logger.hpp"


namespace bounded
{

template < typename, typename >
class Pointer;
template < typename, typename >
class Reference;
template < typename, typename >
class Factory;

namespace detail
{

/** Wrapper type used by Storage objects to store free node list without overhead. */
template < typename Value, typename Index = uint32_t >
union Value_Wrapper
{
    typedef Value val_type;
    typedef Index index_type;

    val_type _val;
    index_type _idx;

    Value_Wrapper() : _idx() {}
    Value_Wrapper(const Value_Wrapper& rhs) : _idx(rhs._idx) {}
    //Value_Wrapper(Value_Wrapper&& rhs) : _idx(rhs._idx) {}
    ~Value_Wrapper() { (&(this->_idx))->~index_type(); }
};

template < typename Value, typename Index = uint32_t >
class Storage
{
public:
    typedef Value val_type;
    typedef Index index_type;
    static_assert(not std::is_const< val_type >::value, "Storage instantiated with const type");
private:
    typedef Value_Wrapper< val_type, index_type > wrapper_type;

public:
    // disable copy & move
    Storage(bool activate = true) : _cont(), _next_free_idx(0), _load(0) { if (activate) { make_active(); } }
    DELETE_COPY_CTOR(Storage)
    DELETE_MOVE_CTOR(Storage)
    DELETE_COPY_ASOP(Storage)
    DELETE_MOVE_ASOP(Storage)

    void make_active() { active_ptr() = this; }

    /** Get and set active storage. */
    static Storage*& active_ptr()
    {
        static Storage* _active_ptr = nullptr;
        return _active_ptr;
    }

    /** Static methods that use active storage. */
    static index_type allocate() { ASSERT(active_ptr()); return active_ptr()->ns_allocate(); }
    static void deallocate(index_type idx) { ASSERT(active_ptr()); active_ptr()->ns_deallocate(idx); }
    static val_type& elem_at(index_type idx) { ASSERT(active_ptr()); return active_ptr()->ns_wrapper_at(idx)._val; }
    static size_t size() { ASSERT(active_ptr()); return active_ptr()->ns_size(); }
    static size_t unused() { ASSERT(active_ptr()); return active_ptr()->ns_unused(); }
    static size_t used() { ASSERT(active_ptr()); return active_ptr()->ns_used(); }
    static std::set< index_type > unused_set() { ASSERT(active_ptr()); return active_ptr()->ns_unused_set(); }

    boost::property_tree::ptree stats() const
    {
        return ptree().put("size", ns_size()).put("unused", ns_unused());
    }

    void save(std::ostream& os) const
    {
        os.write(reinterpret_cast< const char* >(&_next_free_idx), sizeof(_next_free_idx));
        uint64_t sz = _cont.size();
        os.write(reinterpret_cast< const char* >(&sz), sizeof(sz));
        size_t bytes = 0;
        for (const auto& e : _cont)
        {
            os.write(reinterpret_cast< const char* >(&e), sizeof(e));
            bytes += sizeof(e);
        }
        LOG("io", info) << ptree("factory_save")
            .put("type", typeid(Value).name())
            .put("addr", static_cast< const void* >(this))
            .put("next_free_idx", _next_free_idx)
            .put("sz", sz)
            .put("bytes", bytes);
    }

    void load(std::istream& is)
    {
        is.read(reinterpret_cast< char* >(&_next_free_idx), sizeof(_next_free_idx));
        uint64_t sz;
        is.read(reinterpret_cast< char* >(&sz), sizeof(sz));
        _cont.resize(sz);
        size_t bytes = 0;
        for (auto& e : _cont)
        {
            is.read(reinterpret_cast< char* >(&e), sizeof(e));
            bytes += sizeof(e);
        }
        LOG("io", info) << ptree("factory_load")
            .put("type", typeid(Value).name())
            .put("addr", static_cast< const void* >(this))
            .put("next_free_idx", _next_free_idx)
            .put("sz", sz)
            .put("bytes", bytes);
    }

private:
    /** Dereference wrapper. */
    wrapper_type& ns_wrapper_at(index_type idx) const
    {
        return const_cast< wrapper_type& >(
#ifndef NDEBUG
            _cont.at(idx)
#else
            _cont[idx]
#endif
            );
    }

    /** Allocate space for new element. */
    index_type ns_allocate()
    {
        ASSERT(_next_free_idx <= _cont.size());
        index_type idx = _next_free_idx;
        if (_next_free_idx < _cont.size())
        {
            _next_free_idx = ns_wrapper_at(_next_free_idx)._idx;
        }
        else
        {
            _cont.push_back(wrapper_type());
            ++_next_free_idx;
        }
        ++_load;
        LOG("factory", debug1) << ptree("allocate")
            .put("type", typeid(Value).name())
            .put("addr", static_cast< const void* >(this))
            .put("idx", idx);
        return idx;
    }

    /** Deallocate element. */
    void ns_deallocate(index_type idx)
    {
        LOG("factory", debug1) << ptree("deallocate")
            .put("type", typeid(Value).name())
            .put("addr", static_cast< const void* >(this))
            .put("idx", idx);
        ns_wrapper_at(idx)._idx = _next_free_idx;
        _next_free_idx = idx;
        --_load;
    }

    /** Get currently allocated size. */
    size_t ns_size() const { return _cont.size(); }
    /** Get number of unused entries. */
    size_t ns_unused() const { return _cont.size() - _load; }
    /** Get number of unused entries. */
    size_t ns_used() const { return _load; }
    /** Get set of unused entries. */
    std::set< index_type > ns_unused_set() const
    {
        std::set< index_type > res;
        for (auto crt = _next_free_idx; crt < _cont.size(); crt = ns_wrapper_at(crt)._idx)
        {
            res.insert(crt);
        }
        return res;
    }

    std::deque< wrapper_type > _cont;
    //std::vector< wrapper_type > _cont;
    index_type _next_free_idx;
    size_t _load;
}; // class Storage

template < typename Value, typename Index = uint32_t >
struct Identifier
{
    typedef Value val_type;
    typedef Index index_type;
    typedef typename std::is_void< val_type > is_void_t; // true_type iff T is (cv-) void
    typedef typename std::is_const< val_type > is_const_t; // true_type iff T is const-qualified
    static_assert(not is_const_t::value, "Identifier instantiated with const type");
    typedef typename std::add_lvalue_reference< val_type >::type raw_ref_type; // void if T==void
    typedef typename boost::mpl::if_< is_void_t, void, Storage< val_type, index_type > >::type storage_type; // void if T==void

    static index_type const null_val = std::numeric_limits< index_type >::max();

    Identifier(index_type idx = null_val) : _idx(idx) {}

    // allow copy & move
    DEFAULT_COPY_CTOR(Identifier)
    DEFAULT_MOVE_CTOR(Identifier)
    DEFAULT_COPY_ASOP(Identifier)
    DEFAULT_MOVE_ASOP(Identifier)

    /** Dereferencing operator, non-void types only. */
    template < T_ENABLE_IF((not is_void_t::value)) >
    raw_ref_type dereference() const
    {
        ASSERT(*this);
        return storage_type::elem_at(_idx);
    }

    /** Bool conversion: via pointer to private member fcn. */
    BOOL_CONVERSION(Identifier, public, (_idx != null_val))

    /** Increment & decrement operators. */
    Identifier& operator ++ ()    { ++_idx; return *this; }
    Identifier  operator ++ (int) { Identifier res(*this); ++(*this); return res; }
    Identifier& operator -- ()    { --_idx; return *this; }
    Identifier  operator -- (int) { Identifier res(*this); --(*this); return res; }

    /** Comparison operators. */
    friend bool operator == (const Identifier& lhs, const Identifier& rhs) { return lhs._idx == rhs._idx; }
    friend bool operator != (const Identifier& lhs, const Identifier& rhs) { return lhs._idx != rhs._idx; }
    friend bool operator <  (const Identifier& lhs, const Identifier& rhs) { return lhs._idx <  rhs._idx; }
    friend bool operator <= (const Identifier& lhs, const Identifier& rhs) { return lhs._idx <= rhs._idx; }
    friend bool operator >  (const Identifier& lhs, const Identifier& rhs) { return lhs._idx >  rhs._idx; }
    friend bool operator >= (const Identifier& lhs, const Identifier& rhs) { return lhs._idx >= rhs._idx; }

    //friend std::ostream& operator << (std::ostream& os, const Identifier& rhs) { os << rhs._idx; return os; }
    index_type to_int() const { return _idx; }
    boost::property_tree::ptree to_ptree() const
    {
        boost::property_tree::ptree pt;
        if (*this)
        {
            pt.put_value(_idx);
        }
        else
        {
            pt.put_value("null");
        }
        return pt;
    }

    index_type _idx;
}; // struct Identifier

} // namespace detail

template < typename Value, typename Index = uint32_t >
class Pointer
{
public:
    typedef Value val_type;
    typedef Index index_type;
    typedef val_type* raw_ptr_type;
    typedef typename std::remove_const< val_type >::type unqual_val_type;
private:
    typedef typename std::is_void< val_type > is_void_t; // true_type iff T is (cv-) void
    typedef typename std::is_const< val_type > is_const_t; // true_type iff T is const-qualified
    typedef detail::Identifier< unqual_val_type, index_type > idn_type;
public:
    typedef typename boost::mpl::if_< is_void_t, void, Reference< val_type, index_type > >::type reference;

    Pointer(std::nullptr_t = nullptr) {}

    // allow copy & move
    DEFAULT_COPY_CTOR(Pointer)
    DEFAULT_MOVE_CTOR(Pointer)
    DEFAULT_COPY_ASOP(Pointer)
    DEFAULT_MOVE_ASOP(Pointer)

    //CONST_CONVERSIONS(Pointer, val_type)
    // construct const Pointer from non-const Pointer
    template < typename Other_Value,
               T_ENABLE_IF((std::is_same< Value, const unqual_val_type >::value
                            and std::is_same< Other_Value, unqual_val_type >::value)) >
    Pointer(const Pointer< Other_Value, Index >& other) : _idn(other._idn) {}

    /** Remove constness of referred type */
    const Pointer< unqual_val_type, Index >& unconst() const
    { return *reinterpret_cast< const Pointer< unqual_val_type, Index >* >(this); }
    NONCONST_METHOD((Pointer< unqual_val_type, Index >&), unconst)

    /** Constructor from index. */
    static Pointer from_index(const index_type& idx) { Pointer p; p._idn._idx = idx; return p; }

    /** Raw pointer conversion. */
    //raw_ptr_type raw() const { return *this? &_idn.dereference() : nullptr; }
    raw_ptr_type raw() const { return &_idn.dereference(); }

    /** Dereferencing operators. */
    reference operator * () const { ASSERT(*this); return reference(_idn); }
    raw_ptr_type operator -> () const { return raw(); }

    /** Pointer traits interface. */
    static Pointer pointer_to(reference r) { return &r; }
    // workaround for: https://svn.boost.org/trac/boost/ticket/10853
    template < typename Other_Value >
    static Pointer const_cast_from(const Pointer< Other_Value, Index >& other) { return other.unconst(); }

    /** Bool conversion: via pointer to private member fcn. */
    BOOL_CONVERSION(Pointer, public, (_idn))

    /** Increment & decrement operators. */
    Pointer& operator ++ ()    { ++(this->_idn); return *this; }
    Pointer  operator ++ (int) { Pointer res(*this); ++(*this); return res; }
    Pointer& operator -- ()    { --(this->_idn); return *this; }
    Pointer  operator -- (int) { Pointer res(*this); --(*this); return res; }

    /** Comparison operators. */
    friend bool operator == (const Pointer& lhs, const Pointer& rhs) { return lhs._idn == rhs._idn; }
    friend bool operator != (const Pointer& lhs, const Pointer& rhs) { return lhs._idn != rhs._idn; }
    friend bool operator <  (const Pointer& lhs, const Pointer& rhs) { return lhs._idn <  rhs._idn; }
    friend bool operator <= (const Pointer& lhs, const Pointer& rhs) { return lhs._idn <= rhs._idn; }
    friend bool operator >  (const Pointer& lhs, const Pointer& rhs) { return lhs._idn >  rhs._idn; }
    friend bool operator >= (const Pointer& lhs, const Pointer& rhs) { return lhs._idn >= rhs._idn; }

    index_type to_int() const { return _idn.to_int(); }
    boost::property_tree::ptree to_ptree() const { return _idn.to_ptree(); }

private:
    friend class Factory< unqual_val_type, index_type >;
    friend class Pointer< const unqual_val_type, index_type >;

    idn_type _idn;
}; // class Pointer

namespace detail
{

template < typename T, bool = false >
struct to_ptree_caller_impl
{
    static void call(boost::property_tree::ptree&, const T&) {}
};
template < typename T >
struct to_ptree_caller_impl< T, true >
{
    static void call(boost::property_tree::ptree& pt, const T& val) { pt.put_child("deref", val.to_ptree()); }
};
BOOST_TTI_HAS_MEMBER_FUNCTION(to_ptree)
template < typename T >
struct to_ptree_caller
    : public to_ptree_caller_impl< T, has_member_function_to_ptree< const T, boost::property_tree::ptree >::value > {};

} // namespace detail

/** Reference. */
template < typename Value, typename Index = uint32_t >
class Reference
{
public:
    typedef Value val_type;
    typedef Index index_type;
private:
    typedef typename std::is_void< val_type > is_void_t; // true_type iff T is (cv-) void
    typedef typename std::is_const< val_type > is_const_t; // true_type iff T is const-qualified
    typedef detail::to_ptree_caller< val_type > to_ptree_caller_t;
    static_assert(not is_void_t::value, "Reference instantiated with void type");
public:
    typedef val_type& raw_ref_type;
    typedef typename std::remove_const< val_type >::type unqual_val_type;
    typedef Pointer< val_type, index_type > ptr_type;
private:
    typedef detail::Identifier< unqual_val_type, index_type > idn_type;

public:
    // allow construction only
    DELETE_DEF_CTOR(Reference)
    DEFAULT_COPY_CTOR(Reference)
    DEFAULT_MOVE_CTOR(Reference)
    DELETE_MOVE_ASOP(Reference)

    //CONST_CONVERSIONS(Reference, val_type)
    // construct const Reference from non-const Reference
    template < typename Other_Value,
               T_ENABLE_IF((std::is_same< Value, const unqual_val_type >::value
                            and std::is_same< Other_Value, unqual_val_type >::value)) >
    Reference(const Reference< Other_Value, Index >& other) : _idn(other._idn) {}

    /** Remove constness of referred type */
    const Reference< unqual_val_type, Index >& unconst() const
    { return *reinterpret_cast< const Reference< unqual_val_type, Index >* >(this); }
    NONCONST_METHOD((Reference< unqual_val_type, Index >&), unconst)

    // get raw reference
    raw_ref_type raw() const { ASSERT(_idn); return _idn.dereference(); }

    // address-of operator returns bounded pointer
    ptr_type operator & () const { return ptr_type::from_index(_idn._idx); }

    // implicit conversion to raw reference
    operator raw_ref_type () const { return raw(); }

    // allow referred element assignment (non-const only)
    template < T_ENABLE_IF((not is_const_t::value)) >
    const Reference& operator = (const val_type& rhs) const { raw() = rhs; return *this; }
    template < T_ENABLE_IF((not is_const_t::value)) >
    const Reference& operator = (const Reference& rhs) const { raw() = rhs.raw(); return *this; }

    boost::property_tree::ptree to_ptree() const
    {
        ptree pt;
        pt.put("baddr", _idn._idx);
        if (_idn)
        {
            to_ptree_caller_t::call(pt, raw());
        }
        return pt;
    }

private:
    friend class Pointer< val_type, index_type >;
    friend class Reference< const unqual_val_type, index_type >;

    Reference(const idn_type& idn) : _idn(idn) {}

    const idn_type _idn;
}; // class Reference

namespace detail {

template < typename Value, typename Index = uint32_t >
struct Cloner
{
    typedef Value val_type;
    typedef Index index_type;
    static_assert(not std::is_const< val_type >::value, "Cloner instantiated with const type");
    static_assert(std::is_copy_constructible< val_type >::value, "Cloner instantiated with non-copy-constructible type");
    typedef Pointer< val_type, index_type > ptr_type;
    typedef Factory< val_type, index_type > fact_type;

    ptr_type operator () (const val_type& t)
    {
        return fact_type::new_elem(t);
    }
};

template < typename Value, typename Index = uint32_t >
struct Disposer
{
    typedef Value val_type;
    typedef Index index_type;
    static_assert(not std::is_const< val_type >::value, "Disposer instantiated with const type");
    typedef Pointer< val_type, index_type > ptr_type;
    typedef Factory< val_type, index_type > fact_type;

    void operator () (ptr_type ptr)
    {
        fact_type::del_elem(ptr);
    }
};

} // namespace detail

/** Factory class that manages storage for objects of type T. */
template < typename Value, typename Index = uint32_t >
class Factory
    : public detail::Storage< Value, Index >
{
private:
    typedef detail::Storage< Value, Index > Base;
public:
    typedef Value val_type;
    typedef Index index_type;
    static_assert(not std::is_const< val_type >::value, "Factory instantiated with const type");
    typedef Pointer< val_type, index_type > ptr_type;
    typedef Pointer< const val_type, index_type > const_ptr_type;
    typedef Reference< val_type, index_type > ref_type;
    typedef Reference< const val_type, index_type > const_ref_type;
    typedef detail::Cloner< val_type, index_type > cloner_type;
    typedef detail::Disposer< val_type, index_type > disposer_type;
private:
    typedef detail::Identifier< val_type, index_type > idn_type;
    typedef detail::Value_Wrapper< val_type, index_type > wrapper_type;

public:
    /** Default constructor.
     * @param activate Bool; if true, the object is used to resolve all managed pointers.
     */
    Factory(bool activate = true) : Base(activate) {}

    // disable copy & move
    DELETE_COPY_CTOR(Factory)
    DELETE_MOVE_CTOR(Factory)
    DELETE_COPY_ASOP(Factory)
    DELETE_MOVE_ASOP(Factory)

    void make_active() { Base::make_active(); }

    /** Allocate and construct new element. */
    template < typename ...Args >
    static ptr_type new_elem(Args&& ...args)
    {
        ptr_type res = ptr_type::from_index(Base::allocate()); // use active storage
        new (res.raw()) val_type(std::forward<Args>(args)...);
        return res;
    }

    /** Destruct and deallocate given element. */
    static void del_elem(const_ptr_type elem_cptr)
    {
        elem_cptr->~val_type();
        Base::deallocate(elem_cptr._idn._idx); // use active storage
    }
}; // class Factory

/** Stateless allocator that uses the active storage. */
template < typename Value, typename Index = uint32_t >
class Allocator
{
public:
    typedef Value value_type;
    typedef Index index_type;
    typedef Pointer< value_type, index_type > pointer;
    typedef Pointer< const value_type, index_type > const_pointer;
    typedef Reference< value_type, index_type > reference;
    typedef Reference< const value_type, index_type > const_reference;
private:
    typedef detail::Storage< value_type, index_type > storage_type;

public:
    DEFAULT_DEF_CTOR(Allocator)
    DEFAULT_COPY_CTOR(Allocator)
    DEFAULT_MOVE_CTOR(Allocator)
    DEFAULT_COPY_ASOP(Allocator)
    DEFAULT_MOVE_ASOP(Allocator)

    pointer allocate(size_t n, std::allocator< void >::const_pointer = nullptr)
    {
        static_cast< void >(n);
        ASSERT(n == 1);
        return pointer::from_index(storage_type::allocate()); // use active storage
    }
    void deallocate(pointer p, size_t n)
    {
        static_cast< void >(n);
        ASSERT(n == 1);
        storage_type::deallocate(p); // use active storage
    }
}; // class Allocator

template < typename Value, typename Index = uint32_t >
class Pointer_Holder
{
public:
    typedef Value val_type;
    typedef Index index_type;
private:
    typedef typename std::remove_const< val_type >::type unqual_val_type;
    typedef Factory< unqual_val_type, index_type > fact_type;
public:
    typedef typename fact_type::ptr_type ptr_type;
    typedef typename fact_type::const_ptr_type const_ptr_type;

public:
    template <typename... Args>
    Pointer_Holder(Args&& ... args) : _ptr(fact_type::new_elem(std::forward<Args>(args)...)) {}

    DEFAULT_COPY_CTOR(Pointer_Holder)
    DEFAULT_MOVE_CTOR(Pointer_Holder)
    DEFAULT_COPY_ASOP(Pointer_Holder)
    DEFAULT_MOVE_ASOP(Pointer_Holder)

    ~Pointer_Holder() { fact_type::del_elem(_ptr); }

    const_ptr_type get_node () const { return _ptr; }
    ptr_type get_node () { return _ptr; }

private:
    ptr_type _ptr;
}; // class Pointer_Holder

} // namespace bounded

/// Template specialization to allow bounded pointers in unordered_set containers.
///
namespace std
{
    template < typename Value, typename Index >
    struct hash< bounded::Pointer< Value, Index > >
    {
        size_t operator () (const bounded::Pointer< Value, Index >& p) const
        {
            return hash< Index >()(p.to_int());
        }
    }; // struct hash
}


#endif
