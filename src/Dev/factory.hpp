#ifndef __FACTORY_HPP
#define __FACTORY_HPP

#include "shortcuts.hpp" // needed before other includes for is_convertible fix

#include <cstddef>
#include <iostream>
#include <limits>
#include <deque>
#include <type_traits>
#include <boost/intrusive/pointer_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/tti/tti.hpp>
//#include <typeinfo>

#include "global_assert.hpp"
#include "nonconst_methods.hpp"
#include "ptree_tools.hpp"

#define CONST_CONVERSIONS(_type, _t) \
    operator const _type< typename boost::add_const< _t >::type >& () const \
    { return *reinterpret_cast< const _type< typename boost::add_const< _t >::type >* >(this); } \
    operator _type< typename boost::add_const< _t >::type >& () \
    { return *reinterpret_cast< _type< typename boost::add_const< _t >::type >* >(this); } \
    \
    const _type< typename boost::remove_const< _t >::type >& unconst() const \
    { return *reinterpret_cast< const _type< typename boost::remove_const< _t >::type >* >(this); } \
    _type< typename boost::remove_const< _t >::type >& unconst() \
    { return *reinterpret_cast< _type< typename boost::remove_const< _t >::type >* >(this); }

#define T_ENABLE_IF(_cond) \
    bool _unused = true, typename boost::enable_if_c< _unused and (_cond), int >::type = 42

namespace bounded
{
namespace detail
{

template < typename, typename >
struct Identifier;
template < typename, typename >
class Pointer;
template < typename, typename >
class Reference;
template < typename, typename >
struct Cloner;
template < typename, typename >
struct Disposer;
template < typename, typename >
class Storage;
template < typename, typename >
class Factory;
template < typename, typename >
class Static_Allocator;

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

    Identifier() : _idx(null_val) {}

    // allow copy & move
    DEFAULT_COPY_CTOR(Identifier)
    DEFAULT_MOVE_CTOR(Identifier)
    DEFAULT_COPY_ASOP(Identifier)
    DEFAULT_MOVE_ASOP(Identifier)

    //operator Index () const { return _ptr; }
    bool is_null() const { return _idx == null_val; }

    /** Dereferencing operator (non-void types only). */
    template < T_ENABLE_IF((not is_void_t::value)) >
    raw_ref_type dereference() const
    {
        return storage_type::elem_at(*this);
    }

    bool operator == (const Identifier& rhs) const { return _idx == rhs._idx; }
    bool operator < (const Identifier& rhs) const { return _idx < rhs._idx; }
    friend std::ostream& operator << (std::ostream& os, const Identifier& rhs)
    {
        os << "_idx=" << rhs._idx;
        return os;
    }

    index_type _idx;
};

template < typename Value, typename Index >
bool operator != (const Identifier< Value, Index >& lhs, const Identifier< Value, Index >& rhs)
{
    return !(lhs == rhs);
}

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
    typedef Identifier< unqual_val_type, index_type > idn_type;
public:
    typedef typename boost::mpl::if_< is_void_t, void, Reference< val_type, index_type > >::type ref_type;

    template < typename Other_Value >
    struct has_same_base
        : std::is_same< typename std::remove_const< Other_Value >::type,
                        typename std::remove_const< val_type >::type >
    {};

    /** Rebinder.
     * Wrap pointers to T and const T. Anything else gets transformed into a raw pointer.
     */
    template < typename Other_Value >
    struct rebind
    {
        typedef typename boost::mpl::if_< has_same_base< Other_Value >,
                                          Pointer< Other_Value, index_type >,
                                          Other_Value*
                                        >::type type;
    };

    Pointer(std::nullptr_t = nullptr) {}

    // allow copy & move
    DEFAULT_COPY_CTOR(Pointer)
    DEFAULT_MOVE_CTOR(Pointer)
    DEFAULT_COPY_ASOP(Pointer)
    DEFAULT_MOVE_ASOP(Pointer)

    CONST_CONVERSIONS(Pointer, val_type)

    // get raw pointer
    raw_ptr_type raw() const { return _idn.is_null()? nullptr : &_idn.dereference(); }

    // bool conversion: via raw ptr
    operator bool () const { return raw(); }

    //size_t to_int() const { return static_cast< size_t >(_id._ptr); }

    // increment operators
    Pointer& operator ++ () { ++(this->_idn)._idx; return *this; }
    Pointer operator ++ (int) { Pointer res(*this); ++(*this); return res; }

    // dereferencing operators
    ref_type operator * () const { return ref_type(_idn); }
    raw_ptr_type operator -> () const { return raw(); }

    boost::property_tree::ptree to_ptree() const
    {
        std::ostringstream tmp;
        tmp << _idn._idx;
        return boost::property_tree::ptree(tmp.str());
    }

    template < typename Other_Value,
               T_ENABLE_IF((has_same_base< Other_Value >::value)) >
    bool operator == (const Pointer< Other_Value, index_type >& rhs) const { return _idn == rhs._idn; }
    template < typename Other_Value,
               T_ENABLE_IF((has_same_base< Other_Value >::value)) >
    bool operator < (const Pointer< Other_Value, index_type >& rhs) const { return _idn < rhs._idn; }

private:
    template < typename, typename > friend class Pointer;
    friend class Reference< val_type, index_type >;
    //friend class Factory< unqual_val_type, index_type >;
    friend class Storage< unqual_val_type, index_type >;

    Pointer(const idn_type& idn) : _idn(idn) {}

    idn_type _idn;
};

template < typename Value, typename Index >
bool operator == (const Pointer< Value, Index >& lhs, std::nullptr_t)
{
    return lhs == Pointer< Value, Index >(nullptr);
}
template < typename Value, typename Index >
bool operator == (std::nullptr_t, const Pointer< Value, Index >& rhs)
{
    return Pointer< Value, Index >(nullptr) == rhs;
}
template < typename Value_LHS, typename Value_RHS, typename Index >
bool operator != (const Pointer< Value_LHS, Index >& lhs, const Pointer< Value_RHS, Index >& rhs)
{
    return !(lhs == rhs);
}
template < typename Value, typename Index >
bool operator != (const Pointer< Value, Index >& lhs, std::nullptr_t)
{
    return !(lhs == nullptr);
}
template < typename Value, typename Index >
bool operator != (std::nullptr_t, const Pointer< Value, Index >& rhs)
{
    return !(nullptr == rhs);
}

template < typename T, bool = false >
struct to_ptree_caller_impl
{
    static void call(ptree&, const T&) {}
};
template < typename T >
struct to_ptree_caller_impl< T, true >
{
    static void call(ptree& pt, const T& val)
    {
        pt.put("deref", val.to_ptree());
    }
};
BOOST_TTI_HAS_MEMBER_FUNCTION(to_ptree)
template < typename T >
struct to_ptree_caller
    : public to_ptree_caller_impl< T, has_member_function_to_ptree< const T, boost::property_tree::ptree >::value > {};

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
    typedef to_ptree_caller< val_type > to_ptree_caller_t;
    static_assert(not is_void_t::value, "Reference instantiated with void type");
public:
    typedef val_type& raw_ref_type;
    typedef typename std::remove_const< val_type >::type unqual_val_type;
    typedef Pointer< val_type, index_type > ptr_type;
private:
    typedef Identifier< unqual_val_type, index_type > idn_type;

public:
    // allow construction only
    DELETE_DEF_CTOR(Reference)
    DEFAULT_COPY_CTOR(Reference)
    DEFAULT_MOVE_CTOR(Reference)
    DELETE_MOVE_ASOP(Reference)

    CONST_CONVERSIONS(Reference, val_type)

    // get raw reference
    raw_ref_type raw() const { ASSERT(not _idn.is_null()); return _idn.dereference(); }

    // address-of operator returns bounded pointer
    ptr_type operator & () const { return ptr_type(_idn); }

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
        if (not _idn.is_null())
        {
            to_ptree_caller_t::call(pt, raw());
        }
        return pt;
    }

private:
    friend class Pointer< val_type, index_type >;

    Reference(const idn_type& idn) : _idn(idn) {}

    const idn_type _idn;
}; // class Reference

template < typename Value, typename Index = uint32_t >
struct Cloner
{
    typedef Value val_type;
    typedef Index index_type;
    static_assert(not std::is_const< val_type >::value, "Cloner instantiated with const type");
    static_assert(std::is_copy_constructible< val_type >::value, "Cloner instantiated with non-copy-construtible type");
    typedef Identifier< val_type, index_type > idn_type;
    typedef Pointer< val_type, index_type > ptr_type;
    typedef typename idn_type::fact_type fact_type;

    ptr_type operator () (const val_type& t)
    {
        ASSERT(fact_type::get_active_ptr());
        return fact_type::get_active_ptr()->new_elem(t);
    }
};

template < typename Value, typename Index = uint32_t >
struct Disposer
{
    typedef Value val_type;
    typedef Index index_type;
    static_assert(not std::is_const< val_type >::value, "Disposer instantiated with const type");
    typedef Identifier< val_type, index_type > idn_type;
    typedef Pointer< val_type, index_type > ptr_type;
    typedef typename idn_type::fact_type fact_type;

    void operator () (ptr_type ptr)
    {
        ASSERT(fact_type::get_active_ptr());
        fact_type::get_active_ptr()->del_elem(ptr);
    }
};

/** Wrapper type used by Factory objects to store free node list without overhead. */
template < typename Value, typename Index = uint32_t >
union Value_Wrapper
{
    typedef Value val_type;
    typedef Index index_type;
    typedef Identifier< val_type, index_type > idn_type;

    val_type _key;
    idn_type _next_free_idn;

    Value_Wrapper() : _next_free_idn() {}
    Value_Wrapper(const Value_Wrapper& rhs) : _next_free_idn(rhs._next_free_idn) {}
    ~Value_Wrapper() { (&(this->_next_free_idn))->~idn_type(); }
};

template < typename Value, typename Index = uint32_t >
class Storage
{
public:
    typedef Value val_type;
    typedef Index index_type;
    static_assert(not std::is_const< val_type >::value, "Storage instantiated with const type");
    typedef Pointer< val_type, index_type > ptr_type;
    typedef Pointer< const val_type, index_type > const_ptr_type;
private:
    typedef Identifier< val_type, index_type > idn_type;
    typedef Value_Wrapper< val_type, index_type > wrapper_type;

public:
    // disable copy & move
    Storage(bool activate = true) : _cont(), _next_free_idn() { if (activate) { make_active(); } }
    DELETE_COPY_CTOR(Storage)
    DELETE_MOVE_CTOR(Storage)
    DELETE_COPY_ASOP(Storage)
    DELETE_MOVE_ASOP(Storage)

    void make_active() { active_ptr() = this; }

    /** Get and set active storage. */
    static Storage*& active_ptr() { return _active_ptr; }

    /** Static methods that use active storage. */
    static ptr_type allocate() { ASSERT(active_ptr()); return active_ptr()->ns_allocate(); }
    static void deallocate(ptr_type elem_cptr) { ASSERT(active_ptr()); active_ptr()->ns_deallocate(elem_cptr); }
    static val_type& elem_at(const idn_type& idn) { ASSERT(active_ptr()); return active_ptr()->ns_elem_at(idn); }
    static size_t size() { ASSERT(active_ptr()); return active_ptr()->ns_size(); }
    static size_t unused() { ASSERT(active_ptr()); return active_ptr()->ns_unused(); }

private:
    /** Allocate space for new element. */
    ptr_type ns_allocate()
    {
        ptr_type res;
        if (not _next_free_idn.is_null())
        {
            res._idn = _next_free_idn;
            _next_free_idn = _cont.at(_next_free_idn._idx)._next_free_idn;
        }
        else
        {
            _cont.push_back(wrapper_type());
            res._idn._idx = _cont.size() - 1;
        }
        //std::clog << "allocating element at: " << res._idn._idx << "\n";
        return res;
    }

    /** Deallocate element. */
    void ns_deallocate(ptr_type elem_ptr)
    {
        //std::clog << "deallocating element at: " << elem_ptr._idn._idx << "\n";
        _cont.at(elem_ptr._idn._idx)._next_free_idn = _next_free_idn;
        _next_free_idn = elem_ptr._idn;
    }

    /** Dereference element. */
    val_type& ns_elem_at(const idn_type& idn) const
    {
        return ns_wrapper_at(idn)._key;
    }

    /** Get currently allocated size. */
    size_t ns_size() const { return _cont.size(); }

    /** Get number of unused entries. */
    size_t ns_unused() const
    {
        size_t res = 0;
        for (auto crt = _next_free_idn; not crt.is_null(); crt = ns_wrapper_at(crt)._next_free_idn)
        {
            ++res;
        }
        return res;
    }

    /** Dereference wrapper. */
    wrapper_type& ns_wrapper_at(const idn_type& idn) const
    {
        return const_cast< wrapper_type& >(_cont.at(idn._idx));
    }

    std::deque< wrapper_type > _cont;
    idn_type _next_free_idn;

    static Storage* _active_ptr;

    friend std::ostream& operator << (std::ostream& os, const Storage& f)
    {
        os << "allocated=" << f.ns_size() << '\n'
           << "n_free=" << f.ns_unused() << '\n'
           << "next_free=" << f._next_free_idn << '\n';
        size_t idx = 0;
        for (const auto& w : f._cont)
        {
            os << idx++ << ": [key=" << w._key
            << ",next_free_idn=" << w._next_free_idn << "]\n";
        }
        return os;
    }
}; // class Storage

template < typename Value, typename Index >
Storage< Value, Index >* Storage< Value, Index >::_active_ptr = nullptr;

/** Factory class that manages storage for objects of type T. */
template < typename Value, typename Index = uint32_t >
class Factory
    : public Storage< Value, Index >
{
private:
    typedef Storage< Value, Index > Base;
public:
    typedef Value val_type;
    typedef Index index_type;
    static_assert(not std::is_const< val_type >::value, "Factory instantiated with const type");
    typedef Pointer< val_type, index_type > ptr_type;
    typedef Pointer< const val_type, index_type > const_ptr_type;
    typedef Reference< val_type, index_type > ref_type;
    typedef Reference< const val_type, index_type > const_ref_type;
    typedef Cloner< val_type, index_type > cloner_type;
    typedef Disposer< val_type, index_type > disposer_type;
private:
    typedef Identifier< val_type, index_type > idn_type;
    typedef Value_Wrapper< val_type, index_type > wrapper_type;

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
        ptr_type res = Base::allocate(); // use active storage
        new (res.raw()) val_type(std::forward<Args>(args)...);
        return res;
    }

    /** Destruct and deallocate given element. */
    static void del_elem(ptr_type elem_ptr)
    {
        elem_ptr->~val_type();
        Base::deallocate(elem_ptr); // use active storage
    }
}; // class Factory

/** Stateless allocator that uses the active storage. */
template < typename Value, typename Index = uint32_t >
class Static_Allocator
{
public:
    typedef Value value_type;
    typedef Index index_type;
    typedef Pointer< value_type, index_type > pointer;
    typedef Pointer< const value_type, index_type > const_pointer;
    typedef Reference< value_type, index_type > reference;
    typedef Reference< const value_type, index_type > const_reference;
    typedef Storage< value_type, index_type > storage_type;

    DEFAULT_DEF_CTOR(Static_Allocator)
    DEFAULT_COPY_CTOR(Static_Allocator)
    DEFAULT_MOVE_CTOR(Static_Allocator)
    DEFAULT_COPY_ASOP(Static_Allocator)
    DEFAULT_MOVE_ASOP(Static_Allocator)

    pointer allocate(size_t n, std::allocator< void >::const_pointer = nullptr)
    {
        ASSERT(n == 1);
        return storage_type::allocate(); // use active storage
    }
    void deallocate(pointer p, size_t n)
    {
        ASSERT(n == 1);
        storage_type::deallocate(p); // use active storage
    }
};

/*
template < typename Value, typename Index = uint32_t >
class Holder
{
public:
    typedef T val_type;
    typedef typename std::remove_const< val_type >::type unqual_val_type;
    typedef Factory< unqual_val_type, Index > fact_type;
private:
    typedef Identifier< unqual_val_type, Index > idn_type;
public:
    typedef typename fact_type::ptr_type ptr_type;
    typedef typename fact_type::const_ptr_type const_ptr_type;
    typedef typename fact_type::ref_type ref_type;
    typedef typename fact_type::const_ref_type const_ref_type;

public:
    template <typename... Args>
    Holder(Args&& ... args) { alloc(std::forward<Args>(args)...); }

    // disable copy & move
    DELETE_COPY_CTOR(Holder)
    DELETE_MOVE_CTOR(Holder)
    DELETE_COPY_ASOP(Holder)
    DELETE_MOVE_ASOP(Holder)

    ~Holder() { dealloc(); }

    Holder& operator = (const val_type& rhs) { dealloc(); alloc(rhs); return *this; }

    operator const_ref_type () const { return *_val_ptr; }
    operator ref_type () { return *_val_ptr; }

    const_ptr_type operator & () const { return _val_ptr; }
    ptr_type operator & () { return _val_ptr; }

private:
    template <typename... Args>
    void alloc(Args&& ... args)
    {
        ASSERT(fact_type::get_active_ptr());
        _val_ptr = fact_type::get_active_ptr()->new_elem(std::forward<Args>(args)...);
    }
    void dealloc()
    {
        ASSERT(fact_type::get_active_ptr());
        fact_type::get_active_ptr()->del_elem(_val_ptr);
    }

    ptr_type _val_ptr;
}; // class Holder
*/

} // namespace detail

using detail::Pointer;
using detail::Reference;
using detail::Factory;
using detail::Static_Allocator;
//using detail::Holder;

} // namespace bounded

namespace boost
{
namespace intrusive
{

template < typename Value, typename Index >
struct pointer_traits< bounded::Pointer< Value, Index > >
{
    typedef Value element_type;
    typedef Index index_type;
    typedef bounded::Pointer< element_type, index_type > pointer;
    typedef bounded::Pointer< const element_type, index_type > const_pointer;
    typedef ptrdiff_t difference_type;
    typedef bounded::Reference< element_type, index_type > reference;

    template < class Other_Value >
    struct rebind_pointer
    {
        typedef typename bounded::Pointer< element_type, index_type >::template rebind< Other_Value >::type type;
    };

    static pointer pointer_to(reference r) { return &r; }
    static pointer const_cast_from(const_pointer cptr) { return cptr.unconst(); }
}; // struct pointer_traits

} // namespace intrusive
} // namespace boost


#endif
