#ifndef __FACTORY_HPP
#define __FACTORY_HPP

#include <cstddef>
#include <iostream>
#include <deque>
#include <type_traits>
#include <boost/intrusive/pointer_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/utility/enable_if.hpp>
#include <typeinfo>

#include "shortcuts.hpp"
#include "global_assert.hpp"
#include "nonconst_methods.hpp"


namespace detail
{

template <class T, class Base_Ptr>
struct Identifier;
template <class T, class Base_Ptr>
class Bounded_Pointer;
template <class T, class Base_Ptr>
class Bounded_Reference;
template <class T, class Base_Ptr>
struct Cloner;
template <class T, class Base_Ptr>
struct Disposer;
template <class T, class Base_Ptr>
class Factory;

template <class T, class Base_Ptr>
bool operator == (const Identifier< T, Base_Ptr >&,
                  const Identifier< T, Base_Ptr >&);

template <class LHS_T, class RHS_T, class Base_Ptr>
bool operator == (const Bounded_Pointer< LHS_T, Base_Ptr >&,
                  const Bounded_Pointer< RHS_T, Base_Ptr >&);
template <class T, class Base_Ptr>
std::ostream& operator << (std::ostream&, const Bounded_Pointer< T, Base_Ptr >&);

template <class T, class Base_Ptr>
std::ostream& operator << (std::ostream& os, const Bounded_Reference< T, Base_Ptr >&);

template <class T, class Base_Ptr>
std::ostream& operator << (std::ostream& os, const Factory< T, Base_Ptr >&);

template <class T, class Base_Ptr = uint32_t>
struct Identifier
{
private:
    struct enabler {};
public:
    static const bool is_void = std::is_same< T, void >::value;

    typedef T val_type;
    static_assert(std::is_same< val_type, typename std::remove_const< val_type >::type >::value, "Identifier instantiated with const type");
    typedef Base_Ptr base_ptr_type;
    typedef typename boost::mpl::if_c< not is_void, Factory< T, Base_Ptr >, void >::type fact_type;

    Identifier() : _ptr(0) {}

    // allow copy & move
    DEFAULT_COPY_CTOR(Identifier)
    DEFAULT_MOVE_CTOR(Identifier)
    DEFAULT_COPY_ASOP(Identifier)
    DEFAULT_MOVE_ASOP(Identifier)

    operator Base_Ptr () const { return _ptr; }

    template <bool _unused = true>
    typename boost::mpl::if_c< not is_void, val_type, enabler >::type& dereference(
        typename boost::enable_if_c< _unused and not is_void, enabler >::type = enabler()) const
    {
        ASSERT(fact_type::get_active_ptr());
        return fact_type::get_active_ptr()->get_elem(*this);
    }

    Base_Ptr _ptr;
};

template <class T, class Base_Ptr>
bool operator == (const Identifier< T, Base_Ptr >& lhs,
                  const Identifier< T, Base_Ptr >& rhs)
{
    return lhs._ptr == rhs._ptr;
}

template <class T, class Base_Ptr = uint32_t>
class Bounded_Pointer
{
public:
    typedef T val_type;
    typedef Base_Ptr base_ptr_type;
    typedef val_type* real_ptr_type;
    typedef typename std::remove_const< val_type >::type unqual_val_type;
private:
    typedef Identifier< unqual_val_type, Base_Ptr > idn_type;
public:
    typedef typename idn_type::fact_type fact_type;
    typedef typename boost::mpl::if_c< not idn_type::is_void, Bounded_Reference< val_type, Base_Ptr >, void >::type ref_type;

    /** Rebinder.
     * Wrap pointers to T, const T, void, and const void. Anything else gets transformed into a raw pointer.
     */
    template <class U>
    struct rebind
    {
        typedef typename boost::mpl::if_c<
        std::is_same< typename std::remove_const< U >::type, typename std::remove_const< T >::type >::value
#ifdef WRAP_VOID
        or std::is_void< U >::value
#endif
        ,
        Bounded_Pointer< U, Base_Ptr >,
        U*
        >::type type;
    };

public:
    Bounded_Pointer() {}
    Bounded_Pointer(std::nullptr_t) {}

    // allow copy & move
    DEFAULT_COPY_CTOR(Bounded_Pointer)
    DEFAULT_MOVE_CTOR(Bounded_Pointer)
    DEFAULT_COPY_ASOP(Bounded_Pointer)
    DEFAULT_MOVE_ASOP(Bounded_Pointer)

    // implicit conversion to const
    operator const Bounded_Pointer< const unqual_val_type, Base_Ptr >& () const
    { return *reinterpret_cast<const Bounded_Pointer< const unqual_val_type, Base_Ptr >* >(this); }
    DEF_NONCONST_CONVERSION((Bounded_Pointer< const unqual_val_type, Base_Ptr >&))

    // explicit conversion to non-const
    explicit operator const Bounded_Pointer< unqual_val_type, Base_Ptr >& () const
    { return *reinterpret_cast< const Bounded_Pointer< unqual_val_type, Base_Ptr >* >(this); }
    explicit DEF_NONCONST_CONVERSION((Bounded_Pointer< unqual_val_type, Base_Ptr >&))

    // conversion to void*
    operator typename boost::mpl::if_c< std::is_const< val_type >::value, const void*, void* >::type () const
    { return &_id.dereference(); }

    // get raw pointer
    val_type* raw() const { return &_id.dereference(); }

    operator bool() const { return this->_id; }
    Bounded_Pointer& operator ++ () { ++(this->_id)._ptr; return *this; }
    Bounded_Pointer operator ++ (int) { Bounded_Pointer res(*this); ++(*this); return res; }

    // dereferencing operators
    ref_type operator * () const { return ref_type(_id); }
    real_ptr_type operator -> () const { return &_id.dereference(); }

private:
    friend class Bounded_Reference< val_type, Base_Ptr >;
    friend class Factory< unqual_val_type, Base_Ptr >;

    template <class LHS_T, class RHS_T, class Other_Base_Ptr>
    friend bool operator == (const Bounded_Pointer< LHS_T, Other_Base_Ptr >&,
                             const Bounded_Pointer< RHS_T, Other_Base_Ptr >&);
    friend std::ostream& operator << <>(std::ostream&, const Bounded_Pointer&);

    Bounded_Pointer(const idn_type& id) : _id(id) {}

    idn_type _id;
};

template <class LHS_T, class RHS_T, class Other_Base_Ptr>
bool operator == (const Bounded_Pointer< LHS_T, Other_Base_Ptr >& lhs,
                  const Bounded_Pointer< RHS_T, Other_Base_Ptr >& rhs)
{
    return lhs._id == rhs._id;
}
template <class T, class Base_Ptr>
std::ostream& operator <<(std::ostream& os, const Bounded_Pointer< T, Base_Ptr >& rhs)
{
    os << "[Bounded_Pointer:&=" << (void*)&rhs
       << ",_ptr=" << rhs._id._ptr << "]";
    return os;
}

/** Bounded Reference. */
template <class T, class Base_Ptr = uint32_t>
class Bounded_Reference
{
private:
    template <int val> struct enabler {};
public:
    typedef T val_type;
    typedef Base_Ptr base_ptr_type;
    typedef val_type& real_ref_type;
    typedef typename std::remove_const< val_type >::type unqual_val_type;
    typedef Bounded_Pointer< val_type, Base_Ptr > ptr_type;
private:
    typedef Identifier< unqual_val_type, Base_Ptr > idn_type;

public:
    // allow construction only
    DEFAULT_COPY_CTOR(Bounded_Reference)
    DEFAULT_MOVE_CTOR(Bounded_Reference)
    DELETE_COPY_ASOP(Bounded_Reference)
    DELETE_MOVE_ASOP(Bounded_Reference)

    // automatic conversion to const
    operator const Bounded_Reference< const unqual_val_type, Base_Ptr >& () const
    { return *reinterpret_cast< const Bounded_Reference< const unqual_val_type, Base_Ptr >* >(this); }
    DEF_NONCONST_CONVERSION(( Bounded_Reference< const unqual_val_type, Base_Ptr >&))

    // explicit conversion to nonconst
    explicit operator const Bounded_Reference< unqual_val_type, Base_Ptr >& () const
    { return *reinterpret_cast< const Bounded_Reference< unqual_val_type, Base_Ptr >* >(this); }
    explicit DEF_NONCONST_CONVERSION(( Bounded_Reference< unqual_val_type, Base_Ptr >&))

    ptr_type operator & () const { return ptr_type(_id); }

    // get raw reference
    real_ref_type raw() const { return _id.dereference(); }
    operator real_ref_type () const { return raw(); }

    const Bounded_Reference& operator = (
        typename boost::mpl::if_c< std::is_same< val_type, unqual_val_type>::value,
                                   const val_type&,
                                   enabler<0>&
                                 >::type rhs) const
    { _id.dereference() = rhs; return *this; }

    const Bounded_Reference& operator = (
        typename boost::mpl::if_c< std::is_same< val_type, unqual_val_type>::value,
                                   const Bounded_Reference&,
                                   enabler<1>&
                                 >::type rhs) const
    { _id.dereference() = rhs._id.dereference(); return *this; }

private:
    friend class Bounded_Pointer< T, Base_Ptr >;
    friend std::ostream& operator << <>(std::ostream&, const Bounded_Reference&);

    Bounded_Reference(const idn_type& id) : _id(id) {}

    const idn_type _id;
};

template <class T, class Base_Ptr>
std::ostream& operator <<(std::ostream& os, const Bounded_Reference< T, Base_Ptr >& rhs)
{
    os << "[Bounded_Reference:&=" << (void*)&rhs._id
       << ",_ptr=" << rhs._id._ptr
       << ",deref=" << rhs._id.dereference() << "]";
    return os;
}

template <class T, class Base_Ptr = uint32_t>
struct Cloner
{
    typedef T val_type;
    static_assert(std::is_same< val_type, typename std::remove_const< val_type >::type >::value, "Cloner instantiated with const type");
    typedef Identifier< val_type, Base_Ptr > idn_type;
    typedef Bounded_Pointer< val_type, Base_Ptr > ptr_type;
    typedef typename idn_type::fact_type fact_type;

    ptr_type operator () (const T& t)
    {
        ASSERT(fact_type::get_active_ptr());
        return fact_type::get_active_ptr()->new_elem(t);
    }
};

template <class T, class Base_Ptr = uint32_t>
struct Disposer
{
    typedef T val_type;
    static_assert(std::is_same< val_type, typename std::remove_const< val_type >::type >::value, "Disposer instantiated with const type");
    typedef Identifier< val_type, Base_Ptr > idn_type;
    typedef Bounded_Pointer< val_type, Base_Ptr > ptr_type;
    typedef typename idn_type::fact_type fact_type;

    void operator () (ptr_type ptr)
    {
        ASSERT(fact_type::get_active_ptr());
        fact_type::get_active_ptr()->del_elem(ptr);
    }
};

/** Wrapper type used by Factory objects to store free node list without overhead. */
template <class T, class Base_Ptr = uint32_t>
union Factory_Wrapper
{
    typedef T val_type;
    typedef Identifier< T, Base_Ptr > idn_type;

    val_type _key;
    idn_type _next_free_idn;

    Factory_Wrapper() : _next_free_idn() {}
    Factory_Wrapper(const Factory_Wrapper& rhs) : _next_free_idn(rhs._next_free_idn) {}
    ~Factory_Wrapper() { (&(this->_next_free_idn))->~idn_type(); }
};

/** Factory class that manages storage for objects of type T. */
template <class T, class Base_Ptr = uint32_t>
class Factory
{
public:
    /** Type of object stored. */
    typedef T val_type;
    static_assert(std::is_same< val_type, typename std::remove_const< val_type >::type >::value, "Factory instantiated with const type");
    typedef const T const_val_type;
    /** Non-constant pointer. */
    typedef Bounded_Pointer< val_type, Base_Ptr > ptr_type;
    /** Constant pointer. */
    typedef Bounded_Pointer< const_val_type, Base_Ptr > const_ptr_type;
    /** Non-constant reference. */
    typedef Bounded_Reference< val_type, Base_Ptr > ref_type;
    /** Constant reference. */
    typedef Bounded_Reference< const_val_type, Base_Ptr > const_ref_type;
    /** Object cloner. */
    typedef Cloner< T, Base_Ptr > cloner;
    /** Object disposer. */
    typedef Disposer< T, Base_Ptr > disposer;
private:
    typedef Identifier< T, Base_Ptr > idn_type;
    typedef Factory_Wrapper< T, Base_Ptr > wrapper_type;

public:
    /** Default constructor.
     * @param activate Bool; if true, the object is used to resolve all managed pointers.
     */
    Factory(bool activate = true) : _cont(1, wrapper_type()), _next_free_idn() { if (activate) { set_active(); } }

    // disable copy & move
    DELETE_COPY_CTOR(Factory)
    DELETE_MOVE_CTOR(Factory)
    DELETE_COPY_ASOP(Factory)
    DELETE_MOVE_ASOP(Factory)

private:
    /** Dereference identifier. */
    wrapper_type& dereference(const idn_type& idn) const
    {
        return const_cast<wrapper_type&>(_cont.at(idn._ptr));
    }
    /** Get element pointed at. */
    val_type& get_elem(const idn_type& idn) const
    {
        return dereference(idn)._key;
    }

public:
    /** Allocate space for new element and return pointer to it. */
    template <typename... Args>
    ptr_type new_elem_ns(Args&& ... args)
    {
        ptr_type res;
        if (_next_free_idn)
        {
            res._id = _next_free_idn;
            _next_free_idn = _cont.at(_next_free_idn._ptr)._next_free_idn;
        }
        else
        {
            _cont.push_back(wrapper_type());
            res._id._ptr = _cont.size() - 1;
        }
        //std::clog << "allocating element at: " << (void*)&_cont.at(res._id._ptr) << '\n';
        new (&_cont.at(res._id._ptr)) val_type(std::forward<Args>(args)...);
        return res;
    }

    /** Delete element pointed at. */
    void del_elem_ns(ptr_type elem_ptr)
    {
        //std::clog << "deallocating element at: " << (void*)&_cont.at(elem_ptr._id._ptr) << '\n';
        elem_ptr->~val_type();
        _cont.at(elem_ptr._id._ptr)._next_free_idn = _next_free_idn;
        _next_free_idn = elem_ptr._id;
    }

    /** Get currently allocated size. */
    size_t size() const { return _cont.size(); }

    /** Get number of unused entries. */
    size_t unused() const
    {
        size_t res = 0;
        for (auto crt = _next_free_idn; crt; crt = dereference(crt)._next_free_idn)
        {
            ++res;
        }
        return res;
    }

    /** Use this object to resolve all managed pointers. */
    void set_active() { _active_ptr = this; }

    /** Get active factory pointer. */
    static Factory* get_active_ptr() { return _active_ptr; }

    /** Static version of new_elem, using active factory. */
    template <typename... Args>
    static ptr_type new_elem(Args&& ... args)
    {
        ASSERT(get_active_ptr());
        return get_active_ptr()->new_elem_ns(std::forward< Args >(args)...);
    }

    /** Static version of del_elem, using active factory. */
    template <typename... Args>
    static void del_elem(Args&& ... args)
    {
        ASSERT(get_active_ptr());
        get_active_ptr()->del_elem_ns(std::forward< Args >(args)...);
    }

private:
    friend struct Identifier< T, Base_Ptr >;

    static Factory* _active_ptr;

    friend std::ostream& operator << <>(std::ostream&, const Factory&);

    std::deque< wrapper_type > _cont;
    idn_type _next_free_idn;
};

template <class T, class Base_Ptr>
Factory< T, Base_Ptr >* Factory< T, Base_Ptr >::_active_ptr = NULL;

template <class T, class Base_Ptr>
std::ostream& operator << (std::ostream& os, const Factory< T, Base_Ptr >& f)
{
    os << "allocated=" << f._cont.size() << '\n';
    size_t n_free = 0;
    for (auto crt = f._next_free_idn; crt; crt = f.dereference(crt)._next_free_idn)
    {
        ++n_free;
    }
    os << "n_free=" << n_free << '\n'
       << "next_free=" << f._next_free_idn << '\n';
    size_t idx = 0;
    for (const auto& w : f._cont)
    {
        os << idx++ << ": [key=" << w._key
           << ",next_free_idn=" << w._next_free_idn << "]\n";
    }
    return os;
}

template <class T, class Base_Ptr = uint32_t>
class Holder
{
public:
    typedef T val_type;
    typedef typename std::remove_const< val_type >::type unqual_val_type;
    typedef Factory< unqual_val_type, Base_Ptr > fact_type;
private:
    typedef Identifier< unqual_val_type, Base_Ptr > idn_type;
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
};

} // namespace detail

using detail::Bounded_Pointer;
using detail::Bounded_Reference;
using detail::Factory;
using detail::Holder;

namespace boost
{
namespace intrusive
{

template <class T, class Base_Ptr >
struct pointer_traits< Bounded_Pointer< T, Base_Ptr > >
{
    typedef T element_type;
    typedef Bounded_Pointer< T, Base_Ptr > pointer;
    typedef Bounded_Pointer< const T, Base_Ptr > const_pointer;
    typedef ptrdiff_t difference_type;
    typedef Bounded_Reference< T, Base_Ptr > reference;

    template <class U>
    struct rebind_pointer
    {
        typedef typename Bounded_Pointer< T, Base_Ptr >::template rebind<U>::type type;
    };

    static pointer pointer_to(reference r)
    {
        return &r;
    }

    static pointer const_cast_from(const const_pointer& cptr)
    {
        return pointer(cptr);
    }
};

} // namespace intrusive
} // namespace boost


#endif
