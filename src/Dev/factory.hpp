#ifndef __FACTORY_HPP
#define __FACTORY_HPP

#include <iostream>
#include <deque>
#include <type_traits>
#include <boost/intrusive/pointer_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/utility/enable_if.hpp>
#include <typeinfo>

#include "common.hpp"


namespace factory
{
    namespace detail
    {
        template <class LHS_T, class RHS_T>
        bool operator != (const LHS_T& lhs, const RHS_T& rhs)
        {
            return !(lhs == rhs);
        }

        template <class T, class Base_Ptr>
        struct Identifier;
        template <class T, class Base_Ptr>
        class Bounded_Pointer;
        template <class T, class Base_Ptr>
        class Bounded_Reference;
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
            typedef Base_Ptr base_ptr_type;
            typedef typename boost::mpl::if_c< not is_void,
                                               Factory< T, Base_Ptr >,
                                               void >::type fact_type;

            Identifier() : _ptr(0) {}
            Identifier(const Identifier& rhs) : _ptr(rhs._ptr) {}
            template <class U>
            explicit Identifier(const Identifier< U, Base_Ptr >& rhs) : _ptr(rhs._ptr) {}

            operator size_t () const { return _ptr; }

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
        protected:
            typedef Identifier< unqual_val_type, Base_Ptr > idn_type;
            static const bool is_void = idn_type::is_void;
        public:
            typedef typename idn_type::fact_type fact_type;
            typedef typename boost::mpl::if_c< not is_void, Bounded_Reference< val_type, Base_Ptr >, void >::type ref_type;

            // rebinder
            template <class U>
            struct rebind
            {
                typedef typename boost::mpl::if_c<
                    std::is_same< std::remove_const< U >, std::remove_const< T > >::value,
                    Bounded_Pointer< U, Base_Ptr >,
                    U*
                    >::type type;
            };

        public:
            Bounded_Pointer() {}
            Bounded_Pointer(const Bounded_Pointer& rhs) : _id(rhs._id) {}
            //explicit Bounded_Pointer(int zero) : _id() { ASSERT(zero == 0); }

            // implicit conversion to const
            operator Bounded_Pointer< const unqual_val_type, Base_Ptr >& ()
            { return *reinterpret_cast<Bounded_Pointer< const unqual_val_type, Base_Ptr >* >(this); }
            operator const Bounded_Pointer< const unqual_val_type, Base_Ptr >& () const
            { return *reinterpret_cast<const Bounded_Pointer< const unqual_val_type, Base_Ptr >* >(this); }

            // unconst
            explicit operator Bounded_Pointer< unqual_val_type, Base_Ptr >& ()
            { return *reinterpret_cast<Bounded_Pointer< unqual_val_type, Base_Ptr >*>(this); }
            explicit operator const Bounded_Pointer< unqual_val_type, Base_Ptr >& () const
            { return *reinterpret_cast<const Bounded_Pointer< unqual_val_type, Base_Ptr >*>(this); }

            // conversion to void*
            operator typename boost::mpl::if_c<
                std::is_same< val_type, unqual_val_type >::value,
                void*,
                const void* >::type () const
            { return &_id.dereference(); }

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
            Bounded_Reference(const Bounded_Reference& rhs) : _id(rhs._id) {}

            // automatic conversion to const
            operator Bounded_Reference< const unqual_val_type, Base_Ptr >& ()
            { return *reinterpret_cast< Bounded_Reference< const unqual_val_type, Base_Ptr >* >(this); }
            operator const Bounded_Reference< const unqual_val_type, Base_Ptr >& () const
            { return *reinterpret_cast< const Bounded_Reference< const unqual_val_type, Base_Ptr >* >(this); }

            // explicit unconst
            explicit operator Bounded_Reference< unqual_val_type, Base_Ptr >& ()
            { return *reinterpret_cast< Bounded_Reference< unqual_val_type, Base_Ptr >* >(this); }
            explicit operator const Bounded_Reference< unqual_val_type, Base_Ptr >& () const
            { return *reinterpret_cast< const Bounded_Reference< unqual_val_type, Base_Ptr >* >(this); }

            ptr_type operator & () const { return ptr_type(_id); }

            operator real_ref_type () const { return _id.dereference(); }

            const Bounded_Reference& operator = (
                typename boost::mpl::if_c< std::is_same< val_type, unqual_val_type>::value,
                                           const val_type&,
                                           enabler<0>& >::type rhs) const
            { _id.dereference() = rhs; return *this; }

            const Bounded_Reference& operator = (
                typename boost::mpl::if_c< std::is_same< val_type, unqual_val_type>::value,
                                           const Bounded_Reference&,
                                           enabler<1>& >::type rhs) const
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

        /** Wrapper type used by Factory objects to store free node list without overhead. */
        template <class T, class Base_Ptr = uint32_t>
        union Factory_Wrapper
        {
            typedef T val_type;
            typedef Identifier< T, Base_Ptr > idn_type;

            val_type _key;
            idn_type _next_free_idn;

            Factory_Wrapper() : _next_free_idn() {}
            Factory_Wrapper(const Factory_Wrapper& rhs)
                : _next_free_idn(rhs._next_free_idn) {}
        };

        /** Factory class that manages storage for objects of type T. */
        template <class T, class Base_Ptr = uint32_t>
        class Factory
        {
        public:
            /** Type of object stored. */
            typedef T val_type;
            typedef const T const_val_type;
            /** Non-constant pointer. */
            typedef Bounded_Pointer< val_type, Base_Ptr > ptr_type;
            /** Constant pointer. */
            typedef Bounded_Pointer< const_val_type, Base_Ptr > const_ptr_type;
            /** Non-constant reference. */
            typedef Bounded_Reference< val_type, Base_Ptr > ref_type;
            /** Constant reference. */
            typedef Bounded_Reference< const_val_type, Base_Ptr > const_ref_type;
        private:
            typedef Identifier< T, Base_Ptr > idn_type;
            typedef Factory_Wrapper< T, Base_Ptr > wrapper_type;

        public:
            /** Default constructor.
             * @param activate Bool; if true, the object is used to resolve all managed pointers.
             */
            Factory(bool activate = false) : _cont(1, wrapper_type()), _next_free_idn() { if (activate) set_active(); }

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
            ptr_type new_elem(Args&&... args)
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
                std::clog << "allocating element at: " << (void*)&_cont.at(res._id._ptr) << '\n';
                new (&_cont.at(res._id._ptr)) val_type(std::forward<Args>(args)...);
                return res;
            }

            /** Delete element pointed at. */
            void del_elem(const_ptr_type elem_ptr)
            {
                std::clog << "deallocating element at: " << (void*)&_cont.at(elem_ptr._id._ptr) << '\n';
                _cont.at(elem_ptr._id._ptr)._next_free_idn = _next_free_idn;
                _next_free_idn = elem_ptr._id;
            }

            /** Use this object to resolve all managed pointers. */
            void set_active() { _active_ptr = this; }

            /** Get active factory pointer. */
            static Factory* get_active_ptr() { return _active_ptr; }

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
        std::ostream& operator << (std::ostream& os, const Factory< T, Base_Ptr >& rhs)
        {
            os << "allocated=" << rhs._cont.size() << '\n';
            size_t n_free = 0;
            for (auto crt = rhs._next_free_idn;
                 crt;
                 crt = rhs.dereference(crt)._next_free_idn)
            {
                ++n_free;
            }
            os << "n_free=" << n_free << '\n'
               << "next_free=" << rhs._next_free_idn << '\n';
            for (auto it = rhs._cont.begin(); it != rhs._cont.end(); ++it)
            {
                os << "[key=" << it->_key
                   << ",next_free_idn=" << it->_next_free_idn << "]\n";
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
            Holder(Args&&... args) { alloc(std::forward<Args>(args)...); }
            Holder(const Holder& rhs) { alloc(rhs); }
            ~Holder() { dealloc(); }

            Holder& operator = (const val_type& rhs) { dealloc(); alloc(rhs); return *this; }
            Holder& operator = (const Holder& rhs) { dealloc(); alloc(rhs); return *this; }

            operator const_ref_type () const { return *_val_ptr; }
            operator ref_type () { return *_val_ptr; }
            operator const val_type& () const { return (const val_type&)(*_val_ptr); }
            operator val_type& () { return (val_type&)(*_val_ptr); }

            const_ptr_type operator & () const { return _val_ptr; }
            ptr_type operator & () { return _val_ptr; }

        private:
            template <typename... Args>
            void alloc(Args&&... args) { _val_ptr = fact_type::get_active_ptr()->new_elem(std::forward<Args>(args)...); }
            void alloc(const Holder& rhs) { alloc((val_type)rhs); }
            void dealloc() { fact_type::get_active_ptr()->del_elem(_val_ptr); }

            ptr_type _val_ptr;
        };
    }
    using detail::Bounded_Pointer;
    using detail::Bounded_Reference;
    using detail::Factory;
    using detail::Holder;
}


namespace boost { namespace intrusive
{
    template <class T, class Base_Ptr >
    struct pointer_traits< factory::Bounded_Pointer< T, Base_Ptr > >
    {
        typedef T element_type;
        typedef factory::Bounded_Pointer< T, Base_Ptr > pointer;
        typedef ptrdiff_t difference_type;
        typedef factory::Bounded_Reference< T, Base_Ptr > reference;

        template <class U>
        struct rebind_pointer
        {
            typedef typename factory::Bounded_Pointer< T, Base_Ptr >::template rebind<U>::type type;
        };

        static pointer pointer_to(reference r)
        {
            std::clog << "pointer_to from [" << typeid(reference).name()
                      << "] to [Ptr<" << typeid(T).name() << ">]\n";
            return &r;
        }
        static pointer pointer_to(const factory::Holder< T, Base_Ptr >& r)
        {
            std::clog << "pointer_to from ["
                      << typeid(const factory::Holder< T, Base_Ptr >&).name()
                      << "] to [Ptr<" << typeid(T).name() << ">]\n";
            return pointer(&r);
        }

        static pointer const_cast_from(
            const factory::Bounded_Pointer<
                const typename factory::Bounded_Pointer< T, Base_Ptr >::unqual_val_type,
                Base_Ptr >& cptr)
        {
            std::clog << "const cast to [Ptr<" << typeid(T).name() << ">]\n";
            return pointer(cptr);
        }
    };
}}


#endif
