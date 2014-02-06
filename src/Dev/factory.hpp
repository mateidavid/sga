#ifndef __FACTORY_HPP
#define __FACTORY_HPP

#include <iostream>
#include <deque>
#include <boost/mpl/if.hpp>
#include <boost/utility/enable_if.hpp>

#include "common.hpp"


namespace factory
{
    namespace detail
    {
        template <class T, class Base_Ptr>
        class Identifier;
        template <class T, bool is_const, class Base_Ptr>
        class Bounded_Pointer;
        template <class T, bool is_const, class Base_Ptr>
        class Bounded_Reference;
        template <class T, class Base_Ptr>
        class Factory;

        template <class T, class Base_Ptr>
        bool operator == (const Identifier< T, Base_Ptr >&, const Identifier< T, Base_Ptr >&);
        template <class T, class Base_Ptr>
        bool operator != (const Identifier< T, Base_Ptr >&, const Identifier< T, Base_Ptr >&);

        template <class T, bool lhs_is_const, bool rhs_is_const, class Base_Ptr>
        bool operator == (const Bounded_Pointer< T, lhs_is_const, Base_Ptr >&,
                          const Bounded_Pointer< T, rhs_is_const, Base_Ptr >&);
        template <class T, bool lhs_is_const, bool rhs_is_const, class Base_Ptr>
        bool operator != (const Bounded_Pointer< T, lhs_is_const, Base_Ptr >&,
                          const Bounded_Pointer< T, rhs_is_const, Base_Ptr >&);
        template <class T, bool is_const, class Base_Ptr>
        std::ostream& operator << (std::ostream& os, const Bounded_Pointer< T, is_const, Base_Ptr >&);

        template <class T, bool is_const, class Base_Ptr>
        std::ostream& operator << (std::ostream& os, const Bounded_Reference< T, is_const, Base_Ptr >&);

        template <class T, class Base_Ptr>
        std::ostream& operator << (std::ostream& os, const Factory< T, Base_Ptr >&);

        template <class T, class Base_Ptr = uint32_t>
        struct Identifier
        {
            typedef T val_type;
            typedef Base_Ptr base_ptr_type;
            typedef Factory< T, Base_Ptr > fact_type;

            Identifier() : _ptr(0) {}
            Identifier(const Identifier& rhs) : _ptr(rhs._ptr) {}

            operator size_t () const { return _ptr; }

            val_type& dereference(const fact_type* fact_cptr = fact_type::get_active_ptr()) const
            {
                ASSERT(fact_cptr);
                return fact_cptr->get_elem(*this);
            }

            Base_Ptr _ptr;
        };

        template <class T, class Base_Ptr>
        bool operator == (const Identifier< T, Base_Ptr >& lhs, const Identifier< T, Base_Ptr >& rhs)
        {
            return lhs._ptr == rhs._ptr;
        }
        template <class T, class Base_Ptr>
        bool operator != (const Identifier< T, Base_Ptr >& lhs, const Identifier< T, Base_Ptr >& rhs)
        {
            return !(lhs == rhs);
        }

        template <class T, bool is_const, class Base_Ptr = uint32_t>
        class Bounded_Pointer
        {
        private:
            struct enabler {};
        public:
            typedef T val_type;
            typedef val_type value_type;
            typedef typename boost::mpl::if_c< is_const, const val_type&, val_type& >::type qual_real_ref_type;
            typedef Base_Ptr base_ptr_type;

            Bounded_Pointer() : _id() {}
            Bounded_Pointer(const Bounded_Pointer& rhs) : _id(rhs._id) {}
            explicit Bounded_Pointer(int zero) : _id() { ASSERT(zero == 0); }

            template <bool _unused = true>
            explicit Bounded_Pointer(const Bounded_Pointer< T, true, Base_Ptr >& rhs,
                                     typename boost::enable_if_c< _unused and not is_const, enabler >::type = enabler())
                : _id(rhs._id) {}

            template <bool _unused = true>
            Bounded_Pointer(const Bounded_Pointer< T, false, Base_Ptr >& rhs,
                            typename boost::enable_if_c< _unused and is_const, enabler >::type = enabler())
                : _id(rhs._id) {}

            operator bool() const { return _id; }
            Bounded_Pointer& operator ++ () { ++_id._ptr; return *this; }
            Bounded_Pointer operator ++ (int) { Bounded_Pointer res(*this); ++(*this); return res; }

            Bounded_Reference< T, is_const, Base_Ptr > operator * () const
            { return Bounded_Reference< T, is_const, Base_Ptr >(_id); }
            qual_real_ref_type operator -> () const { return _id.dereference(); }

            template <class Other_T, bool lhs_is_const, bool rhs_is_const, class Other_Base_Ptr>
            friend bool operator == (const Bounded_Pointer< Other_T, lhs_is_const, Other_Base_Ptr >&,
                                     const Bounded_Pointer< Other_T, rhs_is_const, Other_Base_Ptr >&);
            template <class Other_T, bool lhs_is_const, bool rhs_is_const, class Other_Base_Ptr>
            friend bool operator != (const Bounded_Pointer< Other_T, lhs_is_const, Other_Base_Ptr >&,
                                     const Bounded_Pointer< Other_T, rhs_is_const, Other_Base_Ptr >&);
            friend std::ostream& operator << <>(std::ostream&, const Bounded_Pointer< T, is_const, Base_Ptr >&);

        private:
            friend class Bounded_Pointer< T, not is_const, Base_Ptr >;
            friend class Bounded_Reference< T, is_const, Base_Ptr >;
            friend class Factory< T, Base_Ptr >;

            Bounded_Pointer(const Identifier< T, Base_Ptr >& id) : _id(id) {}

            Identifier< T, Base_Ptr > _id;
        };

        template <class Other_T, bool lhs_is_const, bool rhs_is_const, class Other_Base_Ptr>
        bool operator == (const Bounded_Pointer< Other_T, lhs_is_const, Other_Base_Ptr >& lhs,
                          const Bounded_Pointer< Other_T, rhs_is_const, Other_Base_Ptr >& rhs)
        {
            return lhs._id == rhs._id;
        }
        template <class Other_T, bool lhs_is_const, bool rhs_is_const, class Other_Base_Ptr>
        bool operator != (const Bounded_Pointer< Other_T, lhs_is_const, Other_Base_Ptr >& lhs,
                          const Bounded_Pointer< Other_T, rhs_is_const, Other_Base_Ptr >& rhs)
        {
            return !(lhs == rhs);
        }

        template <class T, bool is_const, class Base_Ptr>
        std::ostream& operator <<(std::ostream& os, const Bounded_Pointer< T, is_const, Base_Ptr >& rhs)
        {
            os << "[Bounded_Pointer:&=" << (void*)&rhs << ",_ptr=" << rhs._id._ptr << "]";
            return os;
        }

        /** Bounded Reference. */
        template <class T, bool is_const, class Base_Ptr = uint32_t>
        class Bounded_Reference
        {
        public:
            typedef T val_type;
            typedef typename boost::mpl::if_c< is_const, const val_type&, val_type& >::type qual_real_ref_type;
            typedef Bounded_Pointer< T, is_const, Base_Ptr > qual_ptr_type;
        private:
            typedef Identifier< T, Base_Ptr > idn_type;

        public:
            Bounded_Reference(const Bounded_Reference& rhs) : _id(rhs._id) {}

            qual_ptr_type operator & () const { return qual_ptr_type(_id); }

            operator qual_real_ref_type () const { return _id.dereference(); }

            const Bounded_Reference& operator = (const val_type& rhs) const
            { _id.dereference() = rhs; return *this; }
            const Bounded_Reference& operator = (const Bounded_Reference& rhs) const
            { _id.dereference() = rhs._id.dereference(); return *this; }

        private:
            friend class Bounded_Pointer< T, is_const, Base_Ptr >;
            friend class Bounded_Reference< T, true, Base_Ptr >;
            friend std::ostream& operator << <>(std::ostream&, const Bounded_Reference&);

            Bounded_Reference(const idn_type& id) : _id(id) {}

            const idn_type _id;
        };

        template <class T, class Base_Ptr>
        class Bounded_Reference< T, true, Base_Ptr >
        {
        public:
            typedef T val_type;
            typedef const val_type& qual_real_ref_type;
            typedef Bounded_Pointer< T, true, Base_Ptr > qual_ptr_type;
        private:
            typedef Identifier< T, Base_Ptr > idn_type;

        public:
            Bounded_Reference(const Bounded_Reference& rhs) : _id(rhs._id) {}
            Bounded_Reference(const Bounded_Reference< T, false, Base_Ptr >& rhs) : _id(rhs._id) {}

            qual_ptr_type operator & () const { return qual_ptr_type(_id); }

            operator qual_real_ref_type () const { return _id.dereference(); }

        private:
            friend class Bounded_Pointer< T, true, Base_Ptr >;
            friend std::ostream& operator << <>(std::ostream&, const Bounded_Reference&);

            Bounded_Reference(const idn_type& id) : _id(id) {}
            Bounded_Reference& operator = (const Bounded_Reference& rhs) { return *this; }

            const idn_type _id;
        };

        template <class T, bool is_const, class Base_Ptr>
        std::ostream& operator <<(std::ostream& os, const Bounded_Reference< T, is_const, Base_Ptr >& rhs)
        {
            os << "[Bounded_Reference:&=" << (void*)&rhs._id << ",_ptr=" << rhs._id._ptr
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
            /** Non-constant pointer. */
            typedef Bounded_Pointer< T, false, Base_Ptr > ptr_type;
            /** Constant pointer. */
            typedef Bounded_Pointer< T, true, Base_Ptr > const_ptr_type;
            /** Non-constant reference. */
            typedef Bounded_Reference< T, false, Base_Ptr > ref_type;
            /** Constant reference. */
            typedef Bounded_Reference< T, true, Base_Ptr > const_ref_type;
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
                std::clog << "allocating space for new element at: " << (void*)&_cont.at(res._id._ptr) << '\n';
                new (&_cont.at(res._id._ptr)) val_type(std::forward<Args>(args)...);
                return res;
            }

            /** Delete element pointed at. */
            void del_elem(const_ptr_type elem_ptr)
            {
                _cont.at(elem_ptr._id._ptr)._next_free_idn = _next_free_idn;
                _next_free_idn = elem_ptr._id;
            }

            /** Use this object to resolve all managed pointers. */
            void set_active() { _active_ptr = this; }

            /** Get active factory pointer. */
            static Factory* get_active_ptr() { return _active_ptr; }

        private:
            friend class Identifier< T, Base_Ptr >;

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
            typedef Bounded_Pointer< T, false, Base_Ptr > ptr_type;
            typedef Bounded_Pointer< T, true, Base_Ptr > const_ptr_type;
            typedef Bounded_Reference< T, false, Base_Ptr > ref_type;
            typedef Bounded_Reference< T, true, Base_Ptr > const_ref_type;
        private:
            typedef Factory< T, Base_Ptr > fact_type;
            typedef Identifier< T, Base_Ptr > idn_type;

        public:
            template <typename... Args>
            Holder(Args&&... args) { alloc(std::forward<Args>(args)...); }
            Holder(const Holder& rhs) { alloc(rhs); }
            ~Holder() { dealloc(); }

            Holder& operator = (const val_type& rhs) { dealloc(); alloc(rhs); return *this; }
            Holder& operator = (const Holder& rhs) { dealloc(); alloc(rhs); return *this; }

            operator const_ref_type () const { return *_val_ptr; }
            operator ref_type () { return *_val_ptr; }

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

namespace std
{
    template <class T, bool is_const, class Base_Ptr >
    struct iterator_traits< factory::Bounded_Pointer< T, is_const, Base_Ptr > >
    {
        typedef T value_type;
        typedef factory::Bounded_Pointer< T, is_const, Base_Ptr > pointer;
        typedef ptrdiff_t difference_type;
        typedef factory::Bounded_Reference< T, is_const, Base_Ptr > reference;
        typedef std::forward_iterator_tag iterator_category;
    };
}

namespace boost
{
    template <class T, bool is_const, class Base_Ptr >
    struct pointer_to_other< factory::Bounded_Pointer< T, is_const, Base_Ptr >, void >
    {
        typedef factory::Bounded_Pointer< T, is_const, Base_Ptr > type;
    };
}

namespace boost { namespace intrusive { namespace detail
{
    template <class T, bool is_const, class Base_Ptr >
    struct smart_ptr_type< factory::Bounded_Pointer< T, is_const, Base_Ptr > >
    {
        typedef T value_type;
        typedef factory::Bounded_Pointer< T, is_const, Base_Ptr > pointer;
        static pointer get(const pointer& ptr) { return ptr; }
    };
}}}

#endif
