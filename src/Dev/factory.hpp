#ifndef __FACTORY_HPP
#define __FACTORY_HPP

#include <iostream>

#include "common.hpp"


namespace factory
{
    namespace detail
    {
        template <class T>
        class Factory;
        template <class T>
        class Identifier_Base;
        template <class T>
        class Identifier;
        template <class T>
        class Const_Identifier;

        template <class T>
        bool operator == (const Identifier_Base< T >&, const Identifier_Base< T >&);
        template <class T>
        bool operator != (const Identifier_Base< T >&, const Identifier_Base< T >&);
        template <class T>
        std::ostream& operator << (std::ostream& os, const Identifier_Base< T >&);
        template <class T>
        std::ostream& operator << (std::ostream& os, const Factory< T >&);

        /** Pointer replacement class.
         *
         * Objects of type T can be referred to using an Identifier< T > as long as they
         * are stored by a Factory< T >.
         * @param T Type pointed at.
         */
        template <class T>
        class Identifier_Base
        {
        private:
            struct from_const_enabler {};
        public:
            typedef T val_type;
            typedef Factory< T > fact_type;

            /** Object resembling a NULL pointer in comparisons. */
            const Identifier_Base null() const { return Identifier_Base(); }
            /** Bool conversion. */
            operator bool() const { return _val != 0; }
            /** Prefix increment operator. */
            Identifier_Base& operator ++ () { ++_val; return *this; }
            /** Postfix increment operator. */
            Identifier_Base operator ++ (int) { Identifier_Base res(*this); ++(*this); return res; }

        protected:
            /** @name Constructors */
            /**@{*/
            /** Empty constructor. */
            Identifier_Base() : _val(0) {}
            /** Copy constructor. */
            Identifier_Base(const Identifier_Base& rhs) : _val(rhs._val) {}
            /**@}*/

            /** Dereference object. */
            val_type& dereference(const fact_type* fact_cptr = fact_type::get_active_ptr()) const
            {
                ASSERT(fact_cptr);
                return fact_cptr->get_elem(*this);
            }

        private:
            friend class Factory< T >;

            friend bool operator == <>(const Identifier_Base&, const Identifier_Base&);
            friend bool operator != <>(const Identifier_Base&, const Identifier_Base&);
            friend std::ostream& operator << <>(std::ostream&, const Identifier_Base&);

            uint32_t _val;
        };

        template <class T>
        bool operator == (const Identifier_Base< T >& lhs, const Identifier_Base< T >& rhs)
        {
            return lhs._val == rhs._val;
        }

        template <class T>
        bool operator != (const Identifier_Base< T >& lhs, const Identifier_Base< T >& rhs)
        {
            return !(lhs == rhs);
        }

        template <class T>
        std::ostream& operator <<(std::ostream& os, const Identifier_Base< T >& rhs)
        {
            os << "[_val=" << rhs._val << "]";
            return os;
        }

        /** Constant identifier. */
        template <class T>
        class Const_Identifier : public Identifier_Base< T >
        {
        public:
            Const_Identifier() {}
            Const_Identifier(const Const_Identifier& rhs) : Identifier_Base< T >(rhs) {}

            const T& operator * () const { return this->dereference(); }
            const T& operator -> () const { return this->dereference(); }
        };

        /** Non-constant identifier. */
        template <class T>
        class Identifier : public Identifier_Base< T >
        {
        public:
            Identifier() {}
            Identifier(const Identifier& rhs) : Identifier_Base< T >(rhs) {}
            explicit Identifier(const Const_Identifier< T >& rhs) : Identifier_Base< T >(rhs) {}
            operator const Const_Identifier< T >& () const { return *(Const_Identifier< T >*)this; }
            operator Const_Identifier< T >& () { return *(Const_Identifier< T >*)this; }

            T& operator * () const { return this->dereference(); }
            T& operator -> () const { return this->dereference(); }
        };

        /** Wrapper type used by Factory objects to store free node list without overhead. */
        template <class T>
        union Factory_Wrapper
        {
            typedef T val_type;
            typedef Identifier< T > idn_type;

            val_type _key;
            idn_type _next_free_idn;

            Factory_Wrapper() : _next_free_idn() {}
            Factory_Wrapper(const Factory_Wrapper& rhs)
                : _next_free_idn(rhs._next_free_idn) {}
        };

        /** Factory class that manages storage for objects of type T. */
        template <class T>
        class Factory
        {
        public:
            /** Type of object stored. */
            typedef T val_type;
            /** Non-constant pointer. */
            typedef Identifier< T > idn_type;
            /** Constant pointer. */
            typedef Const_Identifier< T > const_idn_type;
        private:
            typedef detail::Factory_Wrapper< T > wrapper_type;

        public:
            /** Default constructor.
             * @param activate Bool; if true, the object is used to resolve all managed pointers.
             */
            Factory(bool activate = false) : _cont(1, wrapper_type()), _next_free_idn() { if (activate) set_active(); }

        private:
            /** Dereference identifier. */
            wrapper_type& dereference(const Identifier_Base< T >& elem_idn_base) const
            {
                return const_cast<wrapper_type&>(_cont.at(elem_idn_base._val));
            }
            /** Get element pointed at. */
            val_type& get_elem(const Identifier_Base< T >& elem_idn_base) const
            {
                return dereference(elem_idn_base)._key;
            }

        public:
            /** Allocate space for new element and return pointer to it. */
            idn_type new_elem()
            {
                idn_type res;
                if (_next_free_idn)
                {
                    res = _next_free_idn;
                    _next_free_idn = _cont.at(_next_free_idn._val)._next_free_idn;
                }
                else
                {
                    _cont.push_back(wrapper_type());
                    res._val = _cont.size() - 1;
                }
                std::clog << "allocating space for new element at: "
                          << (void*)&_cont.at(res._val) << '\n';
                return res;
            }

            /** Delete element pointed at. */
            void del_elem(const_idn_type elem_idn)
            {
                _cont.at(elem_idn._val)._next_free_idn = _next_free_idn;
                _next_free_idn = idn_type(elem_idn);
            }

            /** Use this object to resolve all managed pointers. */
            void set_active() { _active_ptr = this; }

            /** Get active factory pointer. */
            static Factory* get_active_ptr() { return _active_ptr; }

        private:
            friend class Identifier_Base< T > ;

            static Factory* _active_ptr;

            friend std::ostream& operator << <>(std::ostream&, const Factory&);

            std::deque< wrapper_type > _cont;
            idn_type _next_free_idn;
        };

        template <class T>
        Factory< T >* Factory< T >::_active_ptr = NULL;

        template <class T>
        std::ostream& operator << (std::ostream& os, const Factory< T >& rhs)
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
        }
    }
    using detail::Identifier;
    using detail::Const_Identifier;
    using detail::Factory;
}


#endif
