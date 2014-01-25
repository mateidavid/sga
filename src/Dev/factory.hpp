#ifndef __FACTORY_HPP
#define __FACTORY_HPP

#include <iostream>
#include <boost/mpl/if.hpp>
#include <boost/utility/enable_if.hpp>

#include "common.hpp"


namespace factory
{
    namespace detail
    {
        template <class T>
        class Factory;

        /** Pointer replacement class.
         *
         * Objects of type T can be referred to using an Identifier< T > as long as they
         * are stored by a Factory< T >.
         * @param T Type pointed at.
         * @param const_val Bool; true iff object points to a constant value.
         */
        template <class T, bool const_val>
        class Qual_Identifier
        {
        private:
            struct from_const_enabler {};
        public:
            typedef T val_type;
            typedef typename boost::mpl::if_c< const_val,
                                               const val_type,
                                               val_type >::type qual_val_type;
            typedef Factory< T > fact_type;

            /** @name Constructors */
            /**@{*/
            /** Empty constructor. */
            Qual_Identifier() : _val(0) {}
            /** Constructor from another non-constant identifier; always enabled. */
            Qual_Identifier(const Qual_Identifier< T, false >& rhs) : _val(rhs._val) {}
            /** Constructor from another constant identifier; enabled only for constant identifiers. */
            template <bool const_val_other>
            Qual_Identifier(const Qual_Identifier< T, const_val_other >& rhs,
                            typename boost::enable_if_c<
                                const_val_other && const_val,
                                from_const_enabler
                                >::type = from_const_enabler()) : _val(rhs._val) {}
            /**@}*/

            /** Object resembling a NULL pointer in comparisons. */
            const Qual_Identifier null() const { return Qual_Identifier(); }
            /** Main dereference operator; uses active factory pointer to retreive object. */
            qual_val_type& operator * () const
            {
                ASSERT(fact_type::_active_ptr);
                return fact_type::_active_ptr->get_elem(*this);
            }
            /** Secondary dereference operator; NOTE: target type must implement it as well. */
            qual_val_type& operator -> () const { return *(*this); }
            /** Bool conversion. */
            operator bool() const { return _val != 0; }

        private:
            friend class Qual_Identifier< T, not const_val >;

            friend class Factory< T >;
            template <class U, bool const_val_1, bool const_val_2>
            friend bool operator == (const Qual_Identifier<U, const_val_1>&,
                                     const Qual_Identifier<U, const_val_2>&);
            template <class U, bool const_val_1, bool const_val_2>
            friend bool operator != (const Qual_Identifier<U, const_val_1>&,
                                     const Qual_Identifier<U, const_val_2>&);
            template <class U, bool const_val_other>
            friend std::ostream& operator << (std::ostream&,
                                              const Qual_Identifier<U, const_val_other>&);
            template <class U>
            friend std::ostream& operator << (std::ostream& os, const Factory< U >& rhs);

            uint32_t _val;
        };

        template <class T, bool const_1, bool const_2>
        bool operator == (const Qual_Identifier< T, const_1 >& lhs,
                          const Qual_Identifier< T, const_2 >& rhs)
        {
            return lhs._val == rhs._val;
        }

        template <class T, bool const_1, bool const_2>
        bool operator != (const Qual_Identifier< T, const_1 >& lhs,
                          const Qual_Identifier< T, const_2 >& rhs)
        {
            return !(lhs == rhs);
        }

        template <class T, bool const_val>
        std::ostream& operator <<(std::ostream& os,
                                  const Qual_Identifier< T, const_val >& rhs)
        {
            os << "[_val=" << rhs._val << "]";
            return os;
        }

        /** Non-constant identifier. */
        template <class T>
        class Identifier : public Qual_Identifier< T, false >
        {
        public:
            Identifier() {}
            Identifier(const Qual_Identifier< T, false >& rhs)
                : Qual_Identifier< T, false >(rhs) {}
        };

        /** Constant identifier. */
        template <class T>
        class Const_Identifier : public Qual_Identifier< T, true >
        {
        public:
            Const_Identifier() {}
            Const_Identifier(const Qual_Identifier< T, true >& rhs)
                : Qual_Identifier< T, true >(rhs) {}
            Const_Identifier(const Qual_Identifier< T, false >& rhs)
                : Qual_Identifier< T, true >(rhs) {}
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

            /** Get element pointed at. */
            val_type& get_elem(idn_type elem_idn)
            {
                return _cont.at(elem_idn._val)._key;
            }

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
            void del_elem(idn_type elem_idn)
            {
                _cont.at(elem_idn._val)._next_free_idn = _next_free_idn;
                _next_free_idn = elem_idn;
            }

            /** Use this object to resolve all managed pointers. */
            void set_active() { _active_ptr = this; }

        private:
            friend class Qual_Identifier< T, true >;
            friend class Qual_Identifier< T, false >;

            static Factory* _active_ptr;

            template <class U>
            friend std::ostream& operator << (std::ostream&, const Factory< U >&);

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
                 crt = rhs._cont.at(crt._val)._next_free_idn)
            {
                ++n_free;
            }
            os << "n_free=" << n_free << '\n'
               << "next_free=" << rhs._next_free_idn._val << '\n';
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
