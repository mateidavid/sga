#include <iostream>
#include "factory.hpp"
using namespace std;


#ifdef ERRORS
#define ERR(s) s
#else
#define ERR(s) (void)0
#endif

#ifndef FACTORY
struct A {};
struct B {};
typedef A A_holder;
typedef B B_holder;
typedef A* a_ptr_type;
typedef const A* ca_ptr_type;
typedef B* b_ptr_type;
typedef const B* cb_ptr_type;
#define EQ_0 = 0
#else
struct A {};
//struct B {};
typedef void B;
typedef factory::Holder< A > A_holder;
//typedef factory::Holder< B > B_holder;
typedef factory::Bounded_Pointer< A > a_ptr_type;
typedef factory::Bounded_Pointer< const A > ca_ptr_type;
typedef factory::Bounded_Pointer< B > b_ptr_type;
typedef factory::Bounded_Pointer< const B > cb_ptr_type;
factory::Factory< A > A_fact(true);
//factory::Factory< B > B_fact(true);
#define EQ_0
#endif

int main()
{
    A_holder a;
    //B_holder b;
    a_ptr_type a_ptr EQ_0;
    ca_ptr_type ca_ptr EQ_0;
    b_ptr_type b_ptr EQ_0;
    cb_ptr_type cb_ptr EQ_0;
    void* v_ptr = 0;
    const void *cv_ptr = 0;
    (void)a;
    //(void)b;
    (void)a_ptr;
    (void)ca_ptr;
    (void)b_ptr;
    (void)cb_ptr;
    (void)v_ptr;
    (void)cv_ptr;

    v_ptr = a_ptr; // A* -> void*
    cv_ptr = a_ptr; // A* -> const void*
    ca_ptr = a_ptr; // A* -> const A*
    ERR(b_ptr = a_ptr); // A* -> B*
    ERR(cb_ptr = a_ptr); // A* -> const B*

    ERR(v_ptr = ca_ptr); // const A* -> void*
    cv_ptr = ca_ptr; // const A* -> const void*
    ERR(a_ptr = ca_ptr); // const A* -> A*
    ERR(b_ptr = ca_ptr); // const A* -> B*
    ERR(cb_ptr = ca_ptr); // const A* -> const B*

    ERR(a_ptr = v_ptr); // void* -> A*
    ERR(ca_ptr = v_ptr); // void* -> const A*
    cv_ptr = v_ptr; // void* -> const void*

    ERR(a_ptr = cv_ptr); // const void* -> A*
    ERR(ca_ptr = cv_ptr); // const void* -> const A*
    ERR(v_ptr = cv_ptr); // const void* -> void*

    (void)*a_ptr;
    (void)*b_ptr;
    ERR((void)*cv_ptr);

    return 0;
}
