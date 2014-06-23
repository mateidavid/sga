#include "shortcuts.hpp"
#include <iostream>
#include <type_traits>
#include <boost/intrusive/pointer_traits.hpp>
#include <boost/tti/tti.hpp>
#include "factory.hpp"

using namespace std;
namespace bi = boost::intrusive;

struct A {};
struct B {};

typedef const A const_A;

typedef bounded::Pointer< A > a_ptr_t;
typedef bounded::Pointer< const A > a_const_ptr_t;
typedef bounded::Reference< A > a_ref_t;
typedef bounded::Reference< const A > a_const_ref_t;

typedef bi::pointer_traits< a_ptr_t > a_ptr_tr;
typedef bi::pointer_traits< a_const_ptr_t > a_const_ptr_tr;

static_assert(is_same< bounded::Pointer< const_A >, bounded::Pointer< const const_A > >::value, "ptr< const_A > != ptr< const const_A >");

static_assert(is_same< typename a_ptr_tr::element_type, A >::value, "a_ptr_tr::element_type != A");
static_assert(is_same< typename a_ptr_tr::pointer, a_ptr_t >::value, "a_ptr_tr::pointer != a_ptr_t");
//static_assert(is_same< typename a_ptr_tr::const_pointer, a_const_ptr_t >::value, "a_ptr_tr::const_pointer != a_const_ptr_t");
static_assert(is_same< typename a_ptr_tr::reference, a_ref_t >::value, "a_ptr_tr::reference != a_ref_t");

static_assert(is_same< boost::intrusive::detail::type_rebinder< bounded::Pointer< A >, const A >::type, bounded::Pointer< const A > >::value, "rebind error");

static_assert(is_same< typename a_ptr_tr::template rebind_pointer< A >::type, a_ptr_t >::value, "a_ptr_tr::rebind_pointer< A > != a_ptr_t");
//static_assert(is_same< typename a_ptr_tr::template rebind_pointer< void >::type, void* >::value, "a_ptr_tr::rebind_pointer< void > != void*");
//static_assert(is_same< typename a_ptr_tr::template rebind_pointer< B >::type, B* >::value, "a_ptr_tr::rebind_pointer< B > != B*");
static_assert(is_same< typename a_ptr_tr::template rebind_pointer< const A >::type, a_const_ptr_t >::value, "a_ptr_tr::rebind_pointer< const A > != a_const_ptr_t");

static_assert(is_same< typename a_const_ptr_tr::pointer, a_const_ptr_t >::value, "a_const_ptr_tr::pointer != a_const_ptr_t");
//static_assert(is_same< typename a_const_ptr_tr::const_pointer, a_const_ptr_t >::value, "a_const_ptr_tr::const_pointer != a_const_ptr_t");
static_assert(is_same< typename a_const_ptr_tr::template rebind_pointer< A >::type, a_ptr_t >::value, "a_const_ptr_tr::rebind_pointer< A > != a_ptr_t");

int main()
{}
