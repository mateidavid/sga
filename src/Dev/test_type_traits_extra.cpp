#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <vector>
#include <list>
#include <functional>
#include <boost/preprocessor.hpp>

#include "type_traits_extra.hpp"

using namespace std_extra;

struct A
{
    typedef int type_1;
    typedef int type_2;
    void fcn_1_AB_v();
    int fcn_1_AB_i();
    void fcn_2_A_v();
    int fcn_2_A_i();
    int fcn_2_Ac_i() const;
    int fcn_4_A_Bc_i();
    int fcn_5_Ac_B_i() const;
    void fcn_6_A_v_i(int);
    int& fcn_7_A(int);
    float& fcn_7_A(float);
};

struct B : public A
{
    typedef void* type_1;
    typedef void* type_3;
    void fcn_1_AB_v();
    int fcn_1_AB_i();
    void fcn_3_B_v();
    int fcn_3_B_i();
    int fcn_3_Bc_i() const;
    int fcn_4_A_Bc_i() const;
    int fcn_5_Ac_B_i();
};

LOOSE_HAS_MEM_TYPE(type_1)
LOOSE_HAS_MEM_TYPE(type_2)
LOOSE_HAS_MEM_TYPE(type_3)

LOOSE_HAS_MEM_FUN(fcn_1_AB_v)
LOOSE_HAS_MEM_FUN(fcn_1_AB_i)
LOOSE_HAS_MEM_FUN(fcn_2_A_v)
LOOSE_HAS_MEM_FUN(fcn_2_A_i)
LOOSE_HAS_MEM_FUN(fcn_3_B_v)
LOOSE_HAS_MEM_FUN(fcn_3_B_i)
LOOSE_HAS_MEM_FUN(fcn_2_Ac_i)
LOOSE_HAS_MEM_FUN(fcn_3_Bc_i)
LOOSE_HAS_MEM_FUN(fcn_4_A_Bc_i)
LOOSE_HAS_MEM_FUN(fcn_5_Ac_B_i)
LOOSE_HAS_MEM_FUN(fcn_6_A_v_i)
LOOSE_HAS_MEM_FUN(fcn_7_A)

TEST_CASE("has_mem_type")
{
    SECTION("inheritance")
    {
        // A: has type_1 & type_2, but not type_3
        CHECK(( is_same< has_mem_type_type_1< A >::type, true_type >::value ));
        CHECK(( is_same< has_mem_type_type_2< A >::type, true_type >::value ));
        CHECK(( is_same< has_mem_type_type_3< A >::type, false_type >::value ));
        // B: has all 3 types
        CHECK(( is_same< has_mem_type_type_1< B >::type, true_type >::value ));
        CHECK(( is_same< has_mem_type_type_2< B >::type, true_type >::value ));
        CHECK(( is_same< has_mem_type_type_3< B >::type, true_type >::value ));
    }
    SECTION("constness")
    {
        // same types
        CHECK(( is_same< has_mem_type_type_2< const A >::type, true_type >::value ));
        CHECK(( is_same< has_mem_type_type_3< const A >::type, false_type >::value ));
        CHECK(( is_same< has_mem_type_type_2< const B >::type, true_type >::value ));
        CHECK(( is_same< has_mem_type_type_3< const B >::type, true_type >::value ));
    }
    SECTION("lvalue_reference")
    {
        // same types
        CHECK(( is_same< has_mem_type_type_2< A& >::type, true_type >::value ));
        CHECK(( is_same< has_mem_type_type_3< A& >::type, false_type >::value ));
        CHECK(( is_same< has_mem_type_type_2< B& >::type, true_type >::value ));
        CHECK(( is_same< has_mem_type_type_3< B& >::type, true_type >::value ));
    }
    SECTION("rvalue_reference")
    {
        // same types
        CHECK(( is_same< has_mem_type_type_2< A&& >::type, true_type >::value ));
        CHECK(( is_same< has_mem_type_type_3< A&& >::type, false_type >::value ));
        CHECK(( is_same< has_mem_type_type_2< B&& >::type, true_type >::value ));
        CHECK(( is_same< has_mem_type_type_3< B&& >::type, true_type >::value ));
    }
    SECTION("other_arguments")
    {
        CHECK(( is_same< has_mem_type_type_1< void >::type, false_type >::value ));
        CHECK(( is_same< has_mem_type_type_1< int >::type, false_type >::value ));
        CHECK(( is_same< has_mem_type_type_1< int* >::type, false_type >::value ));
        CHECK(( is_same< has_mem_type_type_1< A* >::type, false_type >::value ));
        CHECK(( is_same< has_mem_type_type_1< std::function< void(void)> >::type, false_type >::value ));
    }
}

TEST_CASE("has_mem_fun")
{
    SECTION("inheritance")
    {
        CHECK(( is_same< has_mem_fun_fcn_1_AB_v< A, void >::type, true_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_2_A_v< A, void >::type,  true_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_3_B_v< A, void >::type,  false_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_1_AB_v< B, void >::type, true_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_2_A_v< B, void >::type,  true_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_3_B_v< B, void >::type,  true_type >::value ));
    }
    SECTION("constness")
    {
        CHECK(( is_same< has_mem_fun_fcn_1_AB_i< A, int >::type,       true_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_1_AB_i< const A, int >::type, false_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_2_Ac_i< A, int >::type,       true_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_2_Ac_i< const A, int >::type, true_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_3_Bc_i< const A, int >::type, false_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_3_Bc_i< const B, int >::type, true_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_3_Bc_i< const A, int >::type, false_type >::value ));
        //
        CHECK(( is_same< has_mem_fun_fcn_4_A_Bc_i< A, int >::type,       true_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_4_A_Bc_i< const A, int >::type, false_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_4_A_Bc_i< B, int >::type,       true_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_4_A_Bc_i< const B, int >::type, true_type >::value ));
        //
        CHECK(( is_same< has_mem_fun_fcn_5_Ac_B_i< A, int >::type,       true_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_5_Ac_B_i< const A, int >::type, true_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_5_Ac_B_i< B, int >::type,       true_type >::value ));
         // non-const version hides inherited version
        CHECK(( is_same< has_mem_fun_fcn_5_Ac_B_i< const B, int >::type, false_type >::value ));
    }
    SECTION("reference")
    {
        CHECK(( is_same< has_mem_fun_fcn_2_A_v< A&, void >::type,        true_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_2_A_v< A&&, void >::type,       true_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_2_A_v< const A&, void >::type,  false_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_2_A_v< const A&&, void >::type, false_type >::value ));
    }
    SECTION("return_value_conversion")
    {
        CHECK(( is_same< has_mem_fun_fcn_2_A_i< A, char >::type, true_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_2_A_i< A, long >::type, true_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_2_A_v< A, int >::type, false_type >::value ));
    }
    SECTION("argument_conversion")
    {
        CHECK(( is_same< has_mem_fun_fcn_6_A_v_i< A, void, int >::type,  true_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_6_A_v_i< A, void, char >::type, true_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_6_A_v_i< A, void, long >::type, true_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_6_A_v_i< A, void, void >::type, false_type >::value ));
    }
    SECTION("overload_selection")
    {
        CHECK(( is_same< has_mem_fun_fcn_7_A< A, int&, int >::type,     true_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_7_A< A, float&, float >::type, true_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_7_A< A, int&, float >::type,   false_type >::value ));
        CHECK(( is_same< has_mem_fun_fcn_7_A< A, float&, int >::type,   false_type >::value ));
    }
}
