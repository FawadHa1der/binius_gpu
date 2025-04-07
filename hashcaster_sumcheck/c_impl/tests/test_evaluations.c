#include "test_evaluations.h"
#include "../mle_poly.h"
#include "../field.h"
#include "stdlib.h"
#include "string.h"
#include <time.h>    // for seeding
#include "stdio.h" // for uint8_t, uint64_t
#include "../debug_utils.h" // for f128_print
#include "../evaluations.h"
#include "unity.h"



// Test #1: test_pi
void test_pi(void)
{
    // The Rust code: random r, build orbit of 128 squarings
    // then for i in 0..128 => check pi(i)...

    // We'll do a partial approach:
    // generate a random r
    F128 r= f128_rand();

    // build orbit => data[0]= r; data[1]= r^2; ... data[127]= r^(2^127)
    Evaluations* orbit = malloc(sizeof(Evaluations));
    orbit->len= 128;
    orbit->elems= malloc(sizeof(F128)*128);
    {
        F128 tmp= r;
        for(int i=0;i<128;i++){
            orbit->elems[i]= tmp;
            // square tmp
            tmp= f128_mul(tmp,tmp);
        }
    }

    // Then check each i => bit= orbit.pi(i)
    // check if bit is 0 or 1 => we interpret if bit== ZERO =>0, if bit==ONE=>1
    // compare with (r.into_inner() >> i)%2
    // uint64_t r_inner= F128_into_inner(&r);
    for(int i=0;i<128;i++){
        F128 bit= pi_calc(i, orbit->elems);
        int lhs= -1;
        if(f128_eq(bit, f128_zero())) {
            lhs=0;
        } else if(f128_eq(bit, f128_one())) {
            lhs=1;
        } else {
            TEST_FAIL_MESSAGE("bit is neither 0 nor 1 in field!");
        }

        int rhs= f128_get_bit(r, i);
        TEST_ASSERT_EQUAL_MESSAGE(rhs, lhs, "Mismatch in test_pi at index i");
    }

    free(orbit->elems);
    free(orbit);
}


void test_pi_all_zeroes(void)
{
    // build an orbit of all zeros
    Evaluations* orbit;
    orbit= malloc(sizeof(Evaluations));
    orbit->len= 128;
    orbit->elems= malloc(sizeof(F128)*128);
    for(int i=0;i<128;i++){
        orbit->elems[i]= f128_zero();
    }

    // check orbit.pi(i) => 0
    for(int i=0;i<128;i++){
        F128 got= pi_calc(i, orbit->elems);
        TEST_ASSERT_TRUE_MESSAGE(f128_is_zero(got), "Failed for index i in test_pi_all_zeroes");
    }
    free(orbit->elems);
    free(orbit);
}

void test_pi_single_non_zero(void)
{
    // build an orbit of all zeros
    Evaluations* orbit;
    orbit= malloc(sizeof(Evaluations));
    orbit->len= 128;
    orbit->elems= malloc(sizeof(F128)*128);
    for(int i=0;i<128;i++){
        orbit->elems[i]= f128_zero();
    }
    // set orbit[5]= ONE
    orbit->elems[5]= f128_one();

    // compare with COBASIS_FROBENIUS_TRANSPOSE[i][5]
    for(int i=0;i<128;i++){
        F128 expected = COBASIS_FROBENIUS_TRANSPOSE[i][5];
        F128 got= pi_calc(i, orbit->elems);
        TEST_ASSERT_TRUE_MESSAGE(f128_eq(got,expected), "test_pi_single_non_zero mismatch");
    }
}

void test_pi_alternating(void)
{
    // orbit => data[i]= ONE if i%2==0 else ZERO
    Evaluations* orbit;
    orbit= malloc(sizeof(Evaluations));
    orbit->len= 128;
    orbit->elems= malloc(sizeof(F128)*128);

    for(int i=0;i<128;i++){
        if(i%2==0) orbit->elems[i]= f128_one();
        else       orbit->elems[i]= f128_zero();
    }

    for(int i=0;i<128;i++){
        // build "expected" as sum of basis's for j even
        // from snippet => we do the same logic:
        // expected= sum_{ j even} COBASIS_FROBENIUS_TRANSPOSE[i][j]
        F128 sum= f128_zero();
        for(int j=0;j<128;j++){
            if(j%2==0){
                F128 cval= COBASIS_FROBENIUS_TRANSPOSE[i][j];
                sum= f128_add(sum,cval);
            }
        }
        F128 got= pi_calc(i, orbit->elems);
        TEST_ASSERT_TRUE_MESSAGE(f128_eq(got,sum), "Mismatch in test_pi_alternating");
    }
}


void test_pi_random_orbit(void)
{
    // from snippet => we do a random orbit
    Evaluations* orbit;
    orbit= malloc(sizeof(Evaluations));
    orbit->len= 128;
    orbit->elems= malloc(sizeof(F128)*128);

    for(int i=0;i<128;i++){
        orbit->elems[i]= f128_rand();
    }

    // for i => expected= sum_{j=0..127} COBASIS_FROBENIUS_TRANSPOSE[i][j]* orbit[j]
    for(int i=0;i<128;i++){
        F128 sum= f128_zero();
        for(int j=0;j<128;j++){
            F128 cval= COBASIS_FROBENIUS_TRANSPOSE[i][j];
            // sum += cval* orbit[j]
            F128 tmp= f128_mul(cval, orbit->elems[j]);
            sum= f128_add(sum,tmp);
        }
        F128 got= pi_calc(i, orbit->elems);
        TEST_ASSERT_TRUE_MESSAGE(f128_eq(got,sum), "Mismatch in test_pi_random_orbit");
    }
}

void test_twist_untwist(void)
{
    // We'll define random "lhs"
    Evaluations* lhs;
    lhs= malloc(sizeof(Evaluations));
    lhs->len= 128;
    lhs->elems= malloc(sizeof(F128)*128);
    for(int i=0;i<128;i++){
        lhs->elems[i]= f128_rand();
    }
    // clone => "rhs"
    Evaluations* rhs= malloc(sizeof(Evaluations));
    rhs->len= 128;
    rhs->elems= malloc(sizeof(F128)*128);
    memcpy(rhs->elems, lhs->elems, sizeof(F128)*128);

    twist_evals(rhs);
    untwist_evals(rhs);
    // check equality
    for(int i=0;i<128;i++){
        TEST_ASSERT_TRUE_MESSAGE(f128_eq(lhs->elems[i], rhs->elems[i]), "twist->untwist did not restore original");
    }
    // do untwist->twist
    untwist_evals(rhs);
    twist_evals(rhs);
    for(int i=0;i<128;i++){
        TEST_ASSERT_TRUE_MESSAGE(f128_eq(lhs->elems[i], rhs->elems[i]), "untwist->twist did not restore original");
    }
    free(lhs->elems);
    free(lhs);
    free(rhs->elems);
    free(rhs);
}

void test_twist_all_zeros(void)
{
    Evaluations* evals;
    evals= malloc(sizeof(Evaluations));
    evals->len= 128;
    evals->elems= malloc(sizeof(F128)*128);
    // set all to zero
    for(int i=0;i<128;i++){
        evals->elems[i]= f128_zero();
    }
    twist_evals(evals);

    // check all still zero
    for(int i=0;i<128;i++){
        TEST_ASSERT_TRUE_MESSAGE(f128_is_zero(evals->elems[i]), "twist all zero input failed");
    }
}

void test_untwist_all_zeros(void)
{
    Evaluations* evals;
    evals= malloc(sizeof(Evaluations));
    evals->len= 128;
    evals->elems= malloc(sizeof(F128)*128);
    // set all to zero
    for(int i=0;i<128;i++){
        evals->elems[i]= f128_zero();
    }
    untwist_evals(evals);

    for(int i=0;i<128;i++){
        TEST_ASSERT_TRUE_MESSAGE(f128_is_zero(evals->elems[i]), "untwist all zero input failed");
    }
}


