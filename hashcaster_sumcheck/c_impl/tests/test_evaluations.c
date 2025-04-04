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


