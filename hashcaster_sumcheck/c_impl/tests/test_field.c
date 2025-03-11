#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <arm_neon.h>  // If you need NEON intrinsics
#include "unity.h"

#include "../field.h"

void test_f128_is_field(void)
{
    F128 a = f128_rand(); 
    F128 b = f128_rand();
    F128 c = f128_rand();

    F128 one = f128_one();

    // a * 1 == a
    TEST_ASSERT_TRUE(f128_eq(f128_mul(a, one), a));

    // a + b == b + a
    TEST_ASSERT_TRUE(f128_eq(f128_add(a, b), f128_add(b, a)));

    // a * b == b * a
    TEST_ASSERT_TRUE(f128_eq(f128_mul(a, b), f128_mul(b, a)));

    // (a + b) * c == a*c + b*c
    {
        F128 left  = f128_mul(f128_add(a, b), c);
        F128 right = f128_add(f128_mul(a, c), f128_mul(b, c));
        TEST_ASSERT_TRUE(f128_eq(left, right));
    }

    // (a * b) * c == a * (b * c)
    {
        F128 left  = f128_mul(f128_mul(a, b), c);
        F128 right = f128_mul(a, f128_mul(b, c));
        TEST_ASSERT_TRUE(f128_eq(left, right));
    }

    // fr(x) = x*x
    // check fr(a) + fr(b) == fr(a + b)
    {
        F128 fra  = f128_mul(a, a);     // a^2
        F128 frb  = f128_mul(b, b);     // b^2
        F128 sum  = f128_add(a, b);
        F128 frsum = f128_mul(sum, sum);// (a+b)^2

        F128 left  = f128_add(fra, frb);
        TEST_ASSERT_TRUE(f128_eq(left, frsum));
    }

    // Check a^(2^128) == a
    {
        F128 x = a;
        for (int i = 0; i < 128; i++) {
            // x = x^2
            x = f128_mul(x, x);
        }
        TEST_ASSERT_TRUE(f128_eq(a, x));
    }
}

/*-------------------------------------------------------
 * 2) Test Frobenius
 *
 * Rust logic:
 *   a.frob(i) == apow
 *   apow *= apow
 *   for i in [0..128)
 *------------------------------------------------------*/
void test_frobenius(void)
{
    F128 a = f128_rand();
    F128 apow = a;

    for (int i = 0; i < 128; i++) {
        F128 frob_i = f128_frob(a, i);   // a^(2^i)
        TEST_ASSERT_TRUE(f128_eq(frob_i, apow));
        // apow = apow^2
        apow = f128_mul(apow, apow);
    }
}

/*-------------------------------------------------------
 * 3) Test pi_as_expected
 *
 * Rust logic:
 *   r = random
 *   orbit[0] = r
 *   orbit[i+1] = orbit[i]^2
 *   for i in [0..128):
 *     bit = pi(i, orbit)
 *     check bit == 0 or 1
 *     compare to (r.raw >> i) & 1
 *------------------------------------------------------*/
void test_pi_as_expected(void)
{
    F128 r = f128_rand();

    /* Build the "orbit" array of size 128: each entry is squared from the previous. */
    F128 orbit[128];
    F128 tmp = r;
    for (int i = 0; i < 128; i++) {
        orbit[i] = tmp;
        tmp = f128_mul(tmp, tmp);
    }

    /* We need to check if pi(i, orbit) yields 0 or 1, 
       then compare to the i-th bit of r.raw */
    for (int i = 0; i < 128; i++) {
        F128 bit_val = pi_calc(i, orbit);

        int lhs; // will be 0 or 1
        if (f128_is_zero(bit_val)) {
            lhs = 0;
        } else if (f128_is_one(bit_val)) {
            lhs = 1;
        } else {
            TEST_FAIL_MESSAGE("pi(i, orbit) was neither 0 nor 1!");
            return;
        }

        // Compare to right-hand side: (r.raw >> i) & 1.
        // Suppose you have f128_get_low64(r) / f128_get_high64(r):
        uint64_t low  = r.low;
        uint64_t high = r.high;

        // If i < 64, bit is in low; else in high
        int shift = i;
        int rhs;
        if (shift < 64) {
            rhs = (int)((low >> shift) & 1ULL);
        } else {
            rhs = (int)((high >> (shift - 64)) & 1ULL);
        }

        TEST_ASSERT_EQUAL_INT(lhs, rhs);
    }
}

/*-------------------------------------------------------
 * 4) Test twists_logic_and
 *
 * Rust logic:
 *   a, b random
 *   build orbits a_orbit, b_orbit each of length 128
 *   answer = 0
 *   for i in [0..128):
 *     answer += basis(i)* pi(i,a_orbit)* pi(i,b_orbit)
 *   check answer == (a & b)
 *------------------------------------------------------*/
void test_twists_logic_and(void)
{
    F128 a = f128_rand();
    F128 b = f128_rand();

    // Build orbits:
    F128 a_orbit[128];
    F128 b_orbit[128];

    F128 _a = a;
    F128 _b = b;
    for (int i = 0; i < 128; i++) {
        a_orbit[i] = _a;
        b_orbit[i] = _b;
        // square them
        _a = f128_mul(_a, _a);
        _b = f128_mul(_b, _b);
    }

    // Compute answer = \sum_{i=0..127} basis(i)* pi(i,a_orbit)* pi(i,b_orbit)
    F128 answer = f128_zero();
    for (int i = 0; i < 128; i++) {
        F128 pia = pi_calc(i, a_orbit);
        F128 pib = pi_calc(i, b_orbit);
        F128 term = f128_mul(f128_mul(f128_basis(i), pia), pib);
        answer = f128_add(answer, term);
    }

    // expected = a & b
    F128 expected = f128_bitand(a, b);
    TEST_ASSERT_TRUE(f128_eq(answer, expected));
}
