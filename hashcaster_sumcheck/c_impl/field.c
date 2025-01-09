/*****************************************************************************
 * f128.c
 *****************************************************************************/
#include "field.h"
#include <stdlib.h>
#include <time.h>
#include <arm_neon.h>
/* 
 * A small struct to hold the three vectors (H, M, L) that the
 * Rust version returns as a tuple (h, m, l).
 */
typedef struct {
    uint8x16_t h;
    uint8x16_t m;
    uint8x16_t l;
} karatsuba1_result;


/* We will return two 128-bit vectors (x23 and x01) */
typedef struct {
    uint8x16_t x23;
    uint8x16_t x01;
} karatsuba2_result;

/*
 * Example placeholders for global tables from your Rust code:
 *   extern F128 FROBENIUS[128][128];
 *   extern F128 COBASIS[128];
 *   extern F128 COBASIS_FROBENIUS[128][128];
 *
 * In real usage, you’d either:
 *   1) define them here with their actual values, or
 *   2) load them from a separate .c file.
 */

// For illustration only — not real data:
F128 FROBENIUS[128][128];
F128 COBASIS[128];
F128 COBASIS_FROBENIUS[128][128];

/**
 * A placeholder for the GF(2^128) polynomial multiplication function.
 * Replace with your actual field-multiplication logic, or link from
 * another source.
 */
F128 mul_128(const F128 *a, const F128 *b)
{
    // You must implement polynomial multiplication mod your chosen
    // irreducible polynomial. This is just a placeholder that returns 0.
    (void)a;
    (void)b;
    return f128_zero();
}

/**
 * Example random function returning a pseudo-random F128.
 * In a real implementation, use a secure RNG or wrap `rand()` carefully,
 * and seed it properly, etc.
 */
F128 f128_rand(void)
{
    F128 out;
    out.lo = ((uint64_t)rand() << 32) ^ (uint64_t)rand();
    out.hi = ((uint64_t)rand() << 32) ^ (uint64_t)rand();
    return out;
}

/*
 * One-time initialization (optional). You might seed the RNG here or
 * fill the FROBENIUS, COBASIS, COBASIS_FROBENIUS tables, etc.
 */
void f128_init(void)
{
    srand((unsigned int)time(NULL));

    // Example: fill FROBENIUS, COBASIS, etc. with dummy data if needed.
    // In real code, load or define the actual precomputed tables.
    for (int i = 0; i < 128; i++) {
        // COBASIS[i] = basis(i), for instance
        COBASIS[i] = f128_basis(i);
        for (int j = 0; j < 128; j++) {
            // e.g., FROBENIUS[i][j] = f128_zero();
            // or some real table data
            FROBENIUS[i][j] = f128_zero();
            COBASIS_FROBENIUS[i][j] = f128_zero();
        }
    }
}

/*
 * pmull:
 *   Multiplies the LOW 64 bits of a and b (lane 0).
 */
static inline uint8x16_t pmull(uint8x16_t a, uint8x16_t b)
{
    // Reinterpret 'a' and 'b' as two 64-bit lanes each
    uint64x2_t aa = vreinterpretq_u64_u8(a);
    uint64x2_t bb = vreinterpretq_u64_u8(b);

    // Extract the low 64-bit lane (index = 0)
    uint64_t a0 = vgetq_lane_u64(aa, 0);
    uint64_t b0 = vgetq_lane_u64(bb, 0);

    // Cast to poly64_t for polynomial multiplication
    poly64_t ap = a0;
    poly64_t bp = b0;

    // PMULL - multiply the two 64-bit polynomials, result is 128-bit
    poly128_t result = vmull_p64(ap, bp);

    // Convert back to a 128-bit vector of u8
    return (uint8x16_t)result;
}

/*
 * pmull2:
 *   Multiplies the HIGH 64 bits of a and b (lane 1).
 */
static inline uint8x16_t pmull2(uint8x16_t a, uint8x16_t b)
{
    // Reinterpret 'a' and 'b' as two 64-bit lanes each
    uint64x2_t aa = vreinterpretq_u64_u8(a);
    uint64x2_t bb = vreinterpretq_u64_u8(b);

    // Extract the high 64-bit lane (index = 1)
    uint64_t a1 = vgetq_lane_u64(aa, 1);
    uint64_t b1 = vgetq_lane_u64(bb, 1);

    // Cast to poly64_t for polynomial multiplication
    poly64_t ap = a1;
    poly64_t bp = b1;

    // PMULL - multiply the two 64-bit polynomials, result is 128-bit
    poly128_t result = vmull_p64(ap, bp);

    // Convert back to a 128-bit vector of u8
    return (uint8x16_t)(result);
}


static inline karatsuba1_result karatsuba1(uint8x16_t x, uint8x16_t y)
{
    /* 
     * "hi" is obtained by shifting x (or y) by 8 bytes 
     * (the second parameter to vextq_u8 is just 'x' again for a shift).
     * So x_hi corresponds to the top 64 bits of x, x_lo is the bottom 64 bits.
     * Similarly for y_hi.
     */
    uint8x16_t x_hi = vextq_u8(x, x, 8);  // shift by 8 bytes
    uint8x16_t y_hi = vextq_u8(y, y, 8);

    uint8x16_t x_xor = veorq_u8(x, x_hi); // x.hi ^ x.lo
    uint8x16_t y_xor = veorq_u8(y, y_hi); // y.hi ^ y.lo

    /*
     * M = (x.hi ^ x.lo) * (y.hi ^ y.lo)
     * H = x.hi * y.hi
     * L = x.lo * y.lo
     */
    uint8x16_t m = pmull(x_xor, y_xor);
    uint8x16_t h = pmull2(x, y);
    uint8x16_t l = pmull(x, y);

    karatsuba1_result out;
    out.h = h;
    out.m = m;
    out.l = l;
    return out;
}



/*
 * Karatsuba combine step.
 *
 * Given:
 *   h, m, l (each a 128-bit vector of type uint8x16_t)
 *
 * This function outputs two 128-bit vectors (x23, x01) that combine
 * partial products into a 2n-bit result in GF(2^n).
 */
static inline karatsuba2_result karatsuba2(uint8x16_t h, uint8x16_t m, uint8x16_t l)
{
    /*
     * t0 = m ^ vextq_u8(l, h, 8)
     *   = {m0, m1} ^ {l1, h0}
     *   = {m0^l1, m1^h0}
     */
    uint8x16_t t0 = veorq_u8(m, vextq_u8(l, h, 8));

    /*
     * t1 = h ^ l
     *   = {h0^l0, h1^l1}
     */
    uint8x16_t t1 = veorq_u8(h, l);

    /*
     * t = t0 ^ t1
     *   = {m0^l1^h0^l0, m1^h0^h1^l1}
     */
    uint8x16_t t = veorq_u8(t0, t1);

    /*
     * x01 = vextq_u8( vextq_u8(l, l, 8), t, 8 )
     *   = combine {l1, l0} with t, then shift 8 bytes
     *   = {m0^l1^h0^l0, l0}
     */
    uint8x16_t x01 = vextq_u8(vextq_u8(l, l, 8), t, 8);

    /*
     * x23 = vextq_u8(t, vextq_u8(h, h, 8), 8)
     *   = combine t with {h1, h0}, then shift 8 bytes
     *   = {h1, m1^h0^h1^l1}
     */
    uint8x16_t x23 = vextq_u8(t, vextq_u8(h, h, 8), 8);

    /* Return them in a struct */
    karatsuba2_result out;
    out.x23 = x23;
    out.x01 = x01;
    return out;
}


/*---------------------------
 * mul_128:
 *   1) Convert two 128-bit inputs (x, y) from (hi, lo) to NEON vectors
 *   2) Perform karatsuba1 -> (h, m, l)
 *   3) Perform karatsuba2 -> combine to (h, l)
 *   4) Montgomery reduce
 *   5) Convert final vector back to (uint128_t)
 *---------------------------*/
static inline uint128_t mul_128(uint128_t x, uint128_t y)
{
    /*
     * 1) Convert x, y into NEON vectors
     */
    u128_as_vec X, Y, OUT;
    X.lo = x.lo;
    X.hi = x.hi;
    Y.lo = y.lo;
    Y.hi = y.hi;

    /*
     * 2) karatsuba1 -> (h, m, l)
     */
    karatsuba1_result km = karatsuba1(X.vec, Y.vec);

    /*
     * 3) karatsuba2 -> (h, l) from partial products
     *    The Rust code says: let (h, l) = karatsuba2(h, m, l);
     *    We'll name them x23 and x01, consistent with our karatsuba2_result.
     */
    karatsuba2_result kl = karatsuba2(km.h, km.m, km.l);

    /*
     * 4) mont_reduce( x23, x01 ) -> final 128-bit vector
     */
    uint8x16_t reduced = mont_reduce(kl.x23, kl.x01);

    /*
     * 5) Convert the final 128-bit vector back to a struct {lo, hi}.
     */
    OUT.vec = reduced;
    uint128_t result;
    result.lo = OUT.lo;
    result.hi = OUT.hi;

    return result;
}
