/*****************************************************************************
 * f128.h
 *****************************************************************************/
#ifndef F128_H
#define F128_H

#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <assert.h>
#include <arm_neon.h>
#include "types.h"
#include "cobasis_table.h"
#include "frobenius_table.h"
#include "cobasis_frobenius_table.h"

/**
 * Convert a 128-bit struct into raw 128 bits (e.g., two uint64_ts).
 *
 * For convenience, we’ll represent the raw 128 bits in a pair of 64-bit values.
 * You could also return them via pointer or union, depending on your use-case.
 */
static inline void f128_into_raw(const F128 *x, uint64_t *lo, uint64_t *hi)
{
    *lo = x->low;
    *hi = x->high;
}

/**
 * Create an F128 from two 64-bit words: (hi, lo).
 */
static inline F128 f128_from_raw(uint64_t lo, uint64_t hi)
{
    F128 out;
    out.low = lo;
    out.high = hi;
    return out;
}

/**
 * f128_eq - Checks if two F128 elements are equal.
 *
 * Return:
 *   true if they are identical (both low and high 64-bit halves match),
 *   false otherwise.
 */
static inline bool f128_eq(F128 a, F128 b)
{
    return (a.low == b.low) && (a.high == b.high);
}

/**
 * Return an F128 that is 0 in GF(2^128).
 */
static inline F128 f128_zero(void)
{
    return f128_from_raw(0ULL, 0ULL);
}

static inline F128 f128_one(void)
{
    return f128_from_raw(0x0000000000000001ULL, 0x0000000000000000ULL);
}

/**
 * Checks if F128 is zero in GF(2^128).
 */
static inline bool f128_is_zero(const F128 x)
{
    return (x.low == 0ULL && x.high == 0ULL);
}

/**
 * Checks if F128 is one in GF(2^128).
 */
static inline bool f128_is_one(const F128 x)
{
    return (x.low == 1ULL && x.high == 0ULL);
}


/**
 * GF(2^128) “addition” is just bitwise XOR.
 */
static inline F128 f128_add(const F128 a, const F128 b)
{
    F128 out;
    out.low = a.low ^ b.low;
    out.high = a.high ^ b.high;
    return out;
}

/**
 * In-place GF(2^128) “+=” (XOR).
 */
static inline void f128_add_assign(F128 *a, const F128 *b)
{
    a->low ^= b->low;
    a->high ^= b->high;
}

/**
 * GF(2^128) bitwise AND.
 */
static inline F128 f128_bitand(const F128 a, const F128 b)
{
    F128 out;
    out.low = a.low & b.low;
    out.high = a.high & b.high;
    return out;
}

/**
 * In-place AND.
 */
static inline void f128_bitand_assign(F128 *a, const F128 *b)
{
    a->low &= b->low;
    a->high &= b->high;
}

/**
 * Example placeholder for random generation.
 * You must implement this to produce a random 128-bit value.
 */
extern F128 f128_rand(void);

/**
 * Multiply two GF(2^128) tower elements.  
 */
extern F128 mul_128(const F128 a, const F128 b);

/**
 * F128 * F128 => F128
 */
static inline F128 f128_mul(const F128 a, const F128 b)
{
    return mul_128(a, b);
}

/**
 * In-place multiplication: a *= b
 */
static inline void f128_mul_assign(F128 *a, const F128 *b)
{
    *a = mul_128(*a, *b);
}

extern const F128 FROBENIUS[128][128]; // Defined in frobenius_table.h

/**
 * Convert an F128 into a bool array of length 128, little-endian style:
 * bits[0] = least significant bit of (a->low)
 * bits[127] = most significant bit of (a->high)
 */
static inline void f128_to_bits(const F128 a, bool out_bits[128])
{
    // Fill from a->low:
    for (int i = 0; i < 64; i++) {
        out_bits[i] = (a.low >> i) & 1ULL;
    }
    // Fill from a.high:
    for (int i = 0; i < 64; i++) {
        out_bits[64 + i] = (a.high >> i) & 1ULL;
    }
}

/**
 * Frobenius operator: x -> x^(2^k), with k possibly negative, etc.
 * Mirrors the Rust logic with a shift of k modulo 128.
 */
static inline F128 f128_frob(const F128 x, int k)
{
    // Adjust negative k in the same way the Rust code does:
    if (k < 0) {
        k = -k;
        k %= 128;
        k = -k;
        k += 128;
    } else {
        k %= 128;
    }

    // Retrieve the 2D matrix from your `FROBENIUS`.
    // We'll use the row for k (FROBENIUS[k]) to mix bits.
    bool bits[128];
    f128_to_bits(x, bits);

    // XOR of the appropriate rows:
    F128 result = f128_zero();
    for (int i = 0; i < 128; i++) {
        if (bits[i]) {
            // XOR with FROBENIUS[k][i]
            F128 tmp = f128_add(result, FROBENIUS[k][i]);
            result = tmp;
        }
    }
    return result;
}

/**
 * Return the basis element e_i = 1 << i in your chosen representation.
 * Because we store (lo, hi), if i < 64, it goes in lo; if i >= 64, in hi.
 */
static inline F128 f128_basis(int i)
{
    assert(i >= 0 && i < 128);
    if (i < 64) {
        return f128_from_raw((uint64_t)1 << i, 0ULL);
    } else {
        return f128_from_raw(0ULL, (uint64_t)1 << (i - 64));
    }
}

/**
 * COBASIS table.  In Rust: pub const COBASIS: [u128; 128] = ...
 * In C, we might keep it as array of (lo, hi) pairs or F128 elements.
 */
extern const F128 COBASIS[128]; // Defined in cobasis_table.h

/**
 * Return COBASIS[i].
 */
static inline F128 f128_cobasis(int i)
{
    return COBASIS[i];
}

/**
 * In Rust:  pub fn pi(i: usize, twists: &[F128]) -> F128 { ... }
 * We do the same in C, but using pointers.
 * This computes \sum_j (COBASIS_FROBENIUS[j][i]^(2^j)) * twists[j], or
 * more exactly the version in your Rust code:
 *   \sum_j COBASIS_FROBENIUS[j][i] * twists[j].
 *
 * But from your code snippet:
 *   ret += F128::from_raw(COBASIS_FROBENIUS[j][i]) * twists[j];
 * so we just replicate that.
 */
extern const F128 COBASIS_FROBENIUS[128][128]; // Defined in cobasis_frobenius_table.h

static inline F128 pi_calc(int i, const F128 *twists /* array of size 128 */)
{
    assert(i >= 0 && i < 128);
    F128 ret = f128_zero();
    for (int j = 0; j < 128; j++) {
        // ret += F128::from_raw(COBASIS_FROBENIUS[j][i]) * twists[j];
        F128 tmpmul = f128_mul(COBASIS_FROBENIUS[j][i], twists[j]);
        ret = f128_add(ret, tmpmul);
    }
    return ret;
}

#endif // F128_H