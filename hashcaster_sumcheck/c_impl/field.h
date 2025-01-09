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
// For a GF(2^128) element, we store two 64-bit words:
//   hi = the high 64 bits
//   lo = the low 64 bits
typedef struct {
    uint64_t lo;
    uint64_t hi;
} F128;



/*---------------------------
 * A 128-bit type in C
 *---------------------------*/
typedef struct {
    uint64_t lo;  // Low 64 bits
    uint64_t hi;  // High 64 bits
} uint128_t;

/*---------------------------
 * Helper union to reinterpret
 * between (uint128_t) and (uint8x16_t).
 *---------------------------*/
typedef union {
    uint8x16_t vec;   // 128-bit NEON vector
    uint64x2_t u64x2; // Also 128 bits, but two 64-bit lanes
    struct {
        uint64_t lo;
        uint64_t hi;
    };
} u128_as_vec;


/**
 * Convert a 128-bit struct into raw 128 bits (e.g., two uint64_ts).
 *
 * For convenience, we’ll represent the raw 128 bits in a pair of 64-bit values.
 * You could also return them via pointer or union, depending on your use-case.
 */
static inline void f128_into_raw(const F128 *x, uint64_t *lo, uint64_t *hi)
{
    *lo = x->lo;
    *hi = x->hi;
}

/**
 * Create an F128 from two 64-bit words: (hi, lo).
 */
static inline F128 f128_from_raw(uint64_t lo, uint64_t hi)
{
    F128 out;
    out.lo = lo;
    out.hi = hi;
    return out;
}

/**
 * Return an F128 that is 0 in GF(2^128).
 */
static inline F128 f128_zero(void)
{
    return f128_from_raw(0ULL, 0ULL);
}

/**
 * Return an F128 that is 1 in GF(2^128).
 * In the Rust code, the decimal is: 257870231182273679343338569694386847745
 * You can represent this constant in hex, or split it into two 64-bit words.
 *
 *  257870231182273679343338569694386847745 in hex is
 *  0xC9... but let's do the full split. 
 *
 *  One way to get this is:
 *    257870231182273679343338569694386847745 in hex = 
 *       0xC9BB_350D_B6AB_FD4F_8A5A_C6FF_FB_F001
 *    hi = 0xC9BB350DB6ABFD4F
 *    lo = 0x8A5AC6FFFBF001
 *
 *  The exact irreducible polynomial or representation might differ
 *  in your codebase. Adapt accordingly.
 */
static inline F128 f128_one(void)
{
    // Example split (check that this matches your Rust constant precisely!):
    return f128_from_raw(
        0x8A5AC6FFFBF001ULL,  // lo
        0xC9BB350DB6ABFD4FULL // hi
    );
}

/**
 * Checks if F128 is zero in GF(2^128).
 */
static inline bool f128_is_zero(const F128 *x)
{
    return (x->lo == 0ULL && x->hi == 0ULL);
}

/**
 * GF(2^128) “addition” is just bitwise XOR.
 */
static inline F128 f128_add(const F128 *a, const F128 *b)
{
    F128 out;
    out.lo = a->lo ^ b->lo;
    out.hi = a->hi ^ b->hi;
    return out;
}

/**
 * In-place GF(2^128) “+=” (XOR).
 */
static inline void f128_add_assign(F128 *a, const F128 *b)
{
    a->lo ^= b->lo;
    a->hi ^= b->hi;
}

/**
 * GF(2^128) bitwise AND.
 */
static inline F128 f128_bitand(const F128 *a, const F128 *b)
{
    F128 out;
    out.lo = a->lo & b->lo;
    out.hi = a->hi & b->hi;
    return out;
}

/**
 * In-place AND.
 */
static inline void f128_bitand_assign(F128 *a, const F128 *b)
{
    a->lo &= b->lo;
    a->hi &= b->hi;
}

/**
 * Example placeholder for random generation.
 * You must implement this to produce a random 128-bit value.
 */
extern F128 f128_rand(void);

/**
 * Multiply two GF(2^128) elements.  This must implement polynomial
 * multiplication modulo the chosen irreducible polynomial for GF(2^128).
 *
 * In your Rust code, this is provided by `mul_128()`.
 * You must implement or link your version below.
 */
extern F128 mul_128(const F128 *a, const F128 *b);

/**
 * F128 * F128 => F128
 */
static inline F128 f128_mul(const F128 *a, const F128 *b)
{
    return mul_128(a, b);
}

/**
 * In-place multiplication: a *= b
 */
static inline void f128_mul_assign(F128 *a, const F128 *b)
{
    *a = mul_128(a, b);
}

/**
 * Frobenius precomputation table, same as your Rust:
 *   pub const FROBENIUS: [[u128; 128]; 128] = ...
 *
 * In C, you might store them as a 2D array of F128, or as raw 64-bit pairs.
 * For illustration, let’s assume we store them as raw pairs:
 */
extern F128 FROBENIUS[128][128];

/**
 * Convert an F128 into a bool array of length 128, little-endian style:
 * bits[0] = least significant bit of (a->lo)
 * bits[127] = most significant bit of (a->hi)
 */
static inline void f128_to_bits(const F128 *a, bool out_bits[128])
{
    // Fill from a->lo:
    for (int i = 0; i < 64; i++) {
        out_bits[i] = (a->lo >> i) & 1ULL;
    }
    // Fill from a->hi:
    for (int i = 0; i < 64; i++) {
        out_bits[64 + i] = (a->hi >> i) & 1ULL;
    }
}

/**
 * Frobenius operator: x -> x^(2^k), with k possibly negative, etc.
 * Mirrors the Rust logic with a shift of k modulo 128.
 */
static inline F128 f128_frob(const F128 *x, int k)
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
            F128 tmp = f128_add(&result, &FROBENIUS[k][i]);
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
extern F128 COBASIS[128];

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
extern F128 COBASIS_FROBENIUS[128][128];

static inline F128 pi_calc(int i, const F128 *twists /* array of size 128 */)
{
    assert(i >= 0 && i < 128);
    F128 ret = f128_zero();
    for (int j = 0; j < 128; j++) {
        // ret += F128::from_raw(COBASIS_FROBENIUS[j][i]) * twists[j];
        F128 tmpmul = f128_mul(&COBASIS_FROBENIUS[j][i], &twists[j]);
        ret = f128_add(&ret, &tmpmul);
    }
    return ret;
}

#endif // F128_H