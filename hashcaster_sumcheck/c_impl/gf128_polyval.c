#include "gf128_polyval.h"
#include <arm_neon.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "field.h"

/******************************************************************************
 * Helper inline NEON intrinsics for polynomial multiplication
 *  pmull => multiply the LOW 64 bits of x and y => 128-bit result
 *  pmull2 => multiply the HIGH 64 bits of x and y => 128-bit result
 *
 ******************************************************************************/

/* 
 * pmull - multiply the low 64 bits (lane 0) of a and b in GF(2^64),
 * returning a 128-bit result in a uint8x16_t. 
 */
static inline poly64_t cast_poly64_from_u64(uint64_t x)
{
    return (poly64_t)x;
}


/* Convert a poly128_t to a uint8x16_t */
static inline uint8x16_t cast_u8q_from_poly128(poly128_t x)
{
    /* Apple clang provides "vreinterpretq_u8_p128(...)": */
    return vreinterpretq_u8_p128(x);
}

/******************************************************************************
 * pmull - multiply the LOW 64 bits of a,b => carryless 128-bit
 ******************************************************************************/
static inline uint8x16_t pmull(uint8x16_t a, uint8x16_t b)
{
    /* Reinterpret a,b as two 64-bit lanes each. Then take lane 0. */
    uint64x2_t AA = vreinterpretq_u64_u8(a);
    uint64x2_t BB = vreinterpretq_u64_u8(b);

    uint64_t a0 = vgetq_lane_u64(AA, 0);
    uint64_t b0 = vgetq_lane_u64(BB, 0);

    poly64_t A0 = cast_poly64_from_u64(a0);
    poly64_t B0 = cast_poly64_from_u64(b0);

    poly128_t result = vmull_p64(A0, B0);
    return cast_u8q_from_poly128(result);
}

/******************************************************************************
 * pmull2 - multiply the HIGH 64 bits of a,b => carryless 128-bit
 ******************************************************************************/
static inline uint8x16_t pmull2(uint8x16_t a, uint8x16_t b)
{
    uint64x2_t AA = vreinterpretq_u64_u8(a);
    uint64x2_t BB = vreinterpretq_u64_u8(b);

    uint64_t a1 = vgetq_lane_u64(AA, 1);
    uint64_t b1 = vgetq_lane_u64(BB, 1);

    poly64_t A1 = cast_poly64_from_u64(a1);
    poly64_t B1 = cast_poly64_from_u64(b1);

    poly128_t result = vmull_p64(A1, B1);
    return cast_u8q_from_poly128(result);
}

/******************************************************************************
 * karatsuba1 => (h, m, l)
 *
 *   h = x.hi * y.hi  => pmull2(x, y)
 *   l = x.lo * y.lo  => pmull(x, y)
 *   m = (x.hi ^ x.lo) * (y.hi ^ y.lo)
 ******************************************************************************/
static inline void karatsuba1(uint8x16_t x, uint8x16_t y,
                              uint8x16_t *h, uint8x16_t *m, uint8x16_t *l)
{
    // x.hi^x.lo => veorq_u8(x, vextq_u8(x,x,8))
    uint8x16_t x_hi = vextq_u8(x, x, 8);
    uint8x16_t x_xor = veorq_u8(x, x_hi);

    // y.hi^y.lo
    uint8x16_t y_hi = vextq_u8(y, y, 8);
    uint8x16_t y_xor = veorq_u8(y, y_hi);

    *m = pmull(x_xor, y_xor);
    *h = pmull2(x, y);
    *l = pmull(x, y);
}

/******************************************************************************
 * karatsuba2 => combine (h, m, l) into a 256-bit product [x23, x01]
 ******************************************************************************/
static inline void karatsuba2(uint8x16_t h, uint8x16_t m, uint8x16_t l,
                              uint8x16_t *out_x23, uint8x16_t *out_x01)
{
    // We do the second step, as in your Rust code:
    //   t0 = m ^ vextq_u8(l, h, 8)
    //   t1 = h ^ l
    //   t = t0 ^ t1
    uint8x16_t lh1 = vextq_u8(l, h, 8);
    uint8x16_t t0 = veorq_u8(m, lh1);
    uint8x16_t t1 = veorq_u8(h, l);
    uint8x16_t t  = veorq_u8(t0, t1);

    // x01 = vextq_u8( vextq_u8(l,l,8), t, 8 ) => [m0^l1^h0^l0, l0]
    uint8x16_t ll = vextq_u8(l, l, 8); 
    uint8x16_t x01 = vextq_u8(ll, t, 8);

    // x23 = vextq_u8( t, vextq_u8(h,h,8), 8 ) => [h1, m1^h0^h1^l1]
    uint8x16_t hh = vextq_u8(h, h, 8);
    uint8x16_t x23 = vextq_u8(t, hh, 8);

    *out_x23 = x23;
    *out_x01 = x01;
}

/******************************************************************************
 * mont_reduce => do the final reduction (like your Rust code).
 * 
 * We assume "POLYVAL" uses the polynomial x^127 + x^126 + x^121 + x^63 + x^62 + x^57,
 * or some variant. You might adapt as needed.
 *****************************************************************************/
static inline uint8x16_t mont_reduce(uint8x16_t x23, uint8x16_t x01)
{

    static const poly128_t poly_basis =
        ((unsigned __int128)1 << 127) |
        ((unsigned __int128)1 << 126) |
        ((unsigned __int128)1 << 121) |
        ((unsigned __int128)1 << 63 ) |
        ((unsigned __int128)1 << 62 ) |
        ((unsigned __int128)1 << 57 );

    uint8x16_t poly = vreinterpretq_u8_p128(poly_basis);
    
    // a = pmull(x01, poly)
    uint8x16_t a = pmull(x01, poly);

    // b = x01 ^ vextq_u8(a,a,8)
    uint8x16_t a_hi = vextq_u8(a, a, 8);
    uint8x16_t b = veorq_u8(x01, a_hi);

    // c = pmull2(b, poly)
    uint8x16_t c = pmull2(b, poly);

    // out = x23 ^ ( c ^ b )
    uint8x16_t out = veorq_u8(x23, veorq_u8(c, b));
    return out;
}


// let (h, m, l) = karatsuba1(transmute(x), transmute(y));
// let (h, l) = karatsuba2(h, m, l);
// transmute(mont_reduce(h, l))

/******************************************************************************
 * gf128_mul_polyval(x_lo, x_hi, y_lo, y_hi):
 *
 ******************************************************************************/
F128 gf128_mul_polyval(F128 x_in, F128 y_in){

    // uint8x16_t X = vld1q_u8((uint8_t*)&x_in);

    // // 2) Convert y_lo,y_hi into a uint8x16_t
    uint8_t y_bytes[16], x_bytes[16];
    for (int i=0; i<8; i++) {
        y_bytes[i]   = (uint8_t)((y_in.low >> (8*i)) & 0xFF);
        y_bytes[8+i] = (uint8_t)((y_in.high >> (8*i)) & 0xFF);
    }
    for (int i=0; i<8; i++) {
        x_bytes[i]   = (uint8_t)((x_in.low >> (8*i)) & 0xFF);
        x_bytes[8+i] = (uint8_t)((x_in.high >> (8*i)) & 0xFF);
    }
    uint8x16_t Y = vld1q_u8(y_bytes);
    uint8x16_t X = vld1q_u8(x_bytes);

    // 3) karatsuba1 => (h, m, l)
    uint8x16_t h, m, l;
    karatsuba1(X, Y, &h, &m, &l);

    // 4) karatsuba2 => combine => (x23, x01)
    uint8x16_t x23, x01;
    karatsuba2(h, m, l, &x23, &x01);

    // 5) mont_reduce(x23, x01)
    uint8x16_t reduced = mont_reduce(x23, x01);

    // 6) Convert reduced (uint8x16_t) back to F128 => [lo, hi].
    F128 out = {0ULL, 0ULL}; 
    out.low  = 0ULL;
    out.high = 0ULL;
    uint8_t tmp[16];
    vst1q_u8(tmp, reduced);
    for (int i=0; i<8; i++) {
        out.low  |= ((uint64_t)tmp[i]) << (8*i);
        out.high |= ((uint64_t)tmp[8+i]) << (8*i);
    }
    return out;
}

/******************************************************************************
 * v_movemask_epi8:
 * If you need a direct implementation, here's a possible approach that 
 * extracts the top bit of each byte in 'input' to produce a 16-bit mask.
 ******************************************************************************/
uint16_t v_movemask_epi8(uint8x16_t input)
{
    uint8_t bytes[16];
    vst1q_u8(bytes, input);

    uint16_t mask = 0;
    for (int i = 0; i < 16; i++) {
        // top bit => (bytes[i] & 0x80)
        mask |= (uint16_t)(((bytes[i] & 0x80) >> 7) << i);
    }
    return mask;
}