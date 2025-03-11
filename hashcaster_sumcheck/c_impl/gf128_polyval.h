#ifndef GF128_POLYVAL_H
#define GF128_POLYVAL_H

#include <stdint.h>
#include <stdbool.h>
#include  "field.h"

/*
 * On Apple Silicon (AArch64), we can use <arm_neon.h> for NEON intrinsics.
 */
#include <arm_neon.h>

#ifdef __cplusplus
extern "C" {
#endif


/*
 * This function multiplies two 128-bit values in the GF(2^128) "POLYVAL" field 
 * using Karatsuba + Montgomery reduction. The Rust code does:
 *   pub fn mul_128(x: u128, y:u128) -> u128 { ... }
 * We'll do the same in C, returning F128.
 */
F128 gf128_mul_polyval(F128 x, F128 y);

/*
 * If needed, define v_movemask_epi8 in C. 
 */
uint16_t v_movemask_epi8(uint8x16_t input);

#ifdef __cplusplus
}
#endif

#endif /* GF128_POLYVAL_H */