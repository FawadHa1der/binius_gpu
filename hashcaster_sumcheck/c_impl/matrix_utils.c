#include "matrix_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/******************************************************************************
 * 1) log2_exact
 *
 * Rust:
 *   pub fn log2_exact(mut x: usize) -> usize {
 *       let mut c = 0;
 *       loop {
 *           if x % 2 == 1 {
 *               if x == 1 {return c;} else {panic!("Not an exact power of 2.");}
 *           }
 *           c += 1;
 *           x >>= 1;
 *       }
 *   }
 ******************************************************************************/
unsigned log2_exact(unsigned long long x)
{
    unsigned c = 0;
    while (1) {
        if ((x & 1ULL) == 1ULL) {
            if (x == 1ULL) {
                return c;
            } else {
                fprintf(stderr, "Not an exact power of 2.\n");
                abort(); // or handle error differently
            }
        }
        c++;
        x >>= 1ULL;
    }
    // unreachable
}

/******************************************************************************
 * 2) "F128"
 *
 * We'll store as (low, high). We define the placeholders:
 ******************************************************************************/
F128 f128_new(uint64_t low, uint64_t high)
{
    F128 out;
    out.low = low;
    out.high = high;
    return out;
}

void f128_get(const F128 *x, uint64_t *low, uint64_t *high)
{
    *low = x->low;
    *high = x->high;
}


/******************************************************************************
 * 3) Conversions: 
 *    binary128_to_bits(x: F128) -> bool[128]
 ******************************************************************************/
void binary128_to_bits(F128 x, bool out_bits[128])
{
    // interpret x as 16 bytes, little-endian: the Rust code does:
    //   let bytes = cast::<u128, [u8;16]>( x.val() )
    // We'll do it manually:
    uint8_t bytes[16];
    // fill first 8 from x.low
    for (int i = 0; i < 8; i++) {
        bytes[i] = (uint8_t)((x.low >> (8*i)) & 0xFF);
    }
    // next 8 from x.high
    for (int i = 0; i < 8; i++) {
        bytes[8 + i] = (uint8_t)((x.high >> (8*i)) & 0xFF);
    }

    // Then push bits i in [0..8), j in [0..16)
    int idx = 0;
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 8; j++) {
            out_bits[idx++] = ((bytes[i] & (1U << j)) != 0);
        }
    }
}

/******************************************************************************
 * 4) u128_idx(x: &u128, i: usize) -> bool
 *    We'll do a "f128_idx(...)"
 ******************************************************************************/
bool f128_idx(const F128 *x, unsigned i)
{
    // same approach: interpret x->(low, high) as 16 bytes
    // find the byte index and bit index
    uint8_t bytes[16];
    for (int k = 0; k < 8; k++) {
        bytes[k] = (uint8_t)((x->low >> (8*k)) & 0xFF);
    }
    for (int k = 0; k < 8; k++) {
        bytes[8 + k] = (uint8_t)((x->high >> (8*k)) & 0xFF);
    }
    unsigned byte_idx = i >> 3;      // i / 8
    unsigned bit_idx  = i & 7;       // i mod 8
    return (bytes[byte_idx] & (1U << bit_idx)) != 0;
}

/******************************************************************************
 * 5) _u128_from_bits(bits: &[bool]) -> u128
 *    We'll define "u128_from_bits(bits, l, &low, &high)"
 ******************************************************************************/
void u128_from_bits(const bool *bits, size_t length, uint64_t *low, uint64_t *high)
{
    assert(length <= 128);
    uint64_t lo = 0ULL;
    uint64_t hi = 0ULL;
    for (size_t i = 0; i < length; i++) {
        if (bits[i]) {
            if (i < 64) {
                lo |= (1ULL << i);
            } else {
                hi |= (1ULL << (i - 64));
            }
        }
    }
    *low  = lo;
    *high = hi;
}

/******************************************************************************
 * 6) Matrix structure
 *    pub struct Matrix { pub cols: Vec<F128> }
 *
 *    We'll store them in "F128 cols[128]"
 ******************************************************************************/
Matrix matrix_new(const F128 *cols)
{
    Matrix m;
    for (int i = 0; i < 128; i++) {
        m.cols[i] = cols[i];
    }
    return m;
}

/*
 * matrix_apply(&self, vec)
 * Rust code:
 *   let mut ret = F128::new(0);
 *   let vec_bits = binary128_to_bits(vec);
 *   for i in 0..128 {
 *       if vec_bits[i] { ret += matrix[i] }
 *   }
 *   ret
 */
F128 matrix_apply(const Matrix *m, F128 vec)
{
    bool bits[128];
    binary128_to_bits(vec, bits);

    F128 ret = f128_zero();
    for (int i = 0; i < 128; i++) {
        if (bits[i]) {
            ret = f128_add(ret, m->cols[i]);
        }
    }
    return ret;
}

/* 
 * matrix_compose(&self, b)
 * for i in 0..128 { ret.push(self.apply(b.cols[i])); }
 */
Matrix matrix_compose(const Matrix *m, const Matrix *b)
{
    Matrix out;
    for (int i = 0; i < 128; i++) {
        out.cols[i] = matrix_apply(m, b->cols[i]);
    }
    return out;
}

/* 
 * matrix_diag()
 *  builds an identity matrix in GF(2): column i has only bit i set
 */
Matrix matrix_diag(void)
{
    Matrix m;
    for (int i = 0; i < 128; i++) {
        m.cols[i] = f128_basis(i);
    }
    return m;
}

/* 
 * matrix_swap_cols(&mut self, i, j)
 *  swap self.cols[i], self.cols[j]
 */
void matrix_swap_cols(Matrix *m, size_t i, size_t j)
{
    F128 tmp = m->cols[j];
    m->cols[j] = m->cols[i];
    m->cols[i] = tmp;
}

/* 
 * matrix_triang(&mut self, i, j)
 *   self.cols[j] += self.cols[i]
 */
void matrix_triang(Matrix *m, size_t i, size_t j)
{
    m->cols[j] = f128_add(m->cols[j], m->cols[i]);
}

/******************************************************************************
 * 7) matrix_inverse => Option<Self>
 *    We'll do a C function returning bool. If false => no inverse (None). If true => *out_inv is the inverse (Some).
 *
 * Rust code:
 *    let mut a = self.clone();
 *    let mut b = diag();
 *    for i in 0..128 {
 *       if !u128_idx(a.cols[i].val(), i) => find pivot or return None
 *       ...
 *    }
 ******************************************************************************/
bool matrix_inverse(const Matrix *in_m, Matrix *out_inv)
{
    // Make local copies:
    Matrix a = *in_m;          // 'a = self.clone();'
    Matrix b = matrix_diag();  // 'b = Self::diag();'

    for (int i = 0; i < 128; i++) {
        int j = i;
        // if ! f128_idx(&a.cols[i], i), find pivot below
        if (!f128_idx(&a.cols[i], i)) {
            while (1) {
                j++;
                if (j == 128) {
                    // no pivot found => singular
                    return false;
                }
                if (f128_idx(&a.cols[j], i)) {
                    matrix_swap_cols(&a, i, j);
                    matrix_swap_cols(&b, i, j);
                    break;
                }
            }
        }
        // now pivot in row i => eliminate below
        for (;;) {
            j++;
            if (j == 128) break;
            if (f128_idx(&a.cols[j], i)) {
                matrix_triang(&a, i, j);
                matrix_triang(&b, i, j);
            }
        }
    }

    // "a is now under-diagonal"
    // for i in [1..128], for j in [0..i], if f128_idx(&a.cols[j], i) => triang(i,j) 
    for (int i = 1; i < 128; i++) {
        for (int j = 0; j < i; j++) {
            if (f128_idx(&a.cols[j], i)) {
                matrix_triang(&a, i, j);
                matrix_triang(&b, i, j);
            }
        }
    }

    // success => out_inv = b
    *out_inv = b;
    return true;
}