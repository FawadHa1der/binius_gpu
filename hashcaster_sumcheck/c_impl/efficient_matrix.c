
#include "efficient_matrix.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <memory.h>

// 1) Constants
#define NUM_ROWS 128
#define NUM_COLS 128

// The number of 16-bit chunks in a 128-bit value
#define CHUNK_COUNT (128 / 16) // which is 8

// 2) We'll define a struct for each (idx_u64, shift) pair
typedef struct {
    size_t idx_u64;
    size_t shift;
} IdxU64ShiftPair;

// 3) We'll initialize the static array of length CHUNK_COUNT
//    This replicates the Rust logic:
//      for i in 0..CHUNK_COUNT:
//         table[i] = (i/4, 16*(i%4))
static const IdxU64ShiftPair IDX_U64_SHIFT[CHUNK_COUNT] = {
    // i=0 -> (0/4=0, 16*(0%4)=0)
    {0, 0},
    // i=1 -> (1/4=0, 16*(1%4)=16)
    {0, 16},
    // i=2 -> (2/4=0, 16*(2%4)=32)
    {0, 32},
    // i=3 -> (3/4=0, 16*(3%4)=48)
    {0, 48},

    // i=4 -> (4/4=1, 16*(4%4)=0)
    {1, 0},
    // i=5 -> (5/4=1, 16*(5%4)=16)
    {1, 16},
    // i=6 -> (6/4=1, 16*(6%4)=32)
    {1, 32},
    // i=7 -> (7/4=1, 16*(7%4)=48)
    {1, 48}
};




// A default for "EfficientMatrix"
inline EfficientMatrix efficient_matrix_default(void)
{
    EfficientMatrix em;
    for(int i=0; i<EFFICIENT_MATRIX_SIZE; i++){
        em.data[i]= f128_zero();
    }
    return em;
}


// For "from_frobenius_inv_lc" 
EfficientMatrix from_frobenius_inv_lc(const F128 gammas[NUM_COLS])
{
    // 1) compute minus_indices => (NUM_COLS - i) % NUM_COLS
    uint32_t minus_indices[NUM_COLS];
    for(int i=0; i<NUM_COLS; i++){
        minus_indices[i] = (uint32_t)((NUM_COLS - i) % NUM_COLS);
    }

    // 2) build ret array of size NUM_COLS => each ret[j] = sum_{i} gammas[i]* from(FROBENIUS[minus_i][j])
    F128 ret[NUM_COLS];
    for(int j=0; j<NUM_COLS; j++){
        // start sum=0
        F128 sum= f128_zero();
        for(int i=0; i<NUM_COLS; i++){
            F128 tmpVal= FROBENIUS[ minus_indices[i] ][ j];

            // partial= gammas[i]* tmpVal
            F128 partial= f128_mul(gammas[i], tmpVal);
            // sum += partial
            sum= f128_add(sum, partial);
        }
        ret[j]= sum;
    }

    // 3) call from_cols(ret)
    return from_cols(ret);
}

EfficientMatrix from_rows(const F128 rows[NUM_ROWS])
{
    uint64_t colvals[256];
    for(int i=0;i<256;i++){
        colvals[i]= 0ULL; 
    }

    // 2) reinterpret rows as [NUM_ROWS][16] bytes. We'll do a union or cast:
    const uint8_t (*byte_rows)[16] = (const uint8_t (*)[16]) rows;

    // 3) chunk in sets of 16 row
    // each chunk => chunk_idx in [0..(NUM_ROWS/16)] => 8 if 128 rows
    int chunk_count = NUM_ROWS/16;  // =8
    for(int chunk_idx=0; chunk_idx<chunk_count; chunk_idx++){
        // retrieve (idx_u64, shift)
        size_t idx_u64= IDX_U64_SHIFT[chunk_idx].idx_u64;
        size_t shift  = IDX_U64_SHIFT[chunk_idx].shift;

        // each chunk => 16 rows => chunk => &byte_rows[16*chunk_idx .. 16*(chunk_idx+1)-1]
        const uint8_t (*chunk)[16]= &byte_rows[16*chunk_idx];

        // For i in [0..16], gather a temporary array t of length=16 => t[j]= chunk[j][i]
        for(int i=0; i<16; i++){
            uint8_t t[16];
            for(int j=0;j<16;j++){
                t[j]= chunk[j][i];
            }

            // Now we process 8 bits in t => each iteration => do movemask => shift => update colvals
            for(int b=0; b<8; b++){
                // bits= cpu_v_movemask_epi8(t)
                int bits= cpu_v_movemask_epi8(t);
                uint64_t mask= ((uint64_t)bits)<< shift;

                // col index= 2*(8*i + (7-b)) + idx_u64
                // from your snippet => cols[ 2*(8*i+ 7-b)+ idx_u64 ].fetch_xor(mask)
                // We'll do a direct XOR into colvals:
                int col_index= 2*(8*i + (7 - b)) + (int)idx_u64;
                colvals[col_index]^= mask;

                // shift all bits in t left => v_slli_epi64_c(1, t)
                v_slli_epi64_c(1, t);
            }
        }
    }

    // 4) finalize colvals => interpret them as 128 columns => then call from_cols
    // We'll define an array of F128 => 128 columns
    F128 cols[NUM_COLS];
    memset(cols, 0, sizeof(cols));

    // we interpret colvals => 2 per column?? 
    // In the original logic: colvals has 256 entries => each pair is idxU64=0 or 1?
    // We'll do a naive approach: "for j in [0..128], low= colvals[2*j], high= colvals[2*j+1]"
    for(int j=0; j<128; j++){
        F128 c;
        c.low= colvals[2*j + 0];
        c.high= colvals[2*j + 1];
        cols[j]= c;
    }

    // Now from_cols
    return from_cols(cols);
}

static inline EfficientMatrix efficient_matrix_zero(void)
{
    EfficientMatrix em;
    for(int i=0; i<EFFICIENT_MATRIX_SIZE; i++){
        em.data[i].low=0;
        em.data[i].high=0;
    }
    return em;
}
// from_cols => we produce 256 sums for each group of 8 columns => total groups= NUM_COLS/8=16 => total=16*256=4096
EfficientMatrix from_cols(const F128 cols[NUM_COLS])
{
    EfficientMatrix em= efficient_matrix_zero();

    int group_count= NUM_COLS/8; // typically 16
    for(int g=0; g< group_count; g++){
        // sums => em.data[ g*256 .. g*256+255 ]
        int base= g*256;
        // sums[0]= 0
        em.data[ base+0 ]= (F128){0,0};

        // for i in [1..255], do "drop_top_bit((uint8_t)i, &sum_idx, &row_idx)"
        for(int i=1; i<256; i++){
            uint8_t res, bit_idx;
            drop_top_bit((uint8_t)i, &res, &bit_idx);

            // sums[i]= sums[ res ] + cols[ 8*g + bit_idx ]
            // note that res is a 'uint8_t' => sums[base+ res]
            // bit_idx => which of the 8 columns in group
            F128 partial_sum= f128_add( em.data[ base+ res ], cols[ 8*g + bit_idx ]);
            em.data[ base+ i ]= partial_sum;
        }
    }

    return em;
}

F128 efficient_matrix_apply(const EfficientMatrix *matrix, F128 rhs)
{
    // interpret 'rhs' as 16 bytes => for demonstration, we'll do a union
    union {
        F128 f;
        unsigned char bytes[16];
    } alias;
    alias.f= rhs;

    F128 acc= f128_zero();
    // for i in [0..16], let byte= alias.bytes[i], index= byte + 256*i
    // acc= acc + matrix->data[index]
    for(int i=0;i<16;i++){
        unsigned char byteVal= alias.bytes[i];
        int index= byteVal + 256*i;
        acc= f128_add(acc, matrix->data[index]);
    }
    return acc;
}