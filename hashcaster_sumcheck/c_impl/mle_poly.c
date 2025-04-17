#include "mle_poly.h"
#include "field.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include "debug_utils.h"
/******************************************************************************
 * points_default
 *   Returns an empty Points with len=0, elems=NULL.
 ******************************************************************************/
Points points_default(void)
{
    Points p;
    p.elems = NULL;
    p.len = 0;
    return p;
}


/******************************************************************************
 * points_eq_eval
 *   eq_eval(a,b) = ∏_{i=0..k-1} (1 + a[i] + b[i]),
 *   where k = min(a.len, b.len) or we might require they match in len.
 *   We'll assume we multiply up to the smaller len, or require lens match.
 ******************************************************************************/
F128 points_eq_eval(const Points* a, const Points* b)
{
    // If you want to handle mismatched lens, adapt logic. We'll assume they match for simplicity.
    size_t n = (a->len < b->len) ? a->len : b->len;

    F128 acc = f128_one();
    F128 one = f128_one();

    for (size_t i = 0; i < n; i++) {
        // sum = 1 + a[i] + b[i]
        F128 s = f128_add(one, f128_add(a->elems[i], b->elems[i]));
        // acc *= s
        acc = f128_mul(acc, s);
    }

    // If you want "empty => 1" logic, it's already covered by n=0 => loop doesn't run => acc=1.
    return acc;
}



/* eq_poly:
 *   - pt: array of F128, length = l
 *   - l: length of pt
 *   - out_len: (output) the size of the returned array = 2^l
 * returns newly allocated array of F128
 */
// DYNAMIC/HEAP MEMORY ALLOCATED WITHIN
MLE_POLY* points_to_eq_poly(const Points* points)
{
    // The result length is 2^l
    size_t n = ((size_t)1 << points->len);
    // *out_len = n;

    // allocate ret array
    MLE_POLY *ret = malloc(sizeof(MLE_POLY));
    if (!ret) {
        fprintf(stderr, "eq_poly: out of memory\n");
        exit(1);
    }
    F128 *coeff = malloc(n * sizeof(F128));
    if (!coeff) {
        fprintf(stderr, "eq_poly: out of memory\n");
        exit(1);
    }

    // 1) coeff[0] = 1, and presumably set the rest to 0
    //    but we'll do that in the loop's logic. We'll at least zero them:
    for (size_t i = 0; i < n; i++) {
        coeff[i].low  = 0ULL;
        coeff[i].high = 0ULL;
    }
    // coeff[0] = 1
    coeff[0] = f128_one();

    // 2) for i in [0..l):
    for (size_t i = 0; i < points->len; i++) {
        // half = 1 << i
        size_t half = ((size_t)1 << i);

        // for j in [0..half):
        //   coeff[j + half] = pt[i] * coeff[j];
        //   coeff[j] += coeff[j + half];
        for (size_t j = 0; j < half; j++) {
            F128 temp = f128_mul(points->elems[i], coeff[j]);
            coeff[j + half] = temp;
            coeff[j] = f128_add(coeff[j], temp);
        }
    }

    ret->coeffs = coeff;
    ret->len = n;

    return ret;
}

F128 eq_eval(const F128 *left, const F128 *right, size_t length)
{
    // Start the accumulator at '1' in GF(2^128).
    F128 acc = f128_one();

    for (size_t i = 0; i < length; i++) {
        // sum = 1 + x + y  (where + is XOR in GF(2^128)).
        F128 sum = f128_add(f128_add(f128_one(), left[i]), right[i]);
        // acc = acc * sum
        acc = f128_mul(acc, sum);
    }

    return acc;
}

F128* clone_f128s(const F128 *src, size_t length)
{
    F128* dst;
    // dst.length = src->length;
    dst = (F128 *)malloc(length * sizeof(F128));
    if (!dst) {
        // handle error as needed
        return NULL;
    }

    memcpy(dst, src, length * sizeof(F128));
    return dst;
}

Points* clone_points(const Points *src)
{
    Points* dst;
    dst = (Points *)malloc(sizeof(Points));
    if (!dst) {
        // handle error as needed
        return NULL;
    }
    dst->len = src->len;
    dst->elems = clone_f128s(src->elems, src->len);

    return dst;
}


// a 2D array of F128s. Each row is a snapshot of the orbit
INVERSE_ORBIT_POINTS* to_f128_inv_orbit(const Points* input_points)
{
    INVERSE_ORBIT_POINTS *orbits = (INVERSE_ORBIT_POINTS *)malloc(sizeof(INVERSE_ORBIT_POINTS));
    if (!orbits) {
        // handle error if needed
        return NULL;
    }

    // We'll produce 128 snapshots
    const size_t N = 128;
    orbits->len = N;
    orbits->array_of_points = ( Points *)malloc(N * sizeof(Points));
    if (!orbits->array_of_points) {
        // handle error if needed
        return NULL;
    }


    F128* tmp = clone_f128s(input_points->elems, input_points->len);

    // For i in [0..128):
    for (size_t i = 0; i < N; i++) {
        // Square each point in 'tmp'
        for (size_t j = 0; j < input_points->len; j++) {
            tmp[j] = f128_mul(tmp[j], tmp[j]);
        }

        // orbits->array_of_points[i] = *clone_points(tmp);// mem leak here
        orbits->array_of_points[i].elems = clone_f128s(tmp, input_points->len);
        orbits->array_of_points[i].len = input_points->len;
    }

    // Reverse snapshots in place
    for (size_t i = 0; i < N/2; i++) {
        Points t = orbits->array_of_points[i];
        orbits->array_of_points[i] = orbits->array_of_points[N - 1 - i];
        orbits->array_of_points[N - 1 - i] = t;
    }

    free(tmp); // free the temporary array
    return orbits;
}

// returns a list of polynomials, each of which is a sequence of coefficients
MLE_POLY_SEQUENCE* points_to_eq_poly_sequence(const Points *points)
{
    // The result array will have (points_len + 1) polynomials
    size_t array_size = points->len + 1;
    MLE_POLY_SEQUENCE *seq = (MLE_POLY_SEQUENCE *)malloc(sizeof(MLE_POLY_SEQUENCE));
    if (!seq) {
        // handle allocation failure if needed
        return NULL;
    }
    MLE_POLY *polynomials = (MLE_POLY *)malloc(array_size * sizeof(MLE_POLY));
    if (!polynomials) {
        // handle allocation failure if needed
        return NULL;
    }
    seq->mle_poly = polynomials;
    seq->len = array_size;

    // 1) polynomials[0] = a polynomial with 1 coefficient = [f128_one()]
    polynomials[0].len = 1;
    polynomials[0].coeffs = (F128 *)malloc(sizeof(F128));
    if (!polynomials[0].coeffs) {
        // handle error
        free(polynomials);
        // *out_count = 0;
        return NULL;
    }
    polynomials[0].coeffs[0] = f128_one();

    // 2) Build polynomials[1..points_len], iterating over points in reverse order
    //    i = 0..(points_len-1)
    for (size_t i = 0; i < points->len; i++) {
        // "multiplier" = points[ points_len - 1 - i ]
        F128 multiplier = points->elems[points->len - 1 - i];

        // "previous" = polynomials[i]
        size_t prev_len = polynomials[i].len;
        F128 *prev_coeffs = polynomials[i].coeffs;

        // new polynomial size = 1 << (i+1)
        size_t new_len = (size_t)1 << (i + 1);
        F128 *new_coeffs = (F128 *)calloc(new_len, sizeof(F128));
        if (!new_coeffs) {
            // handle error
            // free everything up to i
            for (size_t k = 0; k <= i; k++) {
                free(polynomials[k].coeffs);
            }
            free(polynomials);
            // *out_count = 0;
            return NULL;
        }

        // Fill new_coeffs
        // for j in [0..prev_len) = [0..(1<<i)), do:
        //   let prev = prev_coeffs[j]
        //   let multiplied = multiplier * prev
        //   new_coeffs[2*j]   = prev + multiplied
        //   new_coeffs[2*j+1] = multiplied
        for (size_t j = 0; j < prev_len; j++) {
            F128 p = prev_coeffs[j];
            F128 mulp = f128_mul(multiplier, p);

            // new_coeffs[2*j] = p + mulp
            // new_coeffs[2*j+1] = mulp
            size_t idx0 = 2*j;
            size_t idx1 = 2*j + 1;
            new_coeffs[idx0] = f128_add(p, mulp);
            new_coeffs[idx1] = mulp;
        }

        polynomials[i+1].len = new_len;
        polynomials[i+1].coeffs = new_coeffs;
    }

    // *out_count = array_size;
    return seq;
}


// implement a function to compute if 2 MLE_POLY are equal
bool mle_poly_eq(const MLE_POLY *a, const MLE_POLY *b)
{
    if (a->len != b->len) {
        return false;
    }

    for (size_t i = 0; i < a->len; i++) {
        if (!f128_eq(a->coeffs[i], b->coeffs[i])) {
            return false;
        }
    }

    return true;
}

void drop_top_bit(uint8_t x, uint8_t *result, uint8_t *bit_index) {
    if (x == 0) {
        *result = 0;
        *bit_index = 0;
        return;
    }

    // Compute leading zeros in 8-bit integer
    uint8_t lz = __builtin_clz(x) - 24; // __builtin_clz operates on 32-bit ints
    *bit_index = 7 - lz;
    *result = x & ~(1 << *bit_index);
}

// The eq_sums function implementation
Points* eq_sums(const  MLE_POLY *eq_poly) {
    // Ensure eq_len is divisible by 16
    if (eq_poly->len % 16 != 0) {
        return NULL;  // or handle as error
    }

    size_t num_blocks = eq_poly->len / 8;
    size_t total_results = 256 * num_blocks;

    // Allocate the result array
    F128* elems = (F128*)malloc(total_results * sizeof(F128));
    if (!elems) {
        return NULL;  // handle malloc failure
    }
    Points* result = (Points*)malloc(sizeof(Points));
    if (!result) {
        free(result);
        return NULL;
    }
    result->len = total_results;
    result->elems = elems;

    // Process the polynomial in blocks of 8 coefficients each
    for (size_t block_idx = 0; block_idx < num_blocks; ++block_idx) {
        const F128* block = &eq_poly->coeffs[block_idx * 8];

        // Pointer to current block results
        F128* block_sums = &result->elems[block_idx * 256];

        // Initialize the first sum (empty subset) to zero
        block_sums[0] = f128_zero();

        // Iterate over subsets from 1 to 255
        for (uint16_t subset = 1; subset < 256; ++subset) {
            uint8_t sum_idx, bit_idx;

            // Decompose subset index
            drop_top_bit((uint8_t)subset, &sum_idx, &bit_idx);

            // Compute sum for the current subset
            block_sums[subset] = f128_add(block[bit_idx], block_sums[sum_idx]);
        }
    }

    return result;
}


// This function takes a 16-byte array (x),
// and returns an int representing the combined MSB mask.
int cpu_v_movemask_epi8(const uint8_t* x) {
    // 'ret' accumulates the final 16-bit pattern of MSBs,
    // stored as a 32-bit integer.
    int ret = 0;

    // Loop from 0..16, shifting 'ret' left by 1 each iteration,
    // then add the top bit of x[15-i].
    for (int i = 0; i < 16; i++) {
        // Shift left by 1 to make room for the next bit.
        ret <<= 1;

        // Extract the MSB of x[15 - i].
        // (x[15 - i] >> 7) yields the top bit (0 or 1),
        // then we add it to 'ret'.
        ret += (x[15 - i] >> 7);
    }

    return ret;
}


void v_slli_epi64_c(int K, const uint8_t *x)
{
    // 1) Load the 16-byte array into a 128-bit NEON register of two 64-bit integers.
    //    We treat x as an array of two int64_t.
    int64x2_t data = vld1q_s64((const int64_t*) x);

    // 2) Shift each 64-bit integer in data left by K bits.
    //    vshlq_n_s64 can only accept a compile-time constant, 
    //    so if K is truly constant at compile time, you do vshlq_n_s64(data, K).
    //    Otherwise, use the more general vshlq_s64 with a vector shift count.
    //    For demonstration, we assume K is a compile-time or small value.
    //    We'll do a switch or so:
    int64x2_t result;
    switch(K) {
        case 0:  result = data;               break;
        case 1:  result = vshlq_n_s64(data,1);  break;
        case 2:  result = vshlq_n_s64(data,2);  break;
        // ...
        // Up to 63. For brevity, not enumerating all cases. 
        // In practice, you might generate code or 
        // use vshlq_s64 for a runtime shift with a vector.
        default: 
            // Fallback approach for a runtime shift count:
        {
            // NEON requires a vector for vshlq_s64 if not constant. 
            // We'll create that vector from K.
            int64x2_t shift_vec = vdupq_n_s64(K);
            // Use vshlq_s64 that shifts each 64-bit int by shift_vec's bits. 
            result = vshlq_s64(data, shift_vec);
            break;
        }
    }

    // 3) Store the shifted 128-bit result back into out[16].
    vst1q_s64((int64_t*) x, result);
}

Evaluations* restrict_polynomials(
    const MLE_POLY *polys, // array of polynomials
    size_t N,                 // how many polynomials
    const Points *challenges, // challenge points
    size_t dims               // dimension
) {
    // Step 1: check all polynomials have len == 1 << dims
    for(size_t i=0; i<N; i++){
        // Check length 
        assert(polys[i].len <= (1ULL << dims));
    }

    // Step 2: number of challenges
    size_t num_challenges = challenges->len;
    assert(num_challenges <= dims);

    size_t chunk_size = (1ULL << num_challenges);
    size_t num_chunks = (1ULL << (dims - num_challenges));

    MLE_POLY* eq = points_to_eq_poly(challenges); 

    // Step 5: eq_len = eq.len and must be divisible by 16
    size_t eq_len = eq->len; 
    assert((eq_len % 16)==0);

    // Step 6: eq_sums
    Points * eq_sums_arr = eq_sums(eq);

    // Step 7: create ret array with size = (num_chunks * 128 * N)
    size_t ret_size = num_chunks * 128ULL * N; 
    // We'll allocate 
    F128 *ret = calloc(ret_size, sizeof(F128));

    // Step 8: the main loop from Rust:
    for(size_t q=0; q<N; q++){
        // ret_chunk_start = q * 128 * num_chunks
        size_t ret_chunk_start = q * 128ULL * num_chunks;
        for(size_t i=0; i<num_chunks; i++){
            // for j in [0.. eq_len/16)
            size_t blocks = eq_len/16;
            for(size_t j=0; j<blocks; j++){
                // b1 = j * 512
                size_t b1 = j*512;
                size_t b2 = b1 + 256;
                // size_t b3 = b2 + 256;
                const F128 *v0 = &(eq_sums_arr->elems[b1]);
                const F128 *v1 = &(eq_sums_arr->elems[b2]);

                // We interpret each F128 as 16 bytes. 
                // Then we process each byte.
                // uint8_t* coeffs = (polys[q].coeffs );
                

                const F128 * poly_block_16 = &(polys[q].coeffs[ i*chunk_size + j*16 ]);


                // We'll do a for s in 0..16
                for(size_t s=0; s<16; s++){
                    // Build array t of 16 bytes
                    uint8_t t[16];
                    for(size_t idx=0; idx<16; idx++){
                        uint8_t* poly_poly_bytes = (uint8_t*)(&poly_block_16[idx]);
                        t[idx] = poly_poly_bytes[s];
                    }

                    // for u in 0..8
                    for(size_t u=0; u<8; u++){
                        size_t bits = cpu_v_movemask_epi8(t);
                        
                        // Now do ret[...] += v0[bits & 255] + v1[(bits>>8)&255]
                        // ret index = ret_chunk_start + (s*8 + 7 - u)*num_chunks + i
                        size_t ret_index = ret_chunk_start + (s*8 + (7ULL - u))*num_chunks + i;
                        
                        // print the ret_index
                        //  printf("ret_index: %zu\n", ret_index);
                        
                        // We'll do a field add 
                        F128 sum1 = f128_add(v0[bits & 0xFF], v1[(bits>>8)&0xFF]);
                        // add it to ret[ret_index]

                        ret[ret_index] = f128_add(ret[ret_index], sum1);

                        // updates the t in place
                        v_slli_epi64_c( 1, t);

                    }
                }
            }
        }
    }

    // Return final evaluations
    Points* res;
    res = (Points*)malloc(sizeof(Points));
    res->elems = ret;
    res->len = ret_size;
    return res;
}


// Evaluate polynomial at given 'points'
F128 evaluate_at(
    const MLE_POLY *poly, 
    const Points *points
)
{
    // 1) ensure we have the right length
    size_t required_len = (1ULL << points->len);
    // In real code, handle if poly->len != required_len
    assert(poly->len == required_len);

    // 2) compute eq polynomial
    MLE_POLY* eq = points_to_eq_poly(points);
    // eq.len should also be required_len
    assert(eq->len == required_len);

    // 3) sum up each product of poly->coeffs[i]* eq.coeffs[i]
    F128 sum = f128_zero(); // start with zero
    for(size_t i = 0; i < required_len; i++){
        // multiply
        F128 prod = f128_mul(poly->coeffs[i], eq->coeffs[i]);
        // add to sum
        sum = f128_add(sum, prod);
    }

    // 4) free eq polynomial if allocated
    free(eq->coeffs);

    // 5) return sum
    return sum;
}


MLE_POLY_SEQUENCE* mle_sequence_new(size_t sequence_len, size_t poly_len)
{
    MLE_POLY_SEQUENCE* seq = (MLE_POLY_SEQUENCE*)malloc(sizeof(MLE_POLY_SEQUENCE));
    if (!seq) {
        // handle error
        return NULL;
    }

    seq->len = sequence_len;
    seq->mle_poly = (MLE_POLY*)malloc(sequence_len * sizeof(MLE_POLY));
    if (!seq->mle_poly) {
        free(seq);
        return NULL;
    }

    for (size_t i = 0; i < sequence_len; i++) {
        seq->mle_poly[i].len = poly_len;
        seq->mle_poly[i].coeffs = (F128*)calloc(poly_len , sizeof(F128));
        if (!seq->mle_poly[i].coeffs) {
            // handle error
            for (size_t j = 0; j < i; j++) {
                free(seq->mle_poly[j].coeffs);
            }
            free(seq->mle_poly);
            free(seq);
            return NULL;
        }
    }

    return seq;
}

void mle_sequence_free(MLE_POLY_SEQUENCE* seq){
    if (seq == NULL) {
        return;
    }
    for (size_t i = 0; i < seq->len; i++) {
        free(seq->mle_poly[i].coeffs);
    }
    free(seq->mle_poly);
    free(seq);
}
