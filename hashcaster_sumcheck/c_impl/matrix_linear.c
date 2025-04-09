#include "matrix_linear.h"
#include <stdlib.h>
#include <string.h>
MatrixLinear matrix_linear_new(size_t n_in, size_t n_out, const F128* entries, size_t len_entries)
{
    // Rust: "assert_eq!(entries.len(), n_in*n_out, 'Invalid matrix dimensions');"
    assert(len_entries == n_in*n_out && "Invalid matrix dimensions in matrix_linear_new");

    MatrixLinear mat;
    mat.n_in = n_in;
    mat.n_out= n_out;

    // allocate memory for mat.entries
    mat.entries= (F128*)malloc(sizeof(F128)* (n_in*n_out));
    // copy from user array
    memcpy(mat.entries, entries, sizeof(F128)* (n_in*n_out));

    return mat;
}

// Freed
void matrix_linear_free(MatrixLinear *mat)
{
    if(mat->entries){
        free(mat->entries);
        mat->entries= NULL;
    }
    // optional: mat->n_in=0; mat->n_out=0;
}

// 2) n_in / n_out
size_t matrix_linear_n_in(const MatrixLinear *mat){
    return mat->n_in;
}
size_t matrix_linear_n_out(const MatrixLinear *mat){
    return mat->n_out;
}

// 3) apply: 
// mat is row-major: mat->entries is size n_out rows × n_in columns
// For row r in [0..n_out):
//     out[r]= ∑( col c in [0..n_in): mat->entries[r*n_in + c] * input[c] )
void matrix_linear_apply(const MatrixLinear *mat,
                         const F128 *input, size_t in_len,
                         F128 *output, size_t out_len)
{
    // dimension checks
    assert(in_len  == mat->n_in);
    assert(out_len == mat->n_out);

    // set output to zero
    for(size_t r=0; r<mat->n_out; r++){
        output[r]= f128_zero();
    }

    // row-major => each row is chunk of length n_in
    // row r => mat->entries[r*n_in + 0..n_in-1]
    for(size_t r=0; r< mat->n_out; r++){
        // pointer to row
        const F128* row_ptr= &mat->entries[r* mat->n_in];
        // out_elem => output[r]
        for(size_t c=0; c<mat->n_in; c++){
            // out[r] += row[r*n_in + c]* input[c]
            F128 partial= f128_mul(row_ptr[c], input[c]);
            output[r] = f128_add(output[r], partial);
        }
    }
}

// 4) apply_transposed:
// mat^T => dimension n_in × n_out
// so input has len= n_out, output= n_in
// out[c]= ∑( row r in [0..n_out): mat->entries[r*n_in + c]* input[r] )
void matrix_linear_apply_transposed(const MatrixLinear *mat,
                                    const F128 *input, size_t in_len,
                                    F128 *output, size_t out_len)
{
    assert(in_len  == mat->n_out);
    assert(out_len == mat->n_in);

    // zero output
    for(size_t c=0; c< mat->n_in; c++){
        output[c]= f128_zero();
    }

    // For each row r => row_ptr => mat->entries[r*n_in.. r*n_in + n_in -1]
    // out[c] += row_ptr[c] * input[r]
    for(size_t r=0; r< mat->n_out; r++){
        // input[r]
        const F128 in_val= input[r];
        const F128* row_ptr= &mat->entries[r* mat->n_in];
        for(size_t c=0; c< mat->n_in; c++){
            F128 partial= f128_mul(row_ptr[c], in_val);
            output[c]   = f128_add(output[c], partial);
        }
    }
}