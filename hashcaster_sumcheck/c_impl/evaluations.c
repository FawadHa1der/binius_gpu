#include "evaluations.h"

// twists them in place
void twist_evals(Evaluations *evals)
{
    const size_t chunk_len = 128;
    // 1) Check we can chunk the data in 128-element groups
    if(evals->len % chunk_len != 0){
        fprintf(stderr, "Evaluations must be a multiple of 128!\n");
        return;
    }

    size_t total = evals->len;
    size_t chunkCount = total / chunk_len;
    F128* array = evals->elems;

    for(size_t c=0; c<chunkCount; c++){
        // pointer to the chunk
        F128* chunk = &array[c * chunk_len];

        // We'll build twisted_evals[128] via "array::from_fn(|_| {...})" logic
        // meaning for i in [0..128), we do:
        //   chunk: square each element
        //   sum= (0..128).fold(ZERO, |acc, j| basis(j)*chunk[j] + acc)
        // store sum in twisted[i]

        F128 twisted[chunk_len];

        for(size_t i=0; i<chunk_len; i++){
            // (a) square all 128 elements of chunk
            for(size_t k=0; k<chunk_len; k++){
                chunk[k] = f128_mul(chunk[k], chunk[k]);
            }

            // (b) do sum_{j=0..127} basis(j)* chunk[j]
            F128 sum = f128_zero();
            for(size_t j=0; j<chunk_len; j++){
                F128 b = f128_basis(j);
                F128 prod = f128_mul(b, chunk[j]);
                sum = f128_add(sum, prod);
            }
            twisted[i] = sum;
        }

        // 3) reverse twisted
        for(size_t i=0; i<64; i++){
            F128 tmp = twisted[i];
            twisted[i] = twisted[127-i];
            twisted[127-i] = tmp;
        }

        // 4) copy twisted back to chunk
        for(size_t i=0; i<chunk_len; i++){
            chunk[i] = twisted[i];
        }
    }
}


void untwist_evals(Evaluations *evals)
{
    const size_t chunk_len = 128;

    // Check multiple of 128
    if(evals->len % chunk_len != 0){
        fprintf(stderr, "Evaluations must be a multiple of 128.\n");
        return;
    }

    size_t chunkCount = evals->len / chunk_len;
    F128 * array = evals->elems;

    for(size_t c=0; c<chunkCount; c++){
        // chunk pointer
        F128 *chunk = &array[c*chunk_len];

        // Step 1: Frobenius transformation x -> x^(2^i)
        // chunk[i] = chunk[i].frobenius(i) for i in [0..127]
        for(size_t i=0; i<chunk_len; i++){
            chunk[i] = f128_frob(chunk[i], i);
        }

        F128 untwisted[chunk_len];
        for(size_t i=0; i<chunk_len; i++){
            // can do pre-computations of pi_calc(i, chunk) like in hashcaster explorations
            untwisted[i] = pi_calc(i, chunk);
        }

        // Step 3: Copy untwisted back
        // chunk.copy_from_slice(&untwisted_chunk);
        for(size_t i=0; i<chunk_len; i++){
            chunk[i] = untwisted[i];
        }
    }
}




