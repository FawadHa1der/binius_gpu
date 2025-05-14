#include "boolcheck.h"
#include <math.h>
#include "mle_poly.h"
#include "debug_utils.h"

BoolCheckBuilder* bool_check_builder_new(
    int C_INITIAL_ROUNDS,
    const Points* points,
    Points* claims,
    MLE_POLY_SEQUENCE* polys,
    const Algebraic_Params* algebraic_params,
    const Algebraic_Functions* algebraic_functions
) {
    if (C_INITIAL_ROUNDS >= points->len) {
        fprintf(stderr, "Phase switch parameter C too large\n");
        return NULL;
    }
// N is poly sequence length 
// M is is the clainms,gamma length
    BoolCheckBuilder* builder = malloc(sizeof(BoolCheckBuilder));
    if (builder == NULL) {
        fprintf(stderr, "Failed to allocate memory for BoolCheckBuilder\n");
        return NULL;
    }
    builder->c_initial_rounds = C_INITIAL_ROUNDS;
    builder->points = points;
    builder->gammas = points_init(claims->len, f128_zero());
    builder->claims = points_copy(claims);
    builder->polys = mle_sequence_copy(polys);
    builder->algebraic_operations = algebraic_params;
    builder->algebraic_functions = algebraic_functions;

    // builder->algebraic_operations = algebraic_operations;

    // memcpy(builder->claims->elems, claims->elems, sizeof(F128) * claims->len);
    // for (size_t i = 0; i < builder->claims->len; ++i)
    //     builder->gammas->elems[i] = f128_zero();

    return builder;
}

// Computes bit_mapping and trit_mapping as described
void compute_trit_mappings(size_t c_rounds, uint16_t** bit_mapping_out, size_t* bit_len,
                          uint16_t** trit_mapping_out, size_t* trit_len) {
    size_t pow3 = 1;
    for (size_t i = 0; i < c_rounds + 1; ++i)
        pow3 *= 3;

    uint16_t* trit_mapping = (uint16_t*)malloc(pow3 * sizeof(uint16_t));
    uint16_t* bit_mapping = (uint16_t*)malloc((1 << (c_rounds + 1)) * sizeof(uint16_t));
    size_t bit_index = 0;

    for (size_t i = 0; i < pow3; ++i) {
        size_t current = i;
        uint16_t bin_value = 0;
        int msd = -1;

        // Convert to base-3 and extract binary equivalent + most significant 2
        for (size_t idx = 0; idx <= c_rounds; ++idx) {
            uint16_t digit = current % 3;
            current /= 3;

            if (digit == 2) {
                msd = idx;
            } else {
                bin_value |= (digit & 1) << idx;
            }
        }

        if (msd != -1) {
            trit_mapping[i] = (uint16_t)pow(3, msd);
        } else {
            trit_mapping[i] = bin_value << 1;
            bit_mapping[bit_index++] = (uint16_t)i;
        }
    }

    *bit_mapping_out = (uint16_t*)realloc(bit_mapping, bit_index * sizeof(uint16_t));
    *bit_len = bit_index;
    *trit_mapping_out = trit_mapping;
    *trit_len = pow3;
}


//N input
//M output

F128 linear_compressed(
    const Points* gammas,
    const Algebraic_Params* params,
    const Algebraic_Functions* alg_funcs,
    const Points* arg
) {
    // F128 lin[M];
    Points* lin = points_init(params->output_size, f128_zero());
    alg_funcs->linear(params, arg, lin);  // linear() fills lin[M]

    F128 acc = f128_zero();
    for (size_t i = 0; i < params->output_size; ++i) {
        acc = f128_add(acc, f128_mul(lin->elems[i], gammas->elems[i]));
    }
    return acc;
}

F128 quadratic_compressed(
    const Points* gammas,
    const Algebraic_Params* params,
    const Algebraic_Functions* alg_funcs,
    const Points* arg
) {
    //F128 quad[M];
    Points* quad = points_init(params->output_size, f128_zero());

    
    alg_funcs->quadratic(params, arg, quad);  // quadratic() fills quad[params->output_size]

    F128 acc = f128_zero();
    for (size_t i = 0; i < params->output_size; ++i) {
        acc = f128_add(acc, f128_mul(quad->elems[i], gammas->elems[i]));
    }
    return acc;
}

void algebraic_compressed(
    const Points* gammas,
    const Algebraic_Params* params,
    const Algebraic_Functions* alg_funcs,
    const Points* data,
    size_t idx_a,
    size_t offset,
    F128 result[3]
) {
    F128 *alg[3];
    for (size_t i = 0; i < 3; ++i) {
        alg[i] = (F128*)malloc(sizeof(F128) * params->output_size);
        for (size_t j = 0; j < params->output_size; ++j) {
            alg[i][j] = f128_zero();
        }
    }
    alg_funcs->algebraic(params, data, idx_a, offset, alg);  // algebraic() fills alg[3][M]

    for (size_t j = 0; j < 3; ++j) {
        result[j] = f128_zero();
        for (size_t i = 0; i < params->output_size; ++i) {
            result[j] = f128_add(result[j], f128_mul(alg[j][i], gammas->elems[i]));
        }
    }
}

F128* extend_n_tables(
    const MLE_POLY_SEQUENCE* polys,  // array of N pointers to F128 arrays
    size_t N,
    size_t dims,
    size_t C,
    const uint16_t* trit_mapping,
    BoolCheckBuilder* builder,
    F128 (*linear_compressed)(
        const Points* gammas,
        const Algebraic_Params* params,
        const Algebraic_Functions* alg_funcs,
        const Points* arg
    ),
    F128 (*quadratic_compressed)(
        const Points* gammas,
        const Algebraic_Params* params,
        const Algebraic_Functions* alg_funcs,
        const Points* arg
    )
) {
    assert(C < dims);

    size_t pow3 = 1;
    for (size_t i = 0; i <= C; ++i) pow3 *= 3;

    size_t pow3_adj = 2 * pow3 / 3;
    size_t pow2 = 1 << (dims - C - 1);
    size_t base_stride = 1 << (C + 1);

    F128* result = malloc(sizeof(F128) * pow3 * pow2);
    assert(result != NULL);

    for (size_t chunk_id = 0; chunk_id < pow2; ++chunk_id) {
        F128 tables_ext[pow3_adj][N];

        size_t base_tab_offset = chunk_id * base_stride;
        F128* result_chunk = &result[chunk_id * pow3];

        for (size_t j = 0; j < pow3; ++j) {
            size_t offset = trit_mapping[j];

            if (j < pow3_adj) {
                if (offset % 2 == 0) {
                    size_t idx = base_tab_offset + (offset >> 1);
                    for (size_t z = 0; z < N; ++z) {
                        tables_ext[j][z] = polys->mle_poly[z].coeffs[idx];
                    }

                    Points  tables_ext_points_shallow_copy ;
                    tables_ext_points_shallow_copy.elems = tables_ext[j];
                    tables_ext_points_shallow_copy.len = N;

                    result_chunk[j] = f128_add( quadratic_compressed(builder->gammas, builder->algebraic_operations,
                        builder->algebraic_functions, &tables_ext_points_shallow_copy) ,
                                      linear_compressed(builder->gammas, builder->algebraic_operations,
                                                        builder->algebraic_functions, &tables_ext_points_shallow_copy));
                    // print result_chunk[j]);
                } else {
                    size_t j1 = j - offset;
                    size_t j2 = j - 2 * offset;
                    for (size_t z = 0; z < N; ++z) {
                        tables_ext[j][z] = f128_add( tables_ext[j1][z] , tables_ext[j2][z]);
                    }

                    Points  tables_ext_points_shallow_copy ;
                    tables_ext_points_shallow_copy.elems = tables_ext[j];
                    tables_ext_points_shallow_copy.len = N;

                    result_chunk[j] = quadratic_compressed(builder->gammas, builder->algebraic_operations,
                        builder->algebraic_functions, &tables_ext_points_shallow_copy);
                }
            } else {
                size_t j1 = j - offset;
                size_t j2 = j - 2 * offset;
                F128 args[N];
                for (size_t z = 0; z < N; ++z) {
                    args[z] = f128_add( tables_ext[j1][z] , tables_ext[j2][z]);
                }
                Points  args_shallow_copy ;
                args_shallow_copy.elems = args;
                args_shallow_copy.len = N;

                result_chunk[j] = quadratic_compressed(builder->gammas, builder->algebraic_operations,
                    builder->algebraic_functions, &args_shallow_copy);
            }
            // printf("result_chunk[%zu] = ", j);
            // f128_print (  "is" , result_chunk[j]);

        }
    }

    return result;
}

// Free the BoolCheckBuilder
void bool_check_builder_free(BoolCheckBuilder* builder) {
    if (builder == NULL) {
        return;
    }
    // points_free(builder->points); TODO better cleanup
    points_free(builder->gammas);
    points_free(builder->claims);
    mle_sequence_free(builder->polys);
    free(builder);
}