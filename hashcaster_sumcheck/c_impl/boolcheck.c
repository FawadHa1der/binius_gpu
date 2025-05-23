#include "boolcheck.h"
#include <math.h>
#include "mle_poly.h"
#include "debug_utils.h"
#include "types.h"

BoolCheckBuilder* bool_check_builder_new(
    size_t C_INITIAL_ROUNDS,
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

Points* extend_n_tables(
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

    //F128* result = malloc(sizeof(F128) * pow3 * pow2);
    Points* result = points_init(pow3 * pow2, f128_zero());

    assert(result != NULL);
    for (size_t chunk_id = 0; chunk_id < pow2; ++chunk_id) {
        F128 tables_ext[pow3_adj][N];

        size_t base_tab_offset = chunk_id * base_stride;
        F128* result_chunk = &result->elems[chunk_id * pow3];

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


// Implements the logic of SumcheckBuilder::build for BoolCheckBuilder.
BoolCheck* boolcheck_new(
    BoolCheckBuilder* builder,
    const F128 gamma 
) {
    // 1. Compute gammas from gamma (using initial_claim as gamma here)
    // For test, we use initial_claim as the gamma input.
    // compute_gammas_folding fills builder->gammas with folding challenges.
    builder->gammas = compute_gammas_folding(gamma, builder->algebraic_operations->output_size);

    // 2. Ensure all polynomials have same length = 2^num_vars
    size_t expected_len = 1ULL << builder->points->len;
    for (size_t i = 0; i < builder->polys->len; ++i) {
        assert(expected_len == builder->polys->mle_poly[i].len);
    }

    // 3. Generate bit and trit mappings
    uint16_t* bit_mapping = NULL;
    uint16_t* trit_mapping = NULL;
    size_t bit_len = 0, trit_len = 0;
    compute_trit_mappings(builder->c_initial_rounds, &bit_mapping, &bit_len, &trit_mapping, &trit_len);

    // 4. Generate eq_sequence (offset by 1)
    Points* points_off_by_one = points_init(builder->points->len - 1, f128_zero());
    for (size_t i = 0; i < points_off_by_one->len; ++i) {
        points_off_by_one->elems[i] = builder->points->elems[i + 1];
    } // nbot sure why its off by one?

    MLE_POLY_SEQUENCE* eq_sequence = points_to_eq_poly_sequence(
        points_off_by_one);

    // 5. Compute compressed claim using FixedUnivariatePolynomial
    // Use builder->claims and builder->gammas->elems[0] as evaluation point
    F128 claim =   univariate_polynomial_evaluate_at (builder->claims, gamma);

    // 6. Extend tables
    Points* extended_table = extend_n_tables(
        builder->polys, builder->polys->len, builder->points->len, builder->c_initial_rounds,
        trit_mapping, builder, linear_compressed, quadratic_compressed
    );

    // 7. Allocate and return BoolCheck
    BoolCheck* result = (BoolCheck*) malloc(sizeof(BoolCheck));
    result->bit_mapping = bit_mapping;
    result->bit_mapping_len = bit_len;
    // result->bit_len = bit_len;
    result->eq_sequence = eq_sequence;
    result->claim = claim;
    result->extended_table =  extended_table;
    result->polys = builder->polys;
    result->points = builder->points;
    result->gammas = builder->gammas;
    result->algebraic_operations = builder->algebraic_operations;
    result->algebraic_functions = builder->algebraic_functions;
    result->poly_coords = points_init(0, f128_zero());
    result->challenges = points_init(0, f128_zero());
    result->round_polys = malloc(sizeof(CompressedPoly) * builder->points->len);
    if (result->round_polys == NULL) {
        fprintf(stderr, "Failed to allocate memory for round_polys\n");
        free(result);
        return NULL;
    }
    result->round_polys_len = 0;
    result->num_vars = builder->points->len;
    result->c_param = builder->c_initial_rounds;
    // result->bit_mapping = bit_mapping;

    return result;
}

CompressedPoly* boolcheck_round_polynomial(BoolCheck* bc) {
    size_t round = bc->challenges->len;
    size_t num_vars = bc->num_vars;
    assert(round < num_vars);

    if (bc->round_polys_len > round) {
        return &(bc->round_polys[round]);  // Already cached
    }

    Points* poly_deg2 = points_init(3, f128_zero());

    if (round <= bc->c_param) {
        // === Phase 1 computation using extended table ===
        size_t dims_tail = num_vars - bc->c_param - 1;
        size_t pow3 = pow(3, bc->c_param - round); // needs helper
        for (size_t i = 0; i < (1UL << dims_tail); i++) {
            size_t base_index = i << (bc->c_param - round);
            size_t base_offset = 3 * (i * pow3);
            for (size_t j = 0; j < (1UL << (bc->c_param - round)); j++) {
                size_t offset = base_offset + 3 * bc->bit_mapping[j];
                // F128 multiplier = bc->eq_sequence[num_vars - round - 1].data[base_index + j];
                F128 multiplier = bc->eq_sequence->mle_poly[num_vars - round - 1].coeffs[base_index + j];

                for (int t = 0; t < 3; t++) {
                    poly_deg2->elems[t] = f128_add(poly_deg2->elems[t],
                        f128_mul(bc->extended_table->elems[offset + t], multiplier));
                }
            }
        }
    } else {
        // === Phase 2 computation using restrict ===
        size_t half = 1UL << (num_vars - round - 1);
        for (size_t i = 0; i < half; i++) {
            F128 acc[3] = {f128_zero(), f128_zero(), f128_zero()};
            // F128* algebraic = compute_algebraic(bc, i, 1UL << (num_vars - bc->c_param - 1));
            algebraic_compressed(
                bc->gammas,
                bc->algebraic_operations,
                bc->algebraic_functions,
                bc->poly_coords,
                i,
                1UL << (num_vars - bc->c_param - 1),
                acc
            );
            F128 multiplier = bc->eq_sequence->mle_poly[bc->points->len - round - 1].coeffs[i];
            for (int k = 0; k < 3; k++) {
                poly_deg2->elems[k] = f128_add(poly_deg2->elems[k], f128_mul(acc[k], multiplier));
            }
        }
    }

    // === Compute final compressed polynomial ===
    F128 eq_y = eq_eval(bc->points->elems, bc->challenges->elems, round);
    
    // F128 eq_y = points_eq_eval_slice(&bc->challenges, bc->points, round);
    UnivariatePolynomial* univariate_poly_eval = from_evaluations_deg2(poly_deg2);
    
    // V(t) = W(t) * eq(r_<i>; q)
    for (int i = 0; i < 3; i++) {
        poly_deg2->elems[i] = f128_mul(univariate_poly_eval->elems[i], eq_y);
    }

    // F128 eq_t[2] = {
    //     f128_add(bc->points->elems[round], f128_one()),
    //     f128_one()
    // };
    UnivariatePolynomial* eq_t = points_init(2, f128_zero());
    eq_t->elems[0] = f128_add(bc->points->elems[round], f128_one());
    eq_t->elems[1] = f128_one();

    UnivariatePolynomial* result = multiply_degree2_by_degree1(poly_deg2, eq_t); // returns coeffs[3]

    CompressedPoly* cp = compress_poly(result);
    // F128 computed_claim = evaluate_univariate(&cp, bc->claim);
    // F128 computed_claim = univariate_polynomial_evaluate_at(result, bc->claim);
    assert(f128_eq(cp->sum, bc->claim));

    bc->round_polys[round] = *cp; // TODO clean this up, not sure this is the best approach
    bc->round_polys_len++;
    return cp;
}

void boolcheck_bind(BoolCheck* bc, const F128* r) {
    size_t round = bc->challenges->len;
    size_t num_vars = bc->num_vars;
    assert(round < num_vars);

    // 1. Compute round polynomial, decompress, evaluate at *r, update claim
    CompressedPoly* round_poly = boolcheck_round_polynomial(bc);
    UnivariatePolynomial* round_poly_coeffs = uncompress_poly(round_poly, bc->claim );
    bc->claim = univariate_polynomial_evaluate_at(round_poly_coeffs, *r);
    // Free temp poly
    free(round_poly_coeffs);

    // bc->challenges->elems[round] = *r;

    // // 2. Push *r to bc->challenges
    points_push(bc->challenges, *r);

    // 3. If round <= bc->c_param, phase 1
    if (round <= bc->c_param) {
        // Compute r^2
        F128 r2 = f128_mul(*r, *r);
        // Number of chunks = bc->extended_table->len / 3
        size_t n_chunks = bc->extended_table->len / 3;
        for (size_t i = 0; i < n_chunks; ++i) {
            F128* chunk = &bc->extended_table->elems[i * 3];
            // chunk[0] = chunk[0] + (chunk[0] + chunk[1] + chunk[2]) * r + chunk[2] * r^2
            F128 sum012 = f128_add(chunk[0], f128_add(chunk[1], chunk[2]));
            F128 t1 = f128_add(chunk[0], f128_mul(sum012, *r));
            chunk[0] = f128_add(t1, f128_mul(chunk[2], r2));
            // chunk[1], chunk[2] are not needed anymore
        }
        // Compact vector: keep every third value (i.e., chunk[0] of each chunk)
        for (size_t i = 0; i < n_chunks; ++i) {
            bc->extended_table->elems[i] = bc->extended_table->elems[i * 3];
        }
        bc->extended_table->len = n_chunks;
    } else {
        // Phase 2
        size_t half = 1UL << (num_vars - round - 1);
        size_t chunk_size = 1UL << (num_vars - bc->c_param - 1);
        size_t n_chunks = bc->poly_coords->len / chunk_size;
        for (size_t chunk_idx = 0; chunk_idx < n_chunks; ++chunk_idx) {
            F128* chunk = &bc->poly_coords->elems[chunk_idx * chunk_size];
            for (size_t j = 0; j < half; ++j) {
                // chunk[j] = chunk[2*j] + (chunk[2*j + 1] + chunk[2*j]) * r
                F128 sum = f128_add(chunk[2*j+1], chunk[2*j]);
                chunk[j] = f128_add(chunk[2*j], f128_mul(sum, *r));
            }
        }
        // bc->poly_coords->len remains the same
    }

    // 4. If bc->challenges->len == bc->c_param + 1, clear extended_table and restrict
    if (bc->challenges->len == bc->c_param + 1) {
        // Clear extended_table
        bc->extended_table->len = 0;
        // Compute restrict and store in poly_coords
        // restrict(polys, challenges, num_vars)
        // polys: bc->polys, n = bc->n_polys
        // challenges: bc->challenges
        // Output length: n_polys * 128 * base_index
        
        Evaluations* restricted_evals = restrict_polynomials(
            bc->polys->mle_poly,
            bc->polys->len,
            bc->challenges,
            num_vars
        );
        if (bc->poly_coords) points_free(bc->poly_coords);
        bc->poly_coords = restricted_evals;
    }
}

BoolCheckOutput* boolcheck_finish(BoolCheck* bc) {
    // 1. Assert that bc->challenges->len == bc->num_vars
    assert(bc->challenges->len == bc->num_vars);
    // 2. Allocate BoolCheckOutput* out, compute base_index
    size_t base_index = 1UL << (bc->num_vars - bc->c_param - 1);
    
    size_t output_size = bc->algebraic_operations->input_size; // or number of polynomials
    // 3. Populate frob_evals_array of length 128 * output_size
    size_t frob_evals_len = 128 * output_size;
    Evaluations* frob_evals = points_init(frob_evals_len, f128_zero());
    // F128* frob_evals_array = (F128*)malloc(frob_evals_len * sizeof(F128));
    for (size_t idx = 0; idx < frob_evals_len; ++idx) {
        size_t poly_idx = idx / 128;
        size_t frob_idx = idx % 128;
        // poly_coords: [n_polys * 128 * base_index]
        frob_evals->elems[idx] = bc->poly_coords->elems[(poly_idx * 128 + frob_idx) * base_index];
    }

    twist_evals(frob_evals);
    BoolCheckOutput* out = malloc(sizeof(BoolCheckOutput));
    if (out == NULL) {
        fprintf(stderr, "Failed to allocate memory for BoolCheckOutput\n");
        points_free(frob_evals);
        return NULL;
    }
    out->frob_evals = frob_evals;
    out->round_polys = bc->round_polys;
    if (out->round_polys == NULL) {
        fprintf(stderr, "Failed to allocate memory for round_polys\n");
        points_free(frob_evals);
        free(out);
        return NULL;
    }
    out->round_polys_len = bc->round_polys_len;
    return out;
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

void boolcheck_free(BoolCheck* bc) {
    if (bc == NULL) {
        return;
    }
    free(bc->bit_mapping);
    mle_sequence_free(bc->eq_sequence);
    points_free(bc->extended_table);
    // points_free(bc->points);
    points_free(bc->gammas);
    points_free(bc->poly_coords);
    points_free(bc->challenges);
    free(bc->round_polys);
    free(bc);
}
void boolcheck_output_free(BoolCheckOutput* out) {
    if (out == NULL) {
        return;
    }
    points_free(out->frob_evals);
    free(out);
}