#include "multi_claim.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "efficient_matrix.h"

// // Helper function to compute gamma powers
// static void compute_gamma_pows(F128 gamma, size_t n, F128* out) {
//     //out[0] = 1.0; // Assuming F128 can be initialized like this
//     out[0] = f128_one(); // Initialize first element to 1.0
//     for (size_t i = 1; i < n; i++) {
//         //out[i] = out[i - 1] * gamma;
//         out[i] = f128_mul(out[i - 1], gamma);
//     }
// }

MulticlaimBuilder* multiclaim_builder_new(
    MLE_POLY_SEQUENCE* polys,
    Points* points,
    Evaluations* openings
) {
    assert(polys != NULL);
    assert(points != NULL);
    assert(openings != NULL);
    assert(polys->len > 0);

    // Check all polys have length == 2^points->len
    size_t expected_len = 1 << points->len;
    for (size_t i = 0; i < polys->len; i++) {
        size_t poly_len =  polys->mle_poly[i].len;
        assert(poly_len == expected_len);
    }

    MulticlaimBuilder* builder = (MulticlaimBuilder*)malloc(sizeof(MulticlaimBuilder));
    if (!builder) {
        // Allocation failed
        return NULL;
    }

    builder->polys = polys;
    builder->points = points;
    builder->openings = openings;

    return builder;
}

MultiClaim* multiclaim_builder_build(
    MulticlaimBuilder* builder,
    F128 gamma
) {
    assert(builder != NULL);
    Points* points = builder->points;
    size_t num_points = points->len;
    size_t poly_len = 1 << num_points;

    // Allocate gamma_pows array: size 128 * num_polys
    size_t gamma_pows_len = 128 * builder->polys->len;
    Points* gamma_pows =  compute_gammas_folding(gamma, gamma_pows_len);
    if (!gamma_pows) {
        return NULL;
    }

    // // For each polynomial, compute gamma powers for 128 elements
    // for (size_t poly_i = 0; poly_i < builder->polys->len; poly_i++) {
    //     F128 base_gamma = gamma;
    //     // Compute powers for this polynomial segment
    //     compute_gamma_pows(base_gamma, 128, &gamma_pows[poly_i * 128]);
    // }

    // Allocate poly_composite array of size poly_len
    F128* poly_composite = (F128*)malloc(sizeof(F128) * poly_len);
    if (!poly_composite) {
        points_free(gamma_pows);
        return NULL;
    }
    // Initialize poly_composite to zero
    memset(poly_composite, 0, sizeof(F128) * poly_len);

    // Folding poly_composite = sum_i gamma^i * polys[i]
    // For each poly element j in [0..poly_len)
    for (size_t j = 0; j < poly_len; j++) {
        F128 acc = f128_zero();
        for (size_t i = 0; i < builder->polys->len; i++) {
            // F128 poly_val = mlp_get(&builder->polys->mle_poly[i], j);
            F128 poly_val = builder->polys->mle_poly[i].coeffs[j] ;
            F128 gamma_pow = gamma_pows->elems[i * 128]; // Use first gamma power for folding
            //acc += gamma_pow * poly_val;
            acc = f128_add(acc, f128_mul(gamma_pow, poly_val));
        }
        poly_composite[j] = acc;
    }

    // Allocate folded_openings array size 128
    F128* folded_openings = (F128*)malloc(sizeof(F128) * 128);
    if (!folded_openings) {
        points_free(gamma_pows);
        free(poly_composite);
        return NULL;
    }
    memset(folded_openings, 0, sizeof(F128) * 128);

    // Folding openings: folded_openings[j] = sum_i gamma^i * Evaluations_get(openings, i, j)
    for (size_t i = 0; i < 128; i++) {
        F128 acc = f128_zero();
        for (size_t j = 0; j < builder->polys->len; j++) {
            F128 open_val = builder->openings->elems[i + j * 128];
            F128 gamma_pow = gamma_pows->elems[128 * j];

            acc = f128_add(acc, f128_mul(gamma_pow, open_val));
        }
        folded_openings[i] = acc;
    }
    MLE_POLY* ml_poly_composite = malloc(sizeof(MLE_POLY));
    ml_poly_composite->len = poly_len;
    ml_poly_composite->coeffs = poly_composite;

    MultiClaim* claim = multi_claim_new(
        ml_poly_composite,
        points,
        folded_openings,
        gamma_pows,
        builder->polys
    );

    // Free temporary arrays if MultiClaim_new made copies, else keep ownership
    // Assuming MultiClaim_new takes ownership, do not free here
    return claim;
}

MultiClaim* multi_claim_new(
    MLE_POLY* poly,
    const Points *points,
    const F128 *openings,
    const Points *gamma_pows,
    MLE_POLY_SEQUENCE *polys
    // size_t N
) {
    MultiClaim *mc = malloc(sizeof(MultiClaim));
    mc->polys = polys;
    // mc->N = N;

    // Frobenius linear matrix
    EfficientMatrix *m = from_frobenius_inv_lc(gamma_pows->elems);

    // Equality poly
    MLE_POLY *eq = points_to_eq_poly(points);
    for (size_t i = 0; i < eq->len; i++) {
        eq->coeffs[i] = efficient_matrix_apply(m, eq->coeffs[i]);
    }

    // Initial claim: sum gamma_pows[i] * openings[i]
    F128 claim = f128_zero();
    for (size_t i = 0; i < 128; i++) {
        claim = f128_add(claim, f128_mul(gamma_pows->elems[i], openings[i]));
    }

    // Gamma^128 if available
    mc->gamma = gamma_pows->elems[128];

    // Wrap eq and poly into ProdCheck
    mc->object = prodcheck_new(poly, eq, 1, claim, 0);

    mle_poly_free(eq);
    efficient_matrix_free(m);

    return mc;
}

void multi_claim_bind(MultiClaim *mc, const F128 challenge, int challenge_index) {
    prodcheck_bind(mc->object, challenge, challenge_index);
}

Evaluations* multi_claim_finish(MultiClaim *mc) {
    Evaluations *out = points_init(mc->polys->len, f128_zero());
    F128 *evs = out->elems;

    // All evals except index 0
    for (size_t i = 1; i < mc->polys->len; i++) {
        evs[i] = mle_poly_evaluate_at(&mc->polys->mle_poly[i], mc->object->challenges);
    }

    // Construct univariate poly with these coeffs
    UnivariatePolynomial *uni = points_init(mc->polys->len, f128_zero());
    for (size_t i = 0; i < mc->polys->len; i++) {
        uni->elems[i] = evs[i];
    }

    // evs[0] = eval_at(gamma) + folded_poly_eval
    F128 eval_at_gamma = univariate_polynomial_evaluate_at(uni, mc->gamma);
    evs[0] = f128_add(eval_at_gamma, mc->object->p_polys[0].coeffs[0]);

    points_free(uni);
    return out;
}

void multi_claim_free(MultiClaim *mc) {
    prodcheck_free(mc->object);
    free(mc);
}

