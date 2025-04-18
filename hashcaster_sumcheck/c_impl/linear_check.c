#include "linear_check.h"

#include <stdlib.h>
#include <assert.h>
#include "field.h"
#include "mle_poly.h"
#include "matrix_linear.h"
#include "univariate_poly.h"

// Constructor equivalent for LinCheckBuilder:
LinCheckBuilder* lincheck_builder_new(
    MLE_POLY_SEQUENCE *polys,
    const Points *points,
    const MatrixLinear *matrix,
    size_t num_active_vars,
    const F128 *initial_claims,
    size_t N,
    size_t M
) {
    // Validate matrix dimensions
    assert(matrix->n_in == N * (1 << num_active_vars) && "Invalid matrix input dimension");
    assert(matrix->n_out == M * (1 << num_active_vars) && "Invalid matrix output dimension");

    // Validate num_vars
    size_t num_vars = points->len;
    assert(num_vars >= num_active_vars && "Number of variables must be >= active variables");

    // Validate polynomial sizes
    for (size_t i = 0; i < N; i++) {
        assert(polys->mle_poly[i].len == (1ULL << num_vars) && "Polynomial size mismatch");
    }

    LinCheckBuilder *builder = (LinCheckBuilder*)malloc(sizeof(LinCheckBuilder));
    builder->matrix = matrix;
    builder->polys = polys;
    builder->points = points;
    builder->num_vars = num_vars;
    builder->num_active_vars = num_active_vars;
    builder->initial_claims = (F128*)malloc(M * sizeof(F128));
    for (size_t i = 0; i < M; i++) {
        builder->initial_claims[i] = initial_claims[i];
    }
    builder->N = N;
    builder->M = M;

    return builder;
}

// Builds the ProdCheck object from LinCheckBuilder with gamma:
ProdCheck* lincheck_builder_build(LinCheckBuilder *builder, const F128 *gamma) {
    size_t chunk_size = 1ULL << builder->num_active_vars;

    // Split points into active and dormant
    Points* pt_active, *pt_dormant;
    points_split_at(builder->points, builder->num_active_vars, &pt_active, &pt_dormant);

    // eq_poly for dormant points
    MLE_POLY* eq_dormant = points_to_eq_poly(pt_dormant);

    // Restrict input polynomials
    MLE_POLY_SEQUENCE* p_polys = mle_sequence_new(builder->N, chunk_size, f128_zero());

    for (size_t i = 0; i < builder->N; i++) {
        for (size_t j = 0; j < eq_dormant->len; j++) {
            for (size_t k = 0; k < chunk_size; k++) {
                size_t idx = j * chunk_size + k;
                F128 prod = f128_mul(eq_dormant->coeffs[j], builder->polys->mle_poly[i].coeffs[idx]);
                p_polys->mle_poly[i].coeffs[k] = f128_add(p_polys->mle_poly[i].coeffs[k], prod);
            }
        }
    }

    // Compute gammas: [1, gamma, gamma^2, ..., gamma^(M-1)]
    F128 *gammas = malloc(builder->M * sizeof(F128));
    gammas[0] = f128_one();
    for (size_t i = 1; i < builder->M; i++) {
        gammas[i] = f128_mul(gammas[i - 1], *gamma);
    }

    // eq_poly for active points
    MLE_POLY* eq_active = points_to_eq_poly(pt_active);

    // Combine gamma powers and active equality polynomial
    size_t gamma_eqs_len = builder->M * chunk_size;
    F128 *gamma_eqs = malloc(gamma_eqs_len * sizeof(F128));
    for (size_t m = 0; m < builder->M; m++) {
        for (size_t i = 0; i < chunk_size; i++) {
            gamma_eqs[m * chunk_size + i] = f128_mul(gammas[m], eq_active->coeffs[i]);
        }
    }

    // Apply transposed matrix to compute q
    MLE_POLY* q = mle_poly_new_zeros(builder->N * chunk_size);
    matrix_linear_apply_transposed(builder->matrix, gamma_eqs, gamma_eqs_len , q->coeffs, q->len);

    // Split q into N polynomials
    MLE_POLY_SEQUENCE* q_polys = mle_sequence_new(builder->N, chunk_size, f128_zero());
    for (size_t i = 0; i < builder->N; i++) {
        memcpy(q_polys->mle_poly[i].coeffs, &q->coeffs[i * chunk_size], chunk_size * sizeof(F128));
    }

    // Verify all chunks processed
    assert(q->len == builder->N * chunk_size);

    UnivariatePolynomial* uni_poly = malloc(sizeof(UnivariatePolynomial));
    uni_poly->elems = malloc(builder->M * sizeof(F128));
    uni_poly->len = builder->M;
    for (size_t i = 0; i < builder->M; i++) {
        uni_poly->elems[i] = builder->initial_claims[i];
    }
    // Evaluate the polynomial at the point gamma
    F128 claim = univariate_polynomial_evaluate_at(uni_poly, *gamma);

    // Create ProdCheck instance
    ProdCheck* prodcheck = prodcheck_new(p_polys->mle_poly, q_polys->mle_poly, builder->N, claim, 1);

    // Cleanup
    mle_poly_free(eq_dormant);
    mle_poly_free(eq_active);
    mle_poly_free(q);
    points_free((Points*)uni_poly);
    free(gammas);
    free(gamma_eqs);
    points_free(pt_active);
    points_free(pt_dormant);

    return prodcheck;
}

// Free LinCheckBuilder resources
void lincheck_builder_free(LinCheckBuilder *builder) {
    free(builder->initial_claims);
    free(builder);
}

