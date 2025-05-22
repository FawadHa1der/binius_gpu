#include "and_package.h"


void and_package_algebraic(
    const Algebraic_Params* params,
    const Points* data,
    size_t idx_a,
    size_t offset,
    F128 *ret[3] // TODO a better data structure for a 2D array ?
) {
    // Initialize indices
    idx_a *= 2;
    size_t idx_b = idx_a + offset * 128;

    // Initialize result to zero
    for (size_t j = 0; j < 3; j++) {
        for (size_t k = 0; k < params->output_size; k++) {
            ret[j][k] = f128_zero();
        }
    }

    // Iterate over 128 basis elements
    for (size_t i = 0; i < 128; i++) {
        F128 basis = f128_basis(i);  // assume this returns the i-th basis element

        F128 a = data->elems[idx_a];
        F128 b = data->elems[idx_b];

        F128 a_next = data->elems[idx_a + 1];  // assume bounds are managed externally
        F128 b_next = data->elems[idx_b + 1];
        if (idx_b + 1 >= data->len) {
            // a_next = f128_zero();
            b_next = f128_zero();
        }

        if (idx_a + 1 >= data->len) {
            // a_next = f128_zero();
            a_next = f128_zero();
        }   
        

        // Σ (ϕ_i * a * b)
        F128 t0 = f128_mul(a, b);
        ret[0][0] = f128_add(ret[0][0], f128_mul(basis, t0));

        // Σ (ϕ_i * a_next * b_next)
        F128 t1 = f128_mul(a_next, b_next);
        ret[1][0] = f128_add(ret[1][0], f128_mul(basis, t1));

        // Σ (ϕ_i * (a + a_next) * (b + b_next))
        F128 t2 = f128_mul(f128_add(a, a_next), f128_add(b, b_next));
        ret[2][0] = f128_add(ret[2][0], f128_mul(basis, t2));

        idx_a += offset;
        idx_b += offset;
    }
}

void and_package_linear(
    const Algebraic_Params* params,
    const Points* data,
    Points* out
) {
    for (size_t i = 0; i < params->output_size; i++) {
        out->elems[i] = f128_zero();
    }
}

void and_package_quadratic(
    const Algebraic_Params* params,
    const Points* data,
    Points* out
) {
    for (size_t i = 0; i < params->output_size; i++) {
        out->elems[i] = f128_bitand(data->elems[0], data->elems[1]);
    }
}