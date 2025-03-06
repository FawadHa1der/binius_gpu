
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <arm_neon.h>  // If you need NEON intrinsics
#include "unity.h"

#include "../field.h"
#include "../matrix_utils.h"

/******************************************************************************
 * 1) test_apply_matrix
 *
 * Rust logic:
 *   - Build a 128-element "matrix" array of random F128.
 *   - Build a random "vec" (u128).
 *   - Convert "vec" to bits, each column of "matrix" to bits => compute a
 *     reference "expected_answer" bit array by the known formula.
 *   - Compare that to the output of matrix_apply().
 ******************************************************************************/
void test_apply_matrix(void)
{
    // 1) Create a random matrix of 128 columns
    F128 matrix_cols[128];
    for (int i = 0; i < 128; i++) {
        matrix_cols[i] = f128_rand();  // a stand-in for your "u128_rand(rng)"
    }

    // 2) Create a random vector
    F128 vec = f128_rand();

    // 3) Convert "vec" to bits
    bool vec_bits[128];
    binary128_to_bits(vec, vec_bits);

    // 4) Convert each matrix column to bits, store them in matrix_bits[i][j]
    bool matrix_bits[128][128];
    for (int i = 0; i < 128; i++) {
        bool tmp[128];
        binary128_to_bits(matrix_cols[i], tmp);
        for (int j = 0; j < 128; j++) {
            matrix_bits[i][j] = tmp[j];
        }
    }

    // 5) Build "expected_answer" as a 128-bit boolean array
    //    For each i in [0..128), for j in [0..128):
    //      expected_answer[j] ^= matrix_bits[i][j] && vec_bits[i]
    bool expected_answer[128];
    for (int j = 0; j < 128; j++) {
        expected_answer[j] = false;
    }
    for (int i = 0; i < 128; i++) {
        if (vec_bits[i]) {
            // XOR each bit j with matrix_bits[i][j]
            for (int j = 0; j < 128; j++) {
                // GF(2) addition => "xor" is the same as "!="
                expected_answer[j] = (expected_answer[j] != matrix_bits[i][j]);
            }
        }
    }

    // 6) Create a Matrix from "matrix_cols" and apply "vec"
    Matrix m = matrix_new(matrix_cols);
    F128 answer = matrix_apply(&m, vec);

    // 7) Convert "answer" to bits
    bool answer_bits[128];
    binary128_to_bits(answer, answer_bits);

    // 8) Compare each bit with expected_answer
    for (int j = 0; j < 128; j++) {
        if (answer_bits[j] != expected_answer[j]) {
            char msg[128];
            sprintf(msg, "Mismatch at bit %d in test_apply_matrix", j);
            TEST_FAIL_MESSAGE(msg);
        }
    }

    // If no mismatch, we pass
    TEST_PASS();
}

/******************************************************************************
 * 2) test_invert_matrix
 *
 * Rust logic:
 *   - Start with matrix = Matrix::diag()
 *   - Do ~100,000 random transformations
 *   - Then pick random test_vector
 *   - Invert the matrix
 *   - Check matrix.apply(inv.apply(test_vector)) == test_vector
 *   - Check inv.compose(matrix) == diag
 ******************************************************************************/
void test_invert_matrix(void)
{
    // 1) matrix = diag
    Matrix matrix = matrix_diag();

    // 2) ~100,000 random transformations
    //    We'll replicate the Rust approach:
    for (int iter = 0; iter < 100000; iter++) {
        // Build a pseudo-random "r"
        // If you have a better RNG in your code, use it.
        uint64_t r = ((uint64_t)rand() << 32) ^ rand();

        // r1 = r % (1<<4)
        // r2 = (r >> 4) & 127
        // r3 = (r >> 32) & 127
        size_t r1 = (size_t)(r & 0xF);
        size_t r2 = (size_t)((r >> 4) & 127U);
        size_t r3 = (size_t)((r >> 32) & 127U);

        if (r1 == 0) {
            matrix_swap_cols(&matrix, r2, r3);
        } else {
            if (r2 != r3) {
                matrix_triang(&matrix, r2, r3);
            }
        }
    }

    // 3) A random test_vector
    F128 test_vector = f128_rand();

    // 4) Invert the matrix
    Matrix inv;
    bool ok = matrix_inverse(&matrix, &inv);
    TEST_ASSERT_TRUE_MESSAGE(ok, "test_invert_matrix: matrix not invertible (unexpected)!");

    // 5) Check matrix.apply(inv.apply(test_vector)) == test_vector
    {
        F128 inner = matrix_apply(&inv, test_vector);
        F128 result = matrix_apply(&matrix, inner);

        // Compare 'result' with 'test_vector'.
        // If you have a helper like "bool bf128_eq(a,b)", you can do:
        //   TEST_ASSERT_TRUE(bf128_eq(result, test_vector));
        // Or compare their bits:
        bool rbits[128], tbits[128];
        binary128_to_bits(result, rbits);
        binary128_to_bits(test_vector, tbits);
        for (int i = 0; i < 128; i++) {
            if (rbits[i] != tbits[i]) {
                TEST_FAIL_MESSAGE("test_invert_matrix: mismatch in apply->apply");
            }
        }
    }

    // 6) Check inv.compose(&matrix) == diag
    {
        Matrix composed = matrix_compose(&inv, &matrix);
        Matrix diag = matrix_diag();

        // Check each column
        for (int i = 0; i < 128; i++) {
            uint64_t clo, chi, dlo, dhi;
            f128_get(&composed.cols[i], &clo, &chi);
            f128_get(&diag.cols[i],     &dlo, &dhi);
            if ((clo != dlo) || (chi != dhi)) {
                TEST_FAIL_MESSAGE("test_invert_matrix: inv.compose(matrix) != diag()");
            }
        }
    }

    TEST_PASS();
}
