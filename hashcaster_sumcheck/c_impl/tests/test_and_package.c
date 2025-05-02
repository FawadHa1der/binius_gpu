#include "test_and_package.h"
#include "unity.h"

void test_exec_alg_and(void) {
    // Step 1: Generate two random field elements
    F128 a1 = f128_rand();
    F128 a2 = f128_rand();

    // Step 2: Define test vectors
    F128 a[2] = { a1, a2 };
    F128 b[2] = {
        f128_shr(a1, 1),
        f128_shr(a2, 1),
    };
    F128 c[2] = {
        a1,
        f128_shr(a2, 1),
    };
    F128 d[2] = {
        f128_shr(a1, 1),
        a2,
    };

    // Step 3: Flatten a into bitsliced coordinates
    F128 input_coords[2 * 128];  // 2 elements × 128 bits
    for (int j = 0; j < 2; j++) {
        uint128_t val = (a[j]);
        for (int i = 0; i < 128; i++) {
            input_coords[j * 128 + i] = f128_bitand(f128_shr(val , i), f128_from_uint64(1));
        }
    }

    // Step 4: Call algebraic implementation
    F128 rhs[3][1];
    and_package_algebraic(input_coords, 0, 1, rhs);

    // Step 5: Compute expected quadratic results
    F128 a_quad[1], b_quad[1], c_quad[1], d_quad[1];
    and_package_quadratic(a, a_quad);
    and_package_quadratic(b, b_quad);
    and_package_quadratic(c, c_quad);
    and_package_quadratic(d, d_quad);

    // Step 6: Check algebraic result matches individual quadratics
    TEST_ASSERT_TRUE(f128_eq(rhs[0][0], a_quad[0]));
    TEST_ASSERT_TRUE(f128_eq(rhs[1][0], b_quad[0]));

    F128 sum = f128_add(f128_add(a_quad[0], b_quad[0]),
                        f128_add(c_quad[0], d_quad[0]));
    TEST_ASSERT_TRUE(f128_eq(rhs[2][0], sum));
}