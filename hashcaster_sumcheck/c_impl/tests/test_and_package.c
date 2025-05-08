#include "test_and_package.h"
#include "unity.h"

void test_exec_alg_and(void) {
    // F128 a1 = f128_from_raw(0b1010101010101010101010101010101010101010101010101010101010101010, 0b1100110011001100110011001100110011001100110011001100110011001100);
    // F128 a2 = f128_from_raw(0b1100110011001100110011001100110011001100110011001100110011001100, 0b1010101010101010101010101010101010101010101010101010101010101010);


    F128 a1 = f128_rand();
    F128 a2 = f128_rand();

    // Step 2: Define test vectors
    //F128 a[2] = { a1, a2 };
    Points* a = points_init(2, f128_zero());
    a->elems[0] = a1;
    a->elems[1] = a2;


    // F128 b[2] = {
    //     f128_shr(a1, 1),
    //     f128_shr(a2, 1),
    // };
    Points* b = points_init(2, f128_zero());
    b->elems[0] = f128_shr(a1, 1);
    b->elems[1] = f128_shr(a2, 1);
    // F128 c[2] = {
    //     a1,
    //     f128_shr(a2, 1),
    // };
    Points* c = points_init(2, f128_zero());
    c->elems[0] = a1;
    c->elems[1] = f128_shr(a2, 1);
    // F128 d[2] = {
    //     f128_shr(a1, 1),
    //     a2,
    // };
    Points* d = points_init(2, f128_zero());
    d->elems[0] = f128_shr(a1, 1);
    d->elems[1] = a2;
    
    Algebraic_Params params;
    params.input_size = 2;
    params.output_size = 1;
    // Step 3: Flatten a into bitsliced coordinates
    Points* input_coords = points_init(2* 128, f128_zero());

    for (int j = 0; j < 2; j++) {
        uint128_t val = (a->elems[j]);
        for (int i = 0; i < 128; i++) {
            input_coords->elems[j * 128 + i] = f128_bitand(f128_shr(val , i), f128_from_uint64(1)).low ? f128_one() : f128_zero();
        }
    }

    // Step 4: Call algebraic implementation
    // F128 rhs[3][1];
    F128* rhs[3];
    for (int i = 0; i < 3; i++) {
        rhs[i] = (F128*)malloc(sizeof(F128) * 1);
        rhs[i][0] = f128_zero();
    }
    and_package_algebraic(&params, input_coords, 0, 1, rhs);

    // Step 5: Compute expected quadratic results
    // F128 a_quad[1], b_quad[1], c_quad[1], d_quad[1];
    Points* a_quad = points_init(1, f128_zero());
    Points* b_quad = points_init(1, f128_zero());
    Points* c_quad = points_init(1, f128_zero());
    Points* d_quad = points_init(1, f128_zero());
    
    and_package_quadratic(&params, a, a_quad);
    and_package_quadratic(&params, b, b_quad);
    and_package_quadratic(&params, c, c_quad);
    and_package_quadratic(&params, d, d_quad);

    // Step 6: Check algebraic result matches individual quadratics
    TEST_ASSERT_TRUE(f128_eq(rhs[0][0], a_quad->elems[0]));
    TEST_ASSERT_TRUE(f128_eq(rhs[1][0], b_quad->elems[0]));

    F128 sum = f128_add(f128_add(a_quad->elems[0], b_quad->elems[0]),
                        f128_add(c_quad->elems[0], d_quad->elems[0]));
    TEST_ASSERT_TRUE(f128_eq(rhs[2][0], sum));

    //free
    points_free(a);
    points_free(b);
    points_free(c);
    points_free(d);
    points_free(a_quad);
    points_free(b_quad);
    points_free(c_quad);    
    points_free(d_quad);
    points_free(input_coords);
    free(rhs[0]);
    free(rhs[1]);
    free(rhs[2]);
    // free(params);
}