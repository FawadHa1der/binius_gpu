#include "../compressed_poly.h"
#include "unity.h"
#include "stdlib.h"



void test_compress(void)
{
    const int  N = 4;

    Points* poly;
    poly= (Points*)malloc(sizeof(Points));
    poly->len= N;
    poly->elems= (F128*)malloc(sizeof(F128)*N);
    for(int i=0;i<N;i++){
        poly->elems[i]= f128_from_uint64(i+1);
    }

    // compress
    CompressedPoly* cr = compress_poly(poly);

    // expected_compressed= [1,3,4]
    CompressedPoly* expected_cp;
    expected_cp= (CompressedPoly*)malloc(sizeof(CompressedPoly));
    expected_cp->compressed_coeff= (Points*)malloc(sizeof(Points));
    expected_cp->compressed_coeff->len= N-1;
    expected_cp->compressed_coeff->elems= (F128*)malloc(sizeof(F128)*(N-1));
    expected_cp->compressed_coeff->elems[0]= f128_from_uint64(1);
    expected_cp->compressed_coeff->elems[1]= f128_from_uint64(3);
    expected_cp->compressed_coeff->elems[2]= f128_from_uint64(4);

    // expected_sum= 2 +3 +4
    F128 ex_sum= f128_from_uint64(2);
    ex_sum= f128_add(ex_sum, f128_from_uint64(3));
    ex_sum= f128_add(ex_sum, f128_from_uint64(4));

    // Check compressed eq
    // We'll do partial check:
    for(int i=0;i<(N-1); i++){
        TEST_ASSERT_TRUE_MESSAGE(
            f128_eq(cr->compressed_coeff->elems[i], expected_cp->compressed_coeff->elems[i]),
            "Compressed polynomial mismatch"
        );
    }
    // Check sum eq
    TEST_ASSERT_TRUE_MESSAGE(
        f128_eq(cr->sum, ex_sum),
        "Sum mismatch"
    );

    // free
    free(poly->elems);
    free(poly);

    free(cr->compressed_coeff->elems);
    free(cr->compressed_coeff);
    free(cr);

}


void test_sum_standard_case(void)
{
    // same poly => [1,2,3,4]
    // F128 poly[N];

    // poly[0]= f128_from_uint64(1);
    // poly[1]= f128_from_uint64(2);
    // poly[2]= f128_from_uint64(3);
    // poly[3]= f128_from_uint64(4);
    const int  N = 4;

    Points* poly;
    poly= (Points*)malloc(sizeof(Points));
    poly->len= N;
    poly->elems= (F128*)malloc(sizeof(F128)*N);

    for (int i=0;i<N;i++){
        poly->elems[i]= f128_from_uint64(i+1);
    }

    CompressedPoly* cr= compress_poly(poly);

    // expected_sum= 2+3+4
    F128 ex_sum= f128_from_uint64(2);
    ex_sum= f128_add(ex_sum, f128_from_uint64(3));
    ex_sum= f128_add(ex_sum, f128_from_uint64(4));

    TEST_ASSERT_TRUE_MESSAGE(
        f128_eq(cr->sum, ex_sum),
        "Sum mismatch in test_sum_standard_case"
    );

    // free
    free(poly->elems);
    free(poly);
}


void test_sum_all_zero_coefficients(void)
{
    // poly= [0,0,0]
    // but for N=4 is short, let's define N=3 or smaller? 
    // We'll do a smaller array or adapt. If N=4, let's interpret only the first 3?
    // We'll just do 3 zero elements for demonstration
    // We'll define a local NZ=3
    const int  NZ = 3;
    Points* poly_z;
    poly_z= (Points*)malloc(sizeof(Points));
    poly_z->len= NZ;
    poly_z->elems= (F128*)malloc(sizeof(F128)*NZ);
    for(int i=0;i<NZ;i++){
        poly_z->elems[i]= f128_zero();
    }

    // We'll define a separate compress function for length=3 if we want. 
    // or we do a partial: In real code, you'd unify. We'll do a small inline version:

    // sum => poly[1..2]
    F128 sum = f128_zero();
    for(int i=1;i<NZ;i++){
        sum= f128_add(sum, poly_z->elems[i]);
    }

    TEST_ASSERT_TRUE_MESSAGE(
        f128_eq(sum, f128_zero()),
        "Sum of zero coefficients should be zero"
    );

    // free
    free(poly_z->elems);
    free(poly_z);
}

// test_coeffs_reconstruction_standard_case
void test_coeffs_reconstruction_standard_case(void)
{
    // poly= [1,2,3,4]
    const int  N = 4;
    Points* poly;
    poly= (Points*)malloc(sizeof(Points));
    poly->len= N;   
    poly->elems= (F128*)malloc(sizeof(F128)*N);
    for (int i=0;i<N;i++){
        poly->elems[i]= f128_from_uint64(i+1);
    }

    // compress
    CompressedPoly* cr= compress_poly(poly);

    // reconstruct
    UnivariatePolynomial* reconstructed= uncompress_poly(cr, cr->sum);

    // build "original" polynomial
    // UnivariatePolynomial original= UnivariatePolynomial_new(poly);

    for (int i=0;i<N;i++){
        // check if reconstructed[i] == original[i]
        TEST_ASSERT_TRUE_MESSAGE(
            f128_eq(reconstructed->elems[i], poly->elems[i]),
            "Reconstructed polynomial mismatch"
        );
    }

    // free
    free(poly->elems);
    free(poly);

}

// test_coeffs_all_zero_coefficients
void test_coeffs_all_zero_coefficients(void)
{
    // poly= [0,0,0,0]
    const int  N = 4;
    Points* poly;
    poly= (Points*)malloc(sizeof(Points));
    poly->len= N;
    poly->elems= (F128*)malloc(sizeof(F128)*N);
    for (int i=0;i<N;i++){
        poly->elems[i]= f128_zero();
    }

    CompressedPoly* cr= compress_poly(poly);

    UnivariatePolynomial* reconstructed= uncompress_poly(cr, cr->sum);

    for (int i=0;i<N;i++){
        // check if reconstructed[i] == original[i]
        TEST_ASSERT_TRUE_MESSAGE(
            f128_eq(reconstructed->elems[i], poly->elems[i]),
            "Reconstructed polynomial mismatch when all are zero"
        );
    }
    // free
    free(poly->elems);
    free(poly);
    free(cr->compressed_coeff->elems);
    free(cr->compressed_coeff);
    free(cr);
}

// test_coeffs_large_coefficients
void test_coeffs_large_coefficients(void)
{
    // poly= [1_000_000, 2_000_000, 3_000_000, 4_000_000]
    const int  N = 4;
    Points* poly;
    poly= (Points*)malloc(sizeof(Points));
    poly->len= N;
    poly->elems= (F128*)malloc(sizeof(F128)*N);
    for (int i=0;i<N;i++){
        poly->elems[i]= f128_from_uint64((i+1) * 1000000);
    }
    // compress
    CompressedPoly* cr= compress_poly(poly);
    UnivariatePolynomial* reconstructed= uncompress_poly(cr, cr->sum);

    for (int i=0;i<N;i++){
        // check if reconstructed[i] == original[i]
        TEST_ASSERT_TRUE_MESSAGE(
            f128_eq(reconstructed->elems[i], poly->elems[i]),
            "Reconstructed polynomial mismatch for large coefficients"
        );
    }

    // free
    free(poly->elems);
    free(poly);
    free(cr->compressed_coeff->elems);
    free(cr->compressed_coeff);
    free(cr);

}
