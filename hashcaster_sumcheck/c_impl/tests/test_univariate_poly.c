#include "test_univariate_poly.h"
#include "stdio.h"
#include "stdlib.h"
#include "unity.h"

void test_univariate_polynomial_evaluate_at(void)
{
    // "Define the polynomial P(x)= 3 + 2x + x^2"
    Points* coeffs;
    coeffs= (Points*)malloc(sizeof(Points));
    coeffs->len= 3;
    coeffs->elems= (F128*)malloc(sizeof(F128)*3);
    coeffs->elems[0]= f128_from_uint64(3); // a0
    coeffs->elems[1]= f128_from_uint64(2); // a1
    coeffs->elems[2]= f128_from_uint64(1); // a2

    // We'll define "evaluate" using a function or direct Horner's method inline:
    // define x=2
    F128 x= f128_from_uint64(2);

    // expected= 3 +2*2 +1*(2^2)
    // for demonstration:
    // let two= f128_from_uint64(2)
    F128 term1= f128_mul(f128_from_uint64(1), f128_mul(f128_from_uint64(2), f128_from_uint64(2))); // 2*2=4
    F128 partA= f128_add(f128_from_uint64(3), f128_mul(f128_from_uint64(2), f128_from_uint64(2))); 
    // i.e. 3 + 2*2= 3 +4=7
    F128 expected= f128_add(partA, term1); // 7 +4= 11 in a real sense

    // result = polynomial_evaluate(...) if we have a function:
    F128 result = polynomial_evaluate_at(coeffs, x);

    TEST_ASSERT_TRUE_MESSAGE(
        f128_eq(result, expected),
        "Wrong evaluation result"
    );
}

void test_from_evaluations_deg2(void)
{
    // evals= [3,9,5]
    // W(0)=3 => c0=3
    // W(1)=9 => c0 + c1 + c2=9
    // W(inf)=5 => c2=5
    // => c1= 9-3-5=1

    Points* evals;
    evals= (Points*)malloc(sizeof(Points));
    evals->len= 3;
    evals->elems= (F128*)malloc(sizeof(F128)*3);
    evals->elems[0]= f128_from_uint64(3);
    evals->elems[1]= f128_from_uint64(9);
    evals->elems[2]= f128_from_uint64(5);

    UnivariatePolynomial* poly= from_evaluations_deg2(evals);

    // expected= [3, 1, 5]
    //  c0=3, c1=9-3-5=1, c2=5
    F128 c1part = f128_add(evals->elems[1], evals->elems[0]); // 9-3=6
    F128 c1= f128_add(c1part, evals->elems[2]);       // 6-5=1

    F128 expected[3];
    expected[0]= evals->elems[0];
    expected[1]= c1;
    expected[2]= evals->elems[2];

    // check coefficients
    for(int i=0;i<3;i++){
        TEST_ASSERT_TRUE_MESSAGE(
            f128_eq(poly->elems[i], expected[i]),
            "Incorrect polynomial coefficients"
        );
    }

    // Also check evaluate_at(0)== W(0), evaluate_at(1)== W(1)
    // We'll define x=0 => P(0)= c0 => 3
    F128 zero= f128_zero();
    F128 r0= polynomial_evaluate_at(poly, zero);
    TEST_ASSERT_TRUE_MESSAGE(
        f128_eq(r0, evals->elems[0]),
        "Evaluation at t=0 failed"
    );

    // x=1 => c0+ c1+ c2 => check eq evals->elems[1]
    F128 one= f128_one();
    F128 r1= polynomial_evaluate_at(poly, one);
    TEST_ASSERT_TRUE_MESSAGE(
        f128_eq(r1,evals->elems[1]),
        "Evaluation at t=1 failed"
    );
}

void test_multiply_degree2_by_degree1(void)
{
    // define poly_deg2 => c0=3, c1=2, c2=1
    F128 c0= f128_from_uint64(3);
    F128 c1= f128_from_uint64(2);
    F128 c2= f128_from_uint64(1);

    Points* deg2 = (Points*)malloc(sizeof(Points));
    deg2->len= 3;
    deg2->elems= (F128*)malloc(sizeof(F128)*3);
    deg2->elems[0]= c0;
    deg2->elems[1]= c1;
    deg2->elems[2]= c2;


    // define poly_deg1 => d0=4, d1=5
    F128 d0= f128_from_uint64(4);
    F128 d1= f128_from_uint64(5);

    Points* deg1 = (Points*)malloc(sizeof(Points));
    deg1->len= 2;
    deg1->elems= (F128*)malloc(sizeof(F128)*2);
    deg1->elems[0]= d0;
    deg1->elems[1]= d1;

    // do multiply => result => 4 coeffs 
    // e0= c0*d0
    F128 e0= f128_mul(c0, d0);

    // e1= c0*d1 + c1*d0
    F128 tmp1= f128_mul(c0, d1);
    F128 tmp2= f128_mul(c1, d0);
    F128 e1= f128_add(tmp1, tmp2);

    // e2= c1*d1 + c2*d0
    tmp1= f128_mul(c1, d1);
    tmp2= f128_mul(c2, d0);
    F128 e2= f128_add(tmp1,tmp2);

    // e3= c2*d1
    F128 e3= f128_mul(c2,d1);

    // // We'll define the "result" array of size 4
    // F128 result[4];
    // result[0]= e0;
    // result[1]= e1;
    // result[2]= e2;
    // result[3]= e3;
    // or we can use a struct
    UnivariatePolynomial* result;
    result = multiply_degree2_by_degree1(deg2, deg1);
    


    // expected => same
    F128 expected[4];
    expected[0]= e0;
    expected[1]= e1;
    expected[2]= e2;
    expected[3]= e3;

    // compare
    for(int i=0; i<4; i++){
        TEST_ASSERT_TRUE_MESSAGE(
            f128_eq(result->elems[i], expected[i]),
            "Multiplication of degree-2 and degree-1 polynomials failed"
        );
    }
}
