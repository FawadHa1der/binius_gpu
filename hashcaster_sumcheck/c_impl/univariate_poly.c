#include "univariate_poly.h"
#include <stdlib.h>
#include "types.h"
#include "field.h"

// Evaluate the polynomial using Horner's method:
// P(x) = coeffs[len-1] + x*( coeffs[len-2] + x*( ... + x*coeffs[0] ) )
F128 polynomial_evaluate_at(const UnivariatePolynomial *poly, const F128 at)
{
    F128 eval = f128_zero(); // Start with zero

    // We iterate from the last coefficient down to the first
    // This matches Rust's .rfold(..) over poly->coeffs
    // for i in [len-1..0]
    for(int i=(int)poly->len -1; i>=0; i--){
        // eval = (eval * at) + poly->coeffs[i];
        F128 tmp = f128_mul(eval, at);
        eval = f128_add(tmp, poly->elems[i]);
    }

    return eval;
}


// Construct a degree-2 polynomial from W(0)=c0, W(1)= c0 + c1 + c2, W(∞)= c2
// evaluations[0] => c0
// evaluations[1] => c0 + c1 + c2
// evaluations[2] => c2
UnivariatePolynomial* from_evaluations_deg2(const Points* evaluations)
{
    assert (evaluations->len == 3);
    // c0 = evaluations[0]
    F128 c0 = evaluations->elems[0];

    // c1 = (evaluations->elems[1] - c0) - evaluations->elems[2]
    // i.e. c1 = [ c0 + c1 + c2 ] - c0 - c2
    F128 tmp = f128_add(evaluations->elems[1], c0);
    F128 c1  = f128_add(tmp, evaluations->elems[2]);

    // c2 = evaluations->elems[2]
    F128 c2 = evaluations->elems[2];

    UnivariatePolynomial* poly;
    poly = (UnivariatePolynomial*)malloc(sizeof(UnivariatePolynomial));
    if (!poly) {
        // handle error
        return NULL;
    }
    poly->len = 3; // degree 2 polynomial
    poly->elems = (F128*)malloc(sizeof(F128) * poly->len);
    if (!poly->elems) {
        // handle error
        free(poly);
        return NULL;
    }
    poly->elems[0] = c0;
    poly->elems[1] = c1;
    poly->elems[2] = c2;
    return poly;
}


// Multiply a degree-2 polynomial by a degree-1 polynomial => degree-3 polynomial
// lhs: c0 + c1*x + c2*x^2
// rhs: d0 + d1*x
// result => e0 + e1*x + e2*x^2 + e3*x^3
UnivariatePolynomial* multiply_degree2_by_degree1(
    const UnivariatePolynomial *lhs,
    const UnivariatePolynomial *rhs
)
{
    // let c0= lhs->coeffs[0], c1= lhs->coeffs[1], c2= lhs->coeffs[2]
    // let d0= rhs->coeffs[0], d1= rhs->coeffs[1]
    F128 c0 = lhs->elems[0];
    F128 c1 = lhs->elems[1];
    F128 c2 = lhs->elems[2];

    F128 d0 = rhs->elems[0];
    F128 d1 = rhs->elems[1];

    // e0= c0*d0
    F128 e0 = f128_mul(c0, d0);

    // e1= c0*d1 + c1*d0
    F128 tmp1= f128_mul(c0, d1);
    F128 tmp2= f128_mul(c1, d0);
    F128 e1  = f128_add(tmp1, tmp2);

    // e2= c1*d1 + c2*d0
    tmp1= f128_mul(c1, d1);
    tmp2= f128_mul(c2, d0);
    F128 e2= f128_add(tmp1, tmp2);

    // e3= c2*d1
    F128 e3= f128_mul(c2, d1);

    UnivariatePolynomial* result;
    result = (UnivariatePolynomial*)malloc(sizeof(UnivariatePolynomial));
    if (!result) {
        // handle error
        return NULL;
    }
    result->len = 4; // degree 3 polynomial
    result->elems = (F128*)malloc(sizeof(F128) * result->len);
    if (!result->elems) {
        // handle error
        free(result);
        return NULL;
    }

    result->elems[0] = e0;
    result->elems[1] = e1;
    result->elems[2] = e2;
    result->elems[3] = e3;
    return result;
}