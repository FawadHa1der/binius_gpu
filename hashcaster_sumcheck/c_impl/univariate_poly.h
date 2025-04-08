#ifndef UNIVARIATE_POLY_H
#define UNIVARIATE_POLY_H
#include "field.h"
#include "types.h"


F128 polynomial_evaluate_at(const UnivariatePolynomial *poly, const F128 at);
UnivariatePolynomial* from_evaluations_deg2(const Points* evaluations);
UnivariatePolynomial* multiply_degree2_by_degree1(
    const UnivariatePolynomial *lhs,
    const UnivariatePolynomial *rhs
);
#endif