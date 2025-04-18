#ifndef COMPRESSED_POLY_H
#define COMPRESSED_POLY_H
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>  
#include <assert.h>
#include <arm_neon.h>
#include "types.h"
#include "field.h"
#include "univariate_poly.h"

typedef struct{
    Points *compressed_coeff; // compressed polynomial
    F128 sum;          // sum of the original polynomial

} CompressedPoly;


UnivariatePolynomial* uncompress_poly(
    const CompressedPoly* compressed_poly );
CompressedPoly* compress_poly(const Points* poly);



#endif