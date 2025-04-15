#ifndef MLE_POLY_H
#define MLE_POLY_H
#include "field.h"
#include "types.h"
#include "stdint.h"
#include <stdbool.h>
#include <stdlib.h>
#include <memory.h>  

// multi linear lagrange poly
typedef struct  {
    F128 *coeffs;
    size_t len;
} MLE_POLY ;

typedef struct  {
    MLE_POLY *mle_poly; // array of MLE_POLY
    size_t len;
} MLE_POLY_SEQUENCE ; // TODO replace all arrays of MLE_POLY with MLE_POLY_SEQUENCE

typedef struct  {
    Points *array_of_points ; 
    size_t len;
} INVERSE_ORBIT_POINTS; // basically a vector of Points. 2D array of F128s


// A helper to get the i-th coefficient
static inline F128 mlp_get(const MLE_POLY* poly, size_t i){
    return poly->coeffs[i];
}

static inline void mlp_set(MLE_POLY* poly, size_t i, F128 val){
    poly->coeffs[i] = val;
}


// A constructor from an existing array
static inline MLE_POLY mlp_with_array(const F128* arr, size_t length)
{
    MLE_POLY poly;
    // poly = (MLE_POLY*)malloc(sizeof(MLE_POLY));
    // if (poly == NULL) {
    //     // handle allocation failure
    //     poly->len = 0;
    //     poly->coeffs = NULL;
    //     return poly;
    // }
    poly.len = length;
    if (length == 0) {
        poly.coeffs = NULL;
        return poly;
    }
    poly.coeffs = (F128*)malloc(sizeof(F128) * length);
    if (poly.coeffs == NULL) {
        // handle allocation failure
        poly.len = 0;
        return poly;
    }
    memcpy(poly.coeffs, arr, sizeof(F128)* length);
    return poly;
}

// A destructor
static inline void mlp_free(MLE_POLY* poly)
{
    if (poly == NULL) {
        return;
    }
    if (poly->coeffs != NULL) {
        free(poly->coeffs);
    }
    free(poly);
}


typedef Points Evaluations;


/// UTILS
void v_slli_epi64_c(int K, const uint8_t *x);
void drop_top_bit(uint8_t x, uint8_t *result, uint8_t *bit_index);
int cpu_v_movemask_epi8(const uint8_t* x);

///



Points* clone_points(const Points *src);
MLE_POLY* points_to_eq_poly(const Points* points);
F128 points_eq_eval(const Points* a, const Points* b);
Points points_default(void);
MLE_POLY_SEQUENCE* points_to_eq_poly_sequence(const Points *points);
bool mle_poly_eq(const MLE_POLY *a, const MLE_POLY *b);
INVERSE_ORBIT_POINTS* to_f128_inv_orbit(const Points* input_points);
Points* eq_sums(const  MLE_POLY *eq_poly);
Evaluations* restrict_polynomials(
    const MLE_POLY *polys, // array of polynomials
    size_t N,                 // how many polynomials
    const Points *challenges, // challenge points
    size_t dims               // dimension
);
F128 evaluate_at(
    const MLE_POLY *poly, 
    const Points *points
);
#endif // MLE_POLY_H