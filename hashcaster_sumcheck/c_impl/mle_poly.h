#ifndef MLE_POLY_H
#define MLE_POLY_H
#include "field.h"
#include "types.h"

// multi linear lagrange poly
typedef struct  {
    F128 *coeffs;
    size_t len;
} MLE_POLY ;

// typedef Vector Points;

typedef struct {
    F128* elems;
    size_t len;
} Points; 

typedef struct  {
    MLE_POLY *mle_poly; // array of MLE_POLY
    size_t len;
} MLE_POLY_SEQUENCE ;

typedef struct  {
    Points *array_of_points ; 
    size_t len;
} INVERSE_ORBIT_POINTS; // basically a vector of Points. 2D array of F128s




/// UTILS
void v_slli_epi64_c(int K, const uint8_t *x);
void drop_top_bit(uint8_t x, uint8_t *result, uint8_t *bit_index);
int cpu_v_movemask_epi8(const uint8_t* x);

///



Points* clone_points(const Points *src);
MLE_POLY* points_to_eq_poly(const Points* points);
F128 points_eq_eval(const Points* a, const Points* b);
Points points_default(void);
MLE_POLY_SEQUENCE* to_eq_poly_sequence(const Points *points);
bool mle_poly_eq(const MLE_POLY *a, const MLE_POLY *b);
INVERSE_ORBIT_POINTS* to_f128_inv_orbit(const Points* input_points);
Points* eq_sums(const  MLE_POLY *eq_poly);
Points* restrict_polynomials(
    const MLE_POLY *polys, // array of polynomials
    size_t N,                 // how many polynomials
    const Points *challenges, // challenge points
    size_t dims               // dimension
);
#endif // MLE_POLY_H