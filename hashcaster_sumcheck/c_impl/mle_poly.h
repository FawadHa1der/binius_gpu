#ifndef MLE_POLY_H
#define MLE_POLY_H
#include "field.h"

// multi linear lagrange poly
typedef struct  {
    F128 *coeffs;
    size_t len;
} MLE_POLY ;

//  points/vector of F128s
typedef struct  {
    F128 *elems;
    size_t len;
} Points ;

typedef struct  {
    MLE_POLY *mle_poly; // array of MLE_POLY
    size_t len;
} MLE_POLY_SEQUENCE ;


MLE_POLY* points_to_eq_poly(const Points* points);

F128 points_eq_eval(const Points* a, const Points* b);
Points points_default(void);
MLE_POLY_SEQUENCE* to_eq_poly_sequence(const Points *points);
bool mle_poly_eq(const MLE_POLY *a, const MLE_POLY *b);
#endif // MLE_POLY_H