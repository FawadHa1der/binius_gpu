#ifndef MULTI_CLAIM_BUILDER_H
#define MULTI_CLAIM_BUILDER_H

#include "field.h"
#include "types.h"
#include "mle_poly.h"
#include "prod_check.h"


typedef struct  {
    MLE_POLY_SEQUENCE* polys; // array of MLE_POLY
    Points* points;
    Evaluations* openings; // array of size num_polys * 128
} MulticlaimBuilder;


typedef struct {
    MLE_POLY_SEQUENCE *polys;              // array of N polynomials (external, not owned)
    F128 gamma;                   // gamma used for Frobenius
    ProdCheck *object;            // the 1-polynomial prodcheck
    // size_t N;
} MultiClaim;


MulticlaimBuilder* multiclaim_builder_new(
    MLE_POLY_SEQUENCE* polys,
    Points* points,
    Evaluations* openings
)  ;

MultiClaim* multiclaim_builder_build(
    MulticlaimBuilder* builder,
    F128 gamma
);


MultiClaim* multi_claim_new(
    MLE_POLY* poly,
    const Points *points,
    const F128 *openings,
    const Points *gamma_pows,
    MLE_POLY_SEQUENCE *polys
);

void multi_claim_bind(MultiClaim *mc, const F128 challenge);


Evaluations* multi_claim_finish(MultiClaim *mc);


void multi_claim_free(MultiClaim *mc) ;
#endif