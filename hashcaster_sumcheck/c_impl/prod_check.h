#ifndef PROD_CHECK_H
#define PROD_CHECK_H
#include "field.h"
#include "types.h"
#include "mle_poly.h"
#include "compressed_poly.h"
#include "evaluations.h"

typedef struct {
    size_t N; // how many polynomials
    MLE_POLY* p_polys; // array of length N
    MLE_POLY* q_polys; // array of length N
    F128 claim;
    Points* challenges;
    size_t num_vars;
    int has_cached; // to replicate "Option<CompressedPoly<2>>"
    CompressedPoly* cached_round_msg; // if has_cached=1, then it's valid
} ProdCheck;

typedef struct {
    Evaluations *p_evaluations;
    Evaluations *q_evaluations;
} ProdCheckOutput;

ProdCheck prodcheck_new(MLE_POLY* p_arr, 
    MLE_POLY* q_arr,
                        size_t N,
                        F128 claim,
                        int check_init_claim);
CompressedPoly* prodcheck_round_polynomial(ProdCheck* pc);
void prodcheck_bind(ProdCheck* pc, F128 r, int challenege_index);
ProdCheckOutput prodcheck_finish(ProdCheck pc);                        

void prodcheck_free(ProdCheck pc);

#endif