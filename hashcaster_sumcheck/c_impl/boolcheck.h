#ifndef BOOLCHECK_H
#define BOOLCHECK_H
#include <stdint.h>
#include "field.h"
#include "evaluations.h"
#include "mle_poly.h"
#include "types.h"
#include "univariate_poly.h"
#include "and_package.h"

// #define C_INITIAL_ROUNDS 5 // referred to as C in the paper

typedef struct {
    int c_initial_rounds; // Number of rounds
    const Points* points;
    Points* gammas; // Folding challenge values
    Points* claims; // Initial claims
    MLE_POLY_SEQUENCE* polys; // Multilinear inputs
    const Algebraic_Params* algebraic_operations; // Input, output sizes
    const Algebraic_Functions* algebraic_functions; // Function pointer set
} BoolCheckBuilder;

BoolCheckBuilder* bool_check_builder_new(
    int C_INITIAL_ROUNDS,
    const Points* points,
    Points* claims,
    MLE_POLY_SEQUENCE* polys,
    const Algebraic_Params* algebraic_operations,
    const Algebraic_Functions* algebraic_functions

);

void bool_check_builder_free(BoolCheckBuilder* builder) ;

void compute_trit_mappings(size_t c_rounds, uint16_t** bit_mapping_out, size_t* bit_len,
                           uint16_t** trit_mapping_out, size_t* trit_len);

F128* extend_n_tables(
    const MLE_POLY_SEQUENCE* polys,  // array of N pointers to F128 arrays
    size_t N,
    size_t dims,
    size_t C,
    const uint16_t* trit_mapping,
    BoolCheckBuilder* builder,
    F128 (*linear_compressed)(
        const Points* gammas,
        const Algebraic_Params* params,
        const Algebraic_Functions* alg_funcs,
        const Points* arg
    ),
    F128 (*quadratic_compressed)(
        const Points* gammas,
        const Algebraic_Params* params,
        const Algebraic_Functions* alg_funcs,
        const Points* arg
    )
) ;                           

F128 quadratic_compressed(
    const Points* gammas,
    const Algebraic_Params* params,
    const Algebraic_Functions* alg_funcs,
    const Points* arg
);

F128 linear_compressed(
    const Points* gammas,
    const Algebraic_Params* params,
    const Algebraic_Functions* alg_funcs,
    const Points* arg
);

void algebraic_compressed(
    const Points* gammas,
    const Algebraic_Params* params,
    const Algebraic_Functions* alg_funcs,
    const Points* data,
    size_t idx_a,
    size_t offset,
    F128 result[3]
) ;
#endif