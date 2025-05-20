#ifndef BOOLCHECK_H
#define BOOLCHECK_H
#include <stdint.h>
#include "field.h"
#include "evaluations.h"
#include "mle_poly.h"
#include "types.h"
#include "univariate_poly.h"
#include "and_package.h"
#include "compressed_poly.h"

// #define C_INITIAL_ROUNDS 5 // referred to as C in the paper

typedef struct {
    size_t c_initial_rounds; // Number of rounds
    const Points* points;
    Points* gammas; // Folding challenge values
    Points* claims; // Initial claims
    MLE_POLY_SEQUENCE* polys; // Multilinear inputs
    const Algebraic_Params* algebraic_operations; // Input, output sizes
    const Algebraic_Functions* algebraic_functions; // Function pointer set
} BoolCheckBuilder;

BoolCheckBuilder* bool_check_builder_new(
    size_t C_INITIAL_ROUNDS,
    const Points* points,
    Points* claims,
    MLE_POLY_SEQUENCE* polys,
    const Algebraic_Params* algebraic_operations,
    const Algebraic_Functions* algebraic_functions

);

// pub struct BoolCheck<'a, const N: usize, const M: usize, const C: usize, A: AlgebraicOps<N, M>> {
//     /// The evaluation points for the multilinear polynomials.
//     pub points: &'a Points,

//     /// An array of multilinear polynomials used in the protocol.
//     pub polys: &'a [MultilinearLagrangianPolynomial; N],

//     /// A vector storing intermediate evaluations of the polynomials.
//     pub extended_table: Vec<BinaryField128b>,

//     /// An optional field storing evaluations on a restricted subset of the hypercube.
//     poly_coords: Evaluations,

//     /// A sequence of random challenges provided by the verifier.
//     pub challenges: Points,

//     /// A vector of bit mappings for optimized indexing of polynomial coefficients.
//     pub bit_mapping: Vec<u16>,

//     /// A sequence of equality polynomials used to verify claims.
//     pub eq_sequence: Vec<MultilinearLagrangianPolynomial>,

//     /// A vector of compressed polynomials computed at each round.
//     pub round_polys: Vec<CompressedPoly<3>>,

//     /// The current claim being verified in the protocol.
//     pub claim: BinaryField128b,

//     /// An array of folding challenges used for polynomial compression.
//     pub gammas: [BinaryField128b; M],

//     /// Abstract algebraic operations.
//     pub algebraic_operations: &'a A,
// }

typedef struct {
    const Points* points;                        // length = num_vars
    const MLE_POLY_SEQUENCE* polys;              // array of N polynomials
    Points* extended_table;                        // flattened extended table
    Evaluations* poly_coords;                     // restricted evaluations
    Points* challenges;                           // folding points from verifier
    uint16_t* bit_mapping;                       // used for table lookups
    size_t bit_mapping_len;
    MLE_POLY_SEQUENCE *eq_sequence;               // equality poly sequence
    CompressedPoly* round_polys;          // array of CompressedPoly<2>
    size_t round_polys_len;
    F128 claim;
    Points* gammas;                              // folding challenges
    const Algebraic_Params* algebraic_operations;
    const Algebraic_Functions* algebraic_functions;
    size_t num_vars;
    size_t c_param;
} BoolCheck;


// pub struct BoolCheckOutput<const N: usize, const M: usize> {
//     /// Evaluations of the polynomials on a Frobenius subdomain.
//     pub frob_evals: FixedEvaluations<N>,

//     /// A vector of compressed polynomials computed during the protocol's rounds.
//     pub round_polys: Vec<CompressedPoly<M>>,
// }

typedef struct {
    Evaluations* frob_evals; // Frobenius evaluations
    CompressedPoly* round_polys;  // round polynomials
    size_t round_polys_len;
} BoolCheckOutput;


void bool_check_builder_free(BoolCheckBuilder* builder) ;

void compute_trit_mappings(size_t c_rounds, uint16_t** bit_mapping_out, size_t* bit_len,
                           uint16_t** trit_mapping_out, size_t* trit_len);

Points* extend_n_tables(
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

BoolCheck* boolcheck_new(
    BoolCheckBuilder* builder,
    const F128 gamma 
);

CompressedPoly* boolcheck_round_polynomial(BoolCheck* bc) ;
void boolcheck_bind(BoolCheck* bc, const F128* r);
BoolCheckOutput* boolcheck_finish(BoolCheck* bc);
void bool_check_builder_free(BoolCheckBuilder* builder);
#endif