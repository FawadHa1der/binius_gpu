#ifndef AND_PACKAGE_H
#define AND_PACKAGE_H


#include "field.h"

// #define I 2
// #define O 1

typedef struct {
    size_t input_size;
    size_t output_size;

} Algebraic_Params;


typedef struct {

    void (*quadratic)(
        const Algebraic_Params* params,
        const Points* data,
        Points* out
    );

    void (*algebraic)(
        const Algebraic_Params* params,
        const Points* data,
        size_t idx_a,
        size_t offset,
        F128 *ret[3] // Always [3][x] dims TODO a better data structure for a 2D array ?
    );

    void (*linear)(
        const Algebraic_Params* params,
        const Points* data,
        Points* out
    );

} Algebraic_Functions;

void and_package_algebraic(
    const Algebraic_Params* params,
    const Points* data,
    size_t idx_a,
    size_t offset,
    F128 *ret[3] // Always [3][x] dims TODO a better data structure for a 2D array ?
);

void and_package_linear(
    const Algebraic_Params* params,
    const Points* data,
    Points* out
);

void and_package_quadratic(
    const Algebraic_Params* params,
    const Points* data,
    Points* out
);


#endif // AND_PACKAGE_H