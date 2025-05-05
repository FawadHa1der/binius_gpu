#ifndef AND_PACKAGE_H
#define AND_PACKAGE_H


#include "field.h"

#define I 2
#define O 1

void and_package_algebraic(
    const F128* data,
    size_t data_len,
    size_t idx_a,
    size_t offset,
    F128 ret[3][O]
);

void and_package_linear(
    const F128 data[I],
    F128 out[O]
);

void and_package_quadratic(
    const F128 data[I],
    F128 out[O]
);


#endif // AND_PACKAGE_H