#include "compressed_poly.h"
#include "field.h"
#include "stdlib.h"
#include "string.h"
#include "stdio.h"

CompressedPoly* compress_poly(const Points* poly)
{
    // 1) sum = poly[1..N-1] => skip(1), fold
    F128 sum = f128_zero();
    for(size_t i=1; i<poly->len; i++){
        sum = f128_add(sum, poly->elems[i]);
    }

    // 2) build compressed array => size (N-1)
    //    compressed[0] = poly[0]
    //    compressed[1..(N-2)] copy from poly[2..(N-1)]
    CompressedPoly* compressed;
    compressed = (CompressedPoly*)malloc(sizeof(CompressedPoly));
    if (!compressed) {
        // handle error
        return NULL;
    }
    compressed->compressed_coeff = (Points*)malloc(sizeof(Points));
    if (!compressed->compressed_coeff) {
        // handle error
        free(compressed);
        return NULL;
    }

    compressed->compressed_coeff->len = poly->len - 1;
    compressed->compressed_coeff->elems = (F128*)malloc(sizeof(F128) * compressed->compressed_coeff->len);
    if (!compressed->compressed_coeff->elems) {
        // handle error
        free(compressed);
        return NULL;
    }
    compressed->compressed_coeff->elems[0] = poly->elems[0];

    //->elems poly[2..N-1] => compressed.data[1..(N-2)]
    // number of elements to copy => (N-1) - 1 => N-2
    // i in [2.. (N-1)]
    for(size_t i=2; i<poly->len; i++){
        compressed->compressed_coeff->elems[i-1] = poly->elems[i];
    }

    // 3) return the tuple => (compressed, sum)
    compressed->sum = sum;
    return compressed;
}

// returns a "UnivariatePolynomial".
UnivariatePolynomial* uncompress_poly(
    const CompressedPoly* compressed_poly )
{
    // step 1) c0= compressed_poly->compressed_coeff->elems[0]
    F128 c0= compressed_poly->compressed_coeff->elems[0];

    // step 2) ev_1= c0 + sum
    F128 ev_1= f128_add(c0, compressed_poly->sum);

    // step 3) c1= (sum of all stored in compressed_poly->compressed_coeff->elems) + ev_1
    //   self.iter().fold(ZERO, |a,b| a+b) => sum all compressed_poly->compressed_coeff->elems
    F128 partial_sum = f128_zero();
    for(size_t i=0; i<compressed_poly->compressed_coeff->len; i++){
        partial_sum = f128_add(partial_sum, compressed_poly->compressed_coeff->elems[i]);
    }
    // c1= partial_sum + ev_1
    F128 c1 = f128_add(partial_sum, ev_1);

    // step 4) build full => size N+1
    UnivariatePolynomial* full;
    full = (UnivariatePolynomial*)malloc(sizeof(UnivariatePolynomial));
    if (!full) {
        // handle error
        return NULL;
    }
    full->len = compressed_poly->compressed_coeff->len + 1;
    full->elems = (F128*)malloc(sizeof(F128) * full->len);
    if (!full->elems) {
        // handle error
        free(full);
        return NULL;
    }
    // full[0]= c0
    full->elems[0] = c0;
    // full[1]= c1
    full->elems[1] = c1;
    // full[2..] = compressed_poly->compressed_coeff->elems[1..(N-1)]
    // we copy from compressed_poly->compressed_coeff->elems[1..(N-1)] => full->elems[2..(N)]
    // that is (N-2) elements
    for(size_t i=1; i<(compressed_poly->compressed_coeff->len ); i++){
        full->elems[i+1] = compressed_poly->compressed_coeff->elems[i];
    }

    // return
    return full;
}


bool compressed_poly_eq(const CompressedPoly *compressed_poly_a, const CompressedPoly *compressed_poly_b)
{
    if (compressed_poly_a->compressed_coeff->len != compressed_poly_b->compressed_coeff->len) {
        return false;
    }

    for (size_t i = 0; i < compressed_poly_a->compressed_coeff->len; i++) {
        if (!f128_eq(compressed_poly_a->compressed_coeff->elems[i], compressed_poly_b->compressed_coeff->elems[i])) {
            return false;
        }
    }
    if (!f128_eq(compressed_poly_a->sum, compressed_poly_b->sum)) {
        return false;
    }

    return true;
}
