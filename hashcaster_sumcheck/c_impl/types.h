#ifndef TYPES_H
#define TYPES_H

#include <stdint.h>
#include <stdbool.h>
#include <arm_neon.h>
#include <stdio.h>
#include <stdlib.h>


// For a GF(2^128) element, we store two 64-bit words:
//   hi = the high 64 bits
//   lo = the low 64 bits
typedef struct {
    uint64_t low;
    uint64_t high;
} uint128_t;

typedef uint128_t F128;

// typedef struct {
//     void* elems;
//     size_t len;
// } Vector; // generic type for all arrays. 

/*---------------------------
 * Helper union to reinterpret
 * between (uint128_t) and (uint8x16_t).
 *---------------------------*/
typedef union {
    uint8x16_t vec;   // 128-bit NEON vector
    uint64x2_t u64x2; // Also 128 bits, but two 64-bit lanes
    struct {
        uint64_t lo;
        uint64_t hi;
    };
} u128_as_vec;



typedef struct {
    F128* elems;
    size_t len;
} Points; 


typedef Points UnivariatePolynomial;


static inline void points_free(Points* points){
    if(points->elems != NULL){
        free(points->elems);
        points->elems = NULL;
    }
    points->len = 0;
    free(points);
    points = NULL;
}

static inline Points* points_init(size_t len){
    Points* pts = (Points*)malloc(sizeof(Points));
    if (pts == NULL) {
        // handle allocation failure
        pts->len = 0;
        pts->elems = NULL;
        return pts;
    }
    pts->len = len;
    if (len == 0) {
        pts->elems = NULL;
        return pts;
    }
    pts->elems = (F128*)malloc(sizeof(F128) * len);
    if (pts->elems == NULL) {
        // handle allocation failure
        pts->len = 0;
        return pts;
    }
    return pts;
}



#endif