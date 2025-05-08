#include "types.h"
#include "field.h"
#include <stdlib.h>
#include <string.h>

void points_free(Points* points){
    if(points->elems != NULL){
        free(points->elems);
        points->elems = NULL;
    }
    points->len = 0;
    free(points);
    points = NULL;
}

// copy points
Points* points_copy(const Points* src) {
    if (src == NULL) {
        return NULL;
    }
    Points* dst = (Points*)malloc(sizeof(Points));
    if (dst == NULL) {
        // handle allocation failure
        return NULL;
    }
    dst->len = src->len;
    dst->elems = (F128*)malloc(sizeof(F128) * src->len);
    if (dst->elems == NULL) {
        // handle allocation failure
        free(dst);
        return NULL;
    }
    memcpy(dst->elems, src->elems, sizeof(F128) * src->len);
    return dst;

}

Points* points_init(size_t len, F128 value){
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
    for (size_t i = 0; i < len; i++) {
        pts->elems[i] = value;
    }
    return pts;
}



void points_split_at(const Points* input, size_t split_index, Points** pt_active, Points** pt_dormant) {
    assert(input != NULL);
    assert(pt_active != NULL && pt_dormant != NULL);
    assert(split_index <= input->len);

    *pt_active = (Points*)malloc(sizeof(Points));
    *pt_dormant = (Points*)malloc(sizeof(Points));
    if (*pt_active == NULL || *pt_dormant == NULL) {
        // handle allocation failure
        return;
    }

    // Allocate memory for pt_active and copy first `split_index` elements
    (*pt_active)->len = split_index;
    (*pt_active)->elems = malloc(split_index * sizeof(F128));
    for (size_t i = 0; i < split_index; i++) {
        (*pt_active)->elems[i] = input->elems[i];
    }

    // Allocate memory for pt_dormant and copy the rest
    size_t rest_len = input->len - split_index;
    (*pt_dormant)->len = rest_len;
    (*pt_dormant)->elems = malloc(rest_len * sizeof(F128));
    for (size_t i = 0; i < rest_len; i++) {
        (*pt_dormant)->elems[i] = input->elems[split_index + i];
    }
}


bool points_equal(const Points* a, const Points* b) {
    if (a == NULL || b == NULL) return false;
    if (a->len != b->len) return false;

    for (size_t i = 0; i < a->len; ++i) {
        if (!f128_eq(a->elems[i], b->elems[i])) {
            return false;
        }
    }
    return true;
}