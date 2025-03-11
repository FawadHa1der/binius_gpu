#include "frob_utils.h"
#include "field.h"

/* eq_poly:
 *   - pt: array of F128, length = l
 *   - l: length of pt
 *   - out_len: (output) the size of the returned array = 2^l
 * returns newly allocated array of F128
 */
F128* eq_poly(const F128 *pt, size_t l, size_t *out_len)
{
    // The result length is 2^l
    size_t n = ((size_t)1 << l);
    *out_len = n;

    // allocate ret array
    F128 *ret = malloc(n * sizeof(F128));
    if (!ret) {
        fprintf(stderr, "eq_poly: out of memory\n");
        exit(1);
    }

    // 1) ret[0] = 1, and presumably set the rest to 0
    //    but we'll do that in the loop's logic. We'll at least zero them:
    for (size_t i = 0; i < n; i++) {
        ret[i].low  = 0ULL;
        ret[i].high = 0ULL;
    }
    // ret[0] = 1
    ret[0] = f128_one();

    // 2) for i in [0..l):
    for (size_t i = 0; i < l; i++) {
        // half = 1 << i
        size_t half = ((size_t)1 << i);

        // for j in [0..half):
        //   ret[j + half] = pt[i] * ret[j];
        //   ret[j] += ret[j + half];
        for (size_t j = 0; j < half; j++) {
            F128 temp = f128_mul(pt[i], ret[j]);
            ret[j + half] = temp;
            ret[j] = f128_add(ret[j], temp);
        }
    }

    return ret;
}