#ifndef TEST_HELPERS_H
#define TEST_HELPERS_H
#include <stdint.h>
#include <stdbool.h>
#include <stddef.h> // for size_t
#include "unity.h"
#include "../field.h"


static inline void TEST_ASSERT_F128_EQUAL_MESSAGE(F128 a, F128 b, const char*msg){
    TEST_ASSERT_TRUE_MESSAGE(f128_eq(a,b), msg);
}


static inline void TEST_ASSERT_F128_ARRAY_EQUAL(
    const F128 *actual, const F128 *expected, size_t length, const char *msgPrefix)
{
    for(size_t i=0; i<length; i++){
        if(!f128_eq(actual[i], expected[i])){
            char msg[128];
            snprintf(msg, sizeof(msg), "%s at index %zu mismatch", msgPrefix, i);
            TEST_FAIL_MESSAGE(msg);
        }
    }
}

#endif