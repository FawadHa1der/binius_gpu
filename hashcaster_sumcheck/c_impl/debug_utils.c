
#include <stdio.h>

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>    // for seeding
#include <assert.h>
#include "debug_utils.h" // Assuming this is where F128 is defined


// Simple implementation of f128_to_string
// ONLY DEBUG purposes. Do not use in production code. Memry Leak here 
void f128_print(char* message, F128 val) {
    // Assuming F128 is a struct with two uint64_t members: high and low
    // Format the string as "F128(0x<high>, 0x<low>)"
    // Note: Adjust the format specifier as per your F128 struct definition
    // snprintf is used to safely format the string

    printf( " %s F128(0x%016llx, 0x%016llx)", message, val.high, val.low);
}

