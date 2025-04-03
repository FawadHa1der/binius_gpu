#ifndef DEBUG_UTILS
#define DEBUG_UTILS

#include <stdio.h>

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>    // for seeding
#include <assert.h>
#include "field.h" // Assuming this is where F128 is defined


// Simple implementation of f128_to_string
// ONLY DEBUG purposes. Do not use in production code. Memry Leak here 
void f128_print(char* message, F128 val) ;

#endif