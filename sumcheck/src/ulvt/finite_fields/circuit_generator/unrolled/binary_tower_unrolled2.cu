//This file is auto generated by multiply_and_generate_circuit.cpp
#include <cstdint>
#include "binary_tower_unrolled.cuh"
template<>
__host__ __device__ void multiply_unrolled<2>(const uint32_t *field_element_a, const uint32_t *field_element_b, uint32_t *destination){
uint32_t v1;
uint32_t v2;
uint32_t v3;
uint32_t v4;
uint32_t v5;
uint32_t v6;
uint32_t v7;
uint32_t v8;
uint32_t v9;
uint32_t v10;
uint32_t v11;
uint32_t v12;
uint32_t v13;
uint32_t v14;
uint32_t v15;
uint32_t v16;
uint32_t v17;
uint32_t v18;
uint32_t v19;
uint32_t v20;
uint32_t v21;
uint32_t v22;
uint32_t v23;
uint32_t v24;
uint32_t v25;
uint32_t v26;
uint32_t v27;

v13 = field_element_a[0];
v17 = field_element_b[0];
v15 = field_element_a[1];
v19 = field_element_b[1];
v6 = field_element_a[2];
v8 = field_element_b[2];
v2 = field_element_a[3];
v3 = field_element_b[3];

v1 = v2 & v3;
v4 = v1;
v5 = v6 ^ v2;
v7 = v8 ^ v3;
v1 ^= v6 & v8;
v9 = v1;
v4 ^= v1;
v4 ^= v5 & v7;
v10 = v4;
v11 = v9;
v11 ^= v4;
v12 = v13 ^ v6;
v14 = v15 ^ v2;
v16 = v17 ^ v8;
v18 = v19 ^ v3;
v20 = v15 & v19;
v4 ^= v20;
v21 = v13 ^ v15;
v22 = v17 ^ v19;
v20 ^= v13 & v17;
v9 ^= v20;
v4 ^= v20;
v4 ^= v21 & v22;
v23 = v9;
v24 = v4;
v10 ^= v9;
v11 ^= v4;
v25 = v14 & v18;
v11 ^= v25;
v26 = v12 ^ v14;
v27 = v16 ^ v18;
v25 ^= v12 & v16;
v10 ^= v25;
v11 ^= v25;
v11 ^= v26 & v27;

destination[0] = v23;
destination[1] = v24;
destination[2] = v10;
destination[3] = v11;

}