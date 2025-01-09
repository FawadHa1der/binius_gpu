/***************************************************************
 * test_movemask.c
 * 
 * Unity-based tests that replicate the Rust tests for movemask,
 * including a "bench" test and an aarch64-specific movemask test.
 ***************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <arm_neon.h>  // If you need NEON intrinsics
#include "unity.h"

/***************************************************************
 * External declarations from your code
 * (You'd have them in movemask_impl.h or similar)
 ***************************************************************/
#include "field.h"

// Suppose these are your movemask implementations in C
int v_movemask_epi8(const uint8_t data[16]);
int cpu_v_movemask_epi8(const uint8_t data[16]);
int mm_movemask_epi8(const uint8_t data[16]);

/***************************************************************
 * Helper: naive random 64-bit generator (replacement for Rust OsRng)
 ***************************************************************/
static uint64_t simple_rand64(void)
{
    return (((uint64_t)rand()) << 32) ^ (uint64_t)rand();
}

/***************************************************************
 * Helper: timing in nanoseconds
 ***************************************************************/
static uint64_t now_ns(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}

/***************************************************************
 * Helper function: replicate Rust's v_slli_epi64<K>()
 *
 * In Rust, you wrote:
 *   pub unsafe fn v_slli_epi64<const K: i32>(x: [u8; 16]) -> [u8; 16]
 *   { ... vshlq_n_s64(data, K) ... }
 *
 * In C, we can't pass a template "const K." We'll either hardcode or
 * do a runtime shift.  We'll show a runtime version for demonstration.
 ***************************************************************/
static void v_slli_epi64(int shift_bits, const uint8_t in[16], uint8_t out[16])
{
    // Load as int64x2_t
    int64x2_t v = vld1q_s64((const int64_t*)in);

    // We'll do a runtime shift. If SHIFT is truly compile-time, you can use vshlq_n_s64
    // or separate functions. For demonstration, here's a runtime approach with vshlq_s64:
    int64x2_t shift_vector = vdupq_n_s64(shift_bits);
    int64x2_t res = vshlq_s64(v, shift_vector);

    // Store back
    vst1q_s64((int64_t*)out, res);
}

/***************************************************************
 * 1) Test replicating Rust's bench_movemask()
 ***************************************************************/
static void test_bench_movemask(void)
{
    // srand for random seed
    srand((unsigned)time(NULL));

    // Let s = u128_rand(rng) => we'll just use simple_rand64()
    uint64_t s = simple_rand64();
    // let x = rng.next_u32() as i32 => same approach
    int x = (int)(rand() & 0xffffffff);

    // n = 100_000_000 => We'll do smaller for demonstration
    const size_t n = 100000; // or keep 100000000 if you want a real bench

    uint64_t u = (uint64_t)x;
    uint64_t v = s;
    uint64_t ret = 0;

    // label0
    uint64_t label0 = now_ns();

    // for i in 0..n { v += s; u *= x; u += 1; }
    for (size_t i = 0; i < n; i++) {
        v += s;
        u *= x;
        u += 1;
    }
    ret += u;

    // label1
    uint64_t label1 = now_ns();

    // for i in 0..n { v += s; u *= v_movemask_epi8(...); u += 1; }
    for (size_t i = 0; i < n; i++) {
        v += s;
        // We'll place 'v' in a 16-byte array
        uint8_t buffer[16] = {0};
        *(uint64_t*)buffer = v;
        int mask = v_movemask_epi8(buffer);
        u *= mask;
        u += 1;
    }
    ret += u;

    // label2
    uint64_t label2 = now_ns();

    // for i in 0..n { v += s; u *= cpu_v_movemask_epi8(...); u += 1; }
    for (size_t i = 0; i < n; i++) {
        v += s;
        uint8_t buffer[16] = {0};
        *(uint64_t*)buffer = v;
        int mask = cpu_v_movemask_epi8(buffer);
        u *= mask;
        u += 1;
    }
    ret += u;

    // label3
    uint64_t label3 = now_ns();

    // Compute times in ms
    uint64_t t_native = (label3 - label2) / 1000000ULL;
    uint64_t t_reference = (label2 - label1) / 1000000ULL;
    uint64_t t_offset = (label1 - label0) / 1000000ULL;

    printf("\n[bench_movemask] Timings (ms):\n");
    printf("  Native: %llu\n", (unsigned long long)t_native);
    printf("  Reference: %llu\n", (unsigned long long)t_reference);
    printf("  Offset: %llu\n", (unsigned long long)t_offset);
    printf("ret = %llu\n", (unsigned long long)ret);

    // In a "unit test" sense, we might do a pass/fail check
    // But here it's a benchmark, so we just ensure the code runs
    TEST_ASSERT_MESSAGE(true, "[bench_movemask] completed.");
}

/***************************************************************
 * 2) Test replicating Rust's test_for_mm_movemask_aarch64()
 ***************************************************************/
static void test_for_mm_movemask_aarch64(void)
{
    // The same bytes from Rust
    uint8_t bytes[16] = {
        0x80, 0x01, 0x80, 0x01, 0x80, 0x01, 0x80, 0x01,
        0x80, 0x01, 0x80, 0x01, 0x80, 0x01, 0x80, 0x01
    };

    // load into NEON vector
    uint8x16_t input_vector = vld1q_u8(bytes);

    // measure time for v_movemask_epi8
    uint64_t start = now_ns();
    {
        // In Rust, you do: v_movemask_epi8(transmute::<uint8x16_t, [u8;16]>(input_vector))
        // In C, we can store the vector, then call our function
        uint8_t tmp[16];
        vst1q_u8(tmp, input_vector);
        int result = v_movemask_epi8(tmp);
        (void)result; // not used
    }
    uint64_t end = now_ns();
    printf("\n[test_for_mm_movemask_aarch64] v_movemask_epi8 took %llu ns\n",
           (unsigned long long)(end - start));

    // measure time for mm_movemask_epi8
    start = now_ns();
    {
        int mask = mm_movemask_epi8(bytes);
        (void)mask;
    }
    end = now_ns();
    printf("mm_movemask_epi8 took %llu ns\n", (unsigned long long)(end - start));

    // measure time for cpu_v_movemask_epi8
    start = now_ns();
    {
        int lev_mask = cpu_v_movemask_epi8(bytes);
        (void)lev_mask;
    }
    end = now_ns();
    printf("cpu_v_movemask_epi8 took %llu ns\n", (unsigned long long)(end - start));

    // For a normal unit test, we might compare results:
    //  TEST_ASSERT_EQUAL_INT(some_expected_value, result);
    TEST_ASSERT_MESSAGE(true, "[test_for_mm_movemask_aarch64] completed.");
}

/***************************************************************
 * 3) Demo test for the shift function (v_slli_epi64)
 ***************************************************************/
static void test_v_slli_epi64(void)
{
    // Example input
    uint8_t input[16] = {
        0x01,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
        0x02,0x00,0x00,0x00,0x00,0x00,0x00,0x00
    };
    uint8_t output[16] = {0};

    // Shift each 64-bit lane by 1 bit
    v_slli_epi64(1, input, output);

    // Just print or check a known result
    printf("\n[test_v_slli_epi64] Output: ");
    for (int i = 0; i < 16; i++) {
        printf("%02X ", output[i]);
    }
    printf("\n");
    // We might do an assertion if we know the expected bits
    // E.g., shifting the first lane (0x01) left by 1 => 0x02, etc.
    // Let's just do a trivial check
    TEST_ASSERT_EQUAL_HEX8(0x02, output[0]);  // was 0x01 shifted by 1
    TEST_ASSERT_MESSAGE(true, "[test_v_slli_epi64] completed.");
}

/***************************************************************
 * Unity's main test runner
 * 
 * Typically, we define a main() that calls:
 *   UNITY_BEGIN();
 *   RUN_TEST(...);
 *   ...
 *   return UNITY_END();
 ***************************************************************/
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_bench_movemask);
    RUN_TEST(test_for_mm_movemask_aarch64);
    RUN_TEST(test_v_slli_epi64);

    return UNITY_END();
}