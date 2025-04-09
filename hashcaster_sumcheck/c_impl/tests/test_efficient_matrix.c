#include "test_efficient_matrix.h"
#include "unity.h"
#include "../matrix_utils.h"
#include "time.h"
#include "test_helpers.h"
// Compare entire EfficientMatrix for equality
static int EfficientMatrix_equal(const EfficientMatrix* a, const EfficientMatrix* b){
    for(int i=0; i<EFFICIENT_MATRIX_SIZE; i++){
        if(!f128_eq(a->data[i], b->data[i])){
            return 0;
        }
    }
    return 1;
}

// A helper macro for Unity
#define TEST_ASSERT_EM_EQUAL_MESSAGE(a,b,msg) \
    TEST_ASSERT_TRUE_MESSAGE(EfficientMatrix_equal((a),(b)),(msg))

    
EfficientMatrix efficient_matrix_zero(void)
{
    EfficientMatrix em;
    for(int i=0; i<EFFICIENT_MATRIX_SIZE; i++){
        em.data[i]= f128_zero();
    }
    return em;
}

// ---------------------------------------------------------------------
// 1) test_from_cols_all_zeros
void test_from_cols_all_zeros(void)
{
    // Rust code snippet:
    // let cols = [BinaryField128b::ZERO; 128];
    // let matrix = EfficientMatrix::from_cols(&cols);
    // let expected= EfficientMatrix([ZERO;256*16]);
    // assert_eq!(matrix, expected)

    // Create an all-zero input array of length 128
    F128 cols[128];
    for(int i=0; i<128; i++){
        cols[i] = f128_zero();
    }

    // Call from_cols
    EfficientMatrix matrix = from_cols(cols);

    // Manually compute the expected result => all zero
    EfficientMatrix expected;
    for(int i=0; i<EFFICIENT_MATRIX_SIZE; i++){
        expected.data[i] = f128_zero();
    }

    // Compare
    TEST_ASSERT_MESSAGE(
        EfficientMatrix_equal(&matrix, &expected),
        "Precomputed sums for all-zero input should be zero."
    );
}

// 2) test_from_cols_single_non_zero_column
void test_from_cols_single_non_zero_column(void)
{
    // Rust snippet:
    // let mut cols=[ZERO;128]; cols[0]= from(1);
    // let matrix= from_cols(&cols);
    // let mut expected_precomp=[ZERO;256*16];
    // for i in [0..256]:
    //   if i is even => zero, else => col[0]
    // let expected= EfficientMatrix(expected_precomp);
    // assert_eq!(matrix, expected)

    // 1) Create an array of 128 elements, all zero except the first
    F128 cols[128];
    for(int i=0; i<128; i++){
        cols[i] = f128_zero();
    }
    cols[0] = f128_from_uint64(1); // single non-zero column

    // 2) Call from_cols
    EfficientMatrix matrix = from_cols(cols);

    // 3) Manually build expected data. For indices [0..256) in the first chunk (the rest are 0):
    //    if i is even => zero, else => col[0].
    //    This pattern repeats for each group => but from your snippet, it sets the first chunk specifically.
    // We'll replicate the approach: fill all 4096 with zero, then for the first 256 we set odd indices to col[0].
    EfficientMatrix expected;
    for(int i=0; i<EFFICIENT_MATRIX_SIZE; i++){
        expected.data[i] = f128_zero();
    }

    // For the first 256 entries, if i%2!=0 => set col[0].
    for(int i=0; i<256; i++){
        if((i % 2) != 0){
            expected.data[i] = cols[0];
        }
    }

    // 4) Compare
    TEST_ASSERT_MESSAGE(
        EfficientMatrix_equal(&matrix, &expected),
        "Precomputed sums for a single non-zero column are incorrect."
    );
}
void test_from_cols_with_queries(void)
{
    // 1) Initialize some non-zero columns
    F128 c0 = f128_from_uint64(568ULL);
    F128 c1 = f128_from_uint64(123ULL);
    F128 c2 = f128_from_uint64(456ULL);
    F128 c3 = f128_from_uint64(789ULL);
    F128 c4 = f128_from_uint64(1011ULL);

    // 2) Create a 128-element array of columns, all zero except c0..c4
    F128 cols[128];
    for(int i=0;i<128;i++){
        cols[i]= f128_zero();
    }
    cols[0] = c0;
    cols[1] = c1;
    cols[2] = c2;
    cols[3] = c3;
    cols[4] = c4;

    // 3) Build the matrix
    EfficientMatrix matrix = from_cols(cols);

    // We'll define a macro or function to check equality quickly
    #define CHECK_CELL(index, expectedVal)                                           \
        TEST_ASSERT_TRUE_MESSAGE(                                                   \
            f128_eq(matrix.data[(index)], (expectedVal)),                      \
            "matrix[" #index "] does not match " #expectedVal );

    // 4) Assertions that replicate the Rust lines:
    //    assert_eq!(matrix[0b00000001 as usize], c0);
    //    etc.

    // Indices in decimal: 0b00000001 => 1, 0b00000010 =>2, 0b00000011 =>3, etc.
    // We'll define them carefully:
    CHECK_CELL(1, c0);
    CHECK_CELL(2, c1);
    {
        // c0 + c1
        F128 sum = f128_add(c0, c1);
        CHECK_CELL(3, sum);
    }
    CHECK_CELL(4, c2);
    {
        // c0 + c2
        F128 sum = f128_add(c0, c2);
        CHECK_CELL(5, sum);
    }
    {
        // c1 + c2
        F128 sum = f128_add(c1, c2);
        CHECK_CELL(6, sum);
    }
    {
        // c0 + c1 + c2
        F128 sum = f128_add(f128_add(c0, c1), c2);
        CHECK_CELL(7, sum);
    }
    CHECK_CELL(8, c3);
    {
        // c0 + c3
        F128 sum = f128_add(c0, c3);
        CHECK_CELL(9, sum);
    }
    {
        // c1 + c3
        F128 sum = f128_add(c1, c3);
        CHECK_CELL(10, sum);
    }
    {
        // c0 + c1 + c3
        F128 sum = f128_add(f128_add(c0, c1), c3);
        CHECK_CELL(11, sum);
    }
    {
        // c2 + c3
        F128 sum = f128_add(c2, c3);
        CHECK_CELL(12, sum);
    }
    {
        // c0 + c2 + c3
        F128 sum = f128_add(f128_add(c0, c2), c3);
        CHECK_CELL(13, sum);
    }
    {
        // c1 + c2 + c3
        F128 sum = f128_add(f128_add(c1, c2), c3);
        CHECK_CELL(14, sum);
    }
    {
        // c0 + c1 + c2 + c3
        F128 sum = f128_add(f128_add(c0, c1), f128_add(c2, c3));
        CHECK_CELL(15, sum);
    }
    CHECK_CELL(16, c4);
    {
        // c0 + c4
        F128 sum = f128_add(c0, c4);
        CHECK_CELL(17, sum);
    }
    {
        // c1 + c4
        F128 sum = f128_add(c1, c4);
        CHECK_CELL(18, sum);
    }
    {
        // c0 + c1 + c4
        F128 sum = f128_add(f128_add(c0, c1), c4);
        CHECK_CELL(19, sum);
    }
    {
        // c2 + c4
        F128 sum = f128_add(c2, c4);
        CHECK_CELL(20, sum);
    }
    {
        // c0 + c2 + c4
        F128 sum = f128_add(f128_add(c0, c2), c4);
        CHECK_CELL(21, sum);
    }
    {
        // c1 + c2 + c4
        F128 sum = f128_add(f128_add(c1, c2), c4);
        CHECK_CELL(22, sum);
    }
    {
        // c0 + c1 + c2 + c4
        F128 sum = f128_add(f128_add(c0, c1), f128_add(c2, c4));
        CHECK_CELL(23, sum);
    }
    {
        // c3 + c4 => index= 24 decimal => 0b00011000
        F128 sum = f128_add(c3, c4);
        CHECK_CELL(24, sum);
    }

    // That's it. The checks replicate the Rust test lines.
    #undef CHECK_CELL
}


void test_from_rows_all_zeros(void)
{
    // 1) rows= all zeros
    F128 rows[128];
    for(int i=0;i<128;i++){
        rows[i] = f128_zero();
    }

    // 2) from_rows
    EfficientMatrix matrix = from_rows(rows);

    // 3) expected => from an all-zero input => columns also all zero => from_cols(all zero columns)
    // or build a direct "all zero" matrix if you know the final data layout:
    F128 zero_cols[128];
    for(int i=0;i<128;i++){
        zero_cols[i]= f128_zero();
    }
    EfficientMatrix expected= from_cols(zero_cols);

    // 4) compare
    TEST_ASSERT_EM_EQUAL_MESSAGE(&matrix, &expected, "Resulting matrix should be all zeros.");
}


// test_from_rows_all_ones
void test_from_rows_all_ones(void)
{
    // rust:
    // let rows=[1;128], from_rows => expect the first column is u128::MAX (in the snippet)
    // We'll replicate:

    F128 rows[128];
    // If "from(1)" just means lo=1, hi=0, let's define that:
    F128 one= f128_from_uint64(1);
    for(int i=0;i<128;i++){
        rows[i]= one;
    }

    // from_rows => matrix
    EfficientMatrix matrix= from_rows(rows);

    // Then we define "cols" => all zero except the first column is u128::MAX
    // We'll assume we have "F128_from_uint128(...)" or do the best we can:
    // for an all-ones 128-bit pattern => 0xFFFF_FFFF... => that might require a bigger approach
    F128 cols[128];
    for(int i=0;i<128;i++){
        cols[i]= f128_zero();
    }
    // If we define "all bits => 1" => let allones = 0xFFFF_FFFF_FFFF_FFFF... => for demonstration:
    // We'll define a 128-bit pattern: (uint64_t)-1 for lo & hi if you want:
    F128 allones = {(uint64_t)-1, (uint64_t)-1};
    cols[0]= allones;

    // from_cols => expected
    EfficientMatrix expected= from_cols(cols);

    // compare
    TEST_ASSERT_EM_EQUAL_MESSAGE(&matrix, &expected, "Resulting matrix does not match the expected result.");
}


void test_from_rows_single_element_in_first_row(void)
{
    // rust:
    // rows= all zero except rows[0]=1 => from_rows => expect col[0]=1
    F128 rows[128];
    for(int i=0;i<128;i++){
        rows[i]= f128_zero();
    }
    rows[0]= f128_from_uint64(1);

    // from_rows
    EfficientMatrix matrix= from_rows(rows);

    // expected => only the first column has first bit => so col[0]= 1
    F128 cols[128];
    for(int i=0;i<128;i++){
        cols[i]= f128_zero();
    }
    cols[0]= f128_from_uint64(1);

    EfficientMatrix expected= from_cols(cols);

    TEST_ASSERT_EM_EQUAL_MESSAGE(&matrix, &expected,
        "Resulting matrix does not match expected with single non-zero element in the first row.");
}

// test_from_rows_alternate_rows_full_128_bits
void test_from_rows_alternate_rows_full_128_bits(void)
{
    // rust:
    // even_pattern= 0xAAAA_AAAA... odd_pattern= 0x55555555...
    // rows[i]= even if i%2==0, odd if i%2!=0
    // then from_rows => compare with from_cols of the same pattern transposed in columns

    // We'll define the patterns as 64-bit halves:
    // 0xAAAA_AAAA_AAAA_AAAA => 0xAAAA_AAAA_AAAA_AAAA for the lo half, maybe hi= the same
    // We'll do a single 128-bit approach:
    F128 rows[128];

    // define "F128 pattern" if you have a "F128_from_uint128(...)"
    // or we do "F128 evenF= {0xAAAA_AAAA_AAAA_AAAA,0xAAAA_AAAA_AAAA_AAAA};"
    // for demonstration, let's define them as:
    F128 evenF = {(uint64_t)0xAAAAAAAAAAAAAAAAULL, (uint64_t)0xAAAAAAAAAAAAAAAAULL};
    F128 oddF  = {(uint64_t)0x5555555555555555ULL, (uint64_t)0x5555555555555555ULL};

    for(int i=0;i<128;i++){
        if((i%2)==0){
            rows[i]= evenF;
        } else {
            rows[i]= oddF;
        }
    }

    // from_rows => matrix
    EfficientMatrix matrix= from_rows(rows);

    // Then define columns => i%2==0 => evenF, else oddF
    F128 cols[128];
    for(int i=0;i<128;i++){
        if((i%2)==0){
            cols[i]= evenF;
        } else {
            cols[i]= oddF;
        }
    }

    EfficientMatrix expected= from_cols(cols);

    TEST_ASSERT_EM_EQUAL_MESSAGE(&matrix,&expected,
        "Resulting matrix does not match expected with full 128-bit alternating row patterns.");
}

// test_from_rows_all_rows_max
void test_from_rows_all_rows_max(void)
{
    // rust: rows= [u128::MAX;128]
    // define a "F128 allones= { -1ULL, -1ULL }"
    F128 allones= {(uint64_t)-1, (uint64_t)-1};
    F128 rows[128];
    for(int i=0;i<128;i++){
        rows[i]= allones;
    }

    // from_rows => matrix
    EfficientMatrix matrix= from_rows(rows);

    // define columns => allones for each => from_cols => expected
    F128 cols[128];
    for(int i=0;i<128;i++){
        cols[i]= allones;
    }
    EfficientMatrix expected= from_cols(cols);

    TEST_ASSERT_EM_EQUAL_MESSAGE(
        &matrix, &expected,
        "Resulting matrix does not match the expected result (all rows = max)."
    );
}

void test_apply_all_zeros(void)
{
    // "Create an EfficientMatrix filled with zeros => default()"
    EfficientMatrix matrix= efficient_matrix_zero();

    // "Input vector filled with zeros => F128_ZERO"
    F128 input= f128_zero();

    // "Applying the matrix => result => matrix.apply(input)"
    F128 result= efficient_matrix_apply(&matrix, input);

    // "Expected => zero"
    TEST_ASSERT_F128_EQUAL_MESSAGE(result, f128_zero(),
        "Applying all-zero matrix to zero vector should result in zero.");
}

// =========== 2) test_apply_single_nonzero_byte ===========
void test_apply_single_nonzero_byte(void)
{
    // "Create an EfficientMatrix with a single non-zero value => matrix_data[1]= from(42)"
    // We'll do a local array of 4096 F128:
    EfficientMatrix matrix;
    for(int i=0; i<EFFICIENT_MATRIX_SIZE; i++){
        matrix.data[i]= f128_zero();
    }
    // set data[1]= from(42)
    matrix.data[1]= f128_from_uint64(42ULL);

    // "Input vector with the first byte set to 1 => i.e. input= from(1)"
    F128 input= f128_from_uint64(1ULL);

    // apply => expect result= matrix.data[1]
    F128 result= efficient_matrix_apply(&matrix, input);

    // compare
    F128 expected= f128_from_uint64(42ULL);
    TEST_ASSERT_F128_EQUAL_MESSAGE(result, expected,
        "Applying matrix with single non-zero entry should return that entry.");
}

// =========== 3) test_apply_against_traditional_matrix ====
void test_apply_against_traditional_matrix(void)
{
    // "Setup random columns of u128 => array of 128"
    // We'll do a simplistic approach
    srand((unsigned)time(NULL));
    uint64_t cols[128];
    for(int i=0; i<128; i++){
        // random 64-bit. Real code might do 128. We'll keep it short
        cols[i]= (uint64_t)rand();
    }

    // "Create an EfficientMatrix from the columns"
    F128 fCols[128];
    for(int i=0;i<128;i++){
        // interpret cols[i] as a small F128
        fCols[i]= f128_from_uint64(cols[i]);
    }
    EfficientMatrix matrix= from_cols(fCols);

    // "Create a traditional matrix from the columns => Matrix::new(cols)"
    Matrix traditional = matrix_new(fCols);

    // "Generate a random input vector => 64 bits"
    F128 inputField= f128_rand();

    F128 result_efficient= efficient_matrix_apply(&matrix, inputField);

    // "Apply the traditional matrix => result_traditional= matrix.apply(rhs)"
    F128 expected= matrix_apply(&traditional, inputField);


    TEST_ASSERT_F128_EQUAL_MESSAGE(result_efficient, expected,
        "EfficientMatrix and traditional matrix gave different results");
}

// =========== 4) test_frobenius_inv_lc_all_zeros ==========
void test_frobenius_inv_lc_all_zeros(void)
{
    // "Initialize gammas with all zeros => gammas[128]"
    F128 gammas[128];
    for(int i=0;i<128;i++){
        gammas[i]= f128_zero();
    }

    // "matrix= from_frobenius_inv_lc(&gammas)"
    EfficientMatrix matrix= from_frobenius_inv_lc(gammas);

    // "expected= EfficientMatrix::default() => all zero"
    EfficientMatrix expected= efficient_matrix_zero();

    // "assert eq"
    // We'll define a helper or do a loop compare
    for(int i=0;i<EFFICIENT_MATRIX_SIZE;i++){
        TEST_ASSERT_F128_EQUAL_MESSAGE(matrix.data[i], expected.data[i],
            "Matrix should be all zeros when gammas are zero.");
    }
}

// =========== 5) test_frobenius_inv_lc_single_gamma =======
void test_frobenius_inv_lc_single_gamma(void)
{
    // "Initialize gammas => all zero except gamma[0]= from(5467)"
    F128 gammas[128];
    for(int i=0;i<128;i++){
        gammas[i]= f128_zero();
    }
    // gamma[0]= 5467
    F128 val5467= f128_from_uint64(5467ULL);
    gammas[0]= val5467;

    // "matrix= from_frobenius_inv_lc(&gammas)"
    EfficientMatrix matrix= from_frobenius_inv_lc(gammas);
    F128 expected_cols[128];
    for(int j=0;j<128;j++){
        F128 frobVal= FROBENIUS[0][j];
        F128 partial= f128_add(val5467, frobVal); // WAIT, the snippet says multiply => actually " * "
        // correct logic => partial= val5467 * frobVal
        // We'll define a F128 F128_mul(...) if not existing
        partial= f128_mul(val5467, frobVal);

        expected_cols[j]= partial;
    }
    // "expected= from_cols(&expected_cols)"
    EfficientMatrix expected= from_cols(expected_cols);

    // compare
    for(int i=0;i<EFFICIENT_MATRIX_SIZE;i++){
        TEST_ASSERT_F128_EQUAL_MESSAGE(matrix.data[i], expected.data[i],
            "Matrix should match the expected result for a single non-zero gamma.");
    }
}

// =========== 6) test_frobenius_inv_lc_multiple_nonzero_gammas ===
void test_frobenius_inv_lc_multiple_nonzero_gammas(void)
{
    // "gammas= [0..] => gammas[0]= from(383), gammas[1]= from(463)"
    F128 gammas[128];
    for(int i=0;i<128;i++){
        gammas[i]= f128_zero();
    }
    F128 val383= f128_from_uint64(383ULL);
    F128 val463= f128_from_uint64(463ULL);
    gammas[0]= val383;
    gammas[1]= val463;

    // "matrix= from_frobenius_inv_lc(&gammas)"
    EfficientMatrix matrix= from_frobenius_inv_lc(gammas);

    F128 expected_cols[128];
    for(int j=0;j<128;j++){
        F128 frob0= FROBENIUS[0][j];
        F128 frob127= FROBENIUS[127][j];

        F128 partial1= f128_mul(val383, frob0);
        F128 partial2= f128_mul(val463, frob127);
        F128 sum= f128_add(partial1, partial2);

        expected_cols[j]= sum;
    }

    EfficientMatrix expected= from_cols(expected_cols);

    // compare
    for(int i=0;i<EFFICIENT_MATRIX_SIZE;i++){
        TEST_ASSERT_F128_EQUAL_MESSAGE(matrix.data[i], expected.data[i],
            "Matrix should match the expected result for multiple non-zero gammas.");
    }
}

// =========== 7) test_frobenius_inv_lc_wraparound_indices ===
void test_frobenius_inv_lc_wraparound_indices(void)
{
    // "gammas[127]= from(5363)"
    F128 gammas[128];
    for(int i=0;i<128;i++){
        gammas[i]= f128_zero();
    }
    F128 val5363= f128_from_uint64(5363ULL);
    gammas[127]= val5363;

    // "matrix= from_frobenius_inv_lc(&gammas)"
    EfficientMatrix matrix= from_frobenius_inv_lc(gammas);

    F128 expected_cols[128];
    for(int j=0;j<128;j++){
        F128 frobVal= FROBENIUS[1][j];
        F128 product= f128_mul(val5363, frobVal);
        expected_cols[j]= product;
    }

    EfficientMatrix expected= from_cols(expected_cols);

    // compare
    for(int i=0;i<EFFICIENT_MATRIX_SIZE;i++){
        TEST_ASSERT_F128_EQUAL_MESSAGE(matrix.data[i], expected.data[i],
            "Matrix should match the expected result for wraparound indices.");
    }
}

