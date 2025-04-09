# Makefile

# The compiler we want to use
CC      = clang

# Compiler flags:
#  -O3       : Optimization level 3
#  -g        : Generate debug info
#  -Wall     : Enable most compiler warnings
#  -Wextra   : Enable extra warnings
#  -std=c11  : Use C11 standard
#  -march=armv8-a+crypto : Example for enabling ARMv8 + crypto instructions
#  -mfpu=crypto-neon-fp-armv8 : Example if needed for NEON + crypto
CFLAGS  = -O3 -g -Wall -Wextra -std=c11 -arch arm64

# List your source files here:
SOURCES = ../field.c ../field.h ../matrix_utils.c ../matrix_utils.h tests_main.c test_field.c ../mle_poly.c ../debug_utils.c ../evaluations.c ../compressed_poly.c ../univariate_poly.c ../efficient_matrix.c test_efficient_matrix.c


# Object files are derived from the source files:
OBJECTS = $(filter %.o,$(SOURCES:.c=.o))

# Unity sources
UNITY_DIR = ./
UNITY_SRC = $(UNITY_DIR)/unity.c 

# Your test file & implementation
TEST_SRC = test_field.c ../field.c ../matrix_utils.c tests_main.c test_matrix_utils.c ../gf128_polyval.c test_mle.c ../mle_poly.c ../debug_utils.c ../evaluations.c test_evaluations.c ../compressed_poly.c test_compressed_poly.c  ../univariate_poly.c test_univariate_poly.c ../efficient_matrix.c test_efficient_matrix.c


TEST_EXE = run_tests

test: $(TEST_EXE)

$(TEST_EXE): $(UNITY_SRC) $(TEST_SRC)
	$(CC) $(CFLAGS) -I$(UNITY_DIR) $(UNITY_SRC) $(TEST_SRC) -o $@


# Clean up intermediates and the final binary
clean:
	 rm -rf $(TEST_EXE) $(TEST_EXE).dSYM