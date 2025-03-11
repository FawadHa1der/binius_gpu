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
SOURCES = field.c matrix_utils.c  gf128_polyval.c

# Object files are derived from the source files:
OBJECTS = $(SOURCES:.c=.o)

# Your test file & implementation
TEST_SRC =  field.c matrix_utils.c gf128_polyval.c field.h matrix_utils.h 

TEST_EXE = hashcaster

test: $(TEST_EXE)

$(TEST_EXE):  $(SOURCES)
	$(CC) $(CFLAGS) -I  $(SOURCES) -o $@


# The final executable name:
TARGET  = mul_128_test

# The default make rule:
all: $(TARGET)

# Link rule: combine all object files into one executable
$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS)

# Compile rule: how to compile .c into .o
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up intermediates and the final binary
clean:
	# rm -f $(OBJECTS) $(TARGET)