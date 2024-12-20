import re

# Function to rename the variables based on the rules
def rename_variables(c_code):
    # Replace x0 to x127 with x[0] to x[127]
    for i in range(32):
        c_code = re.sub(rf'\bx{i}\b', f'x[{i}]', c_code)

    # Replace x16 to x31 with y[0] to y[127]
    for i in range(32, 64):
        c_code = re.sub(rf'\bx{i}\b', f'y[{i - 32}]', c_code)

    # Replace y0 to 127 with z[0] to z[127]
    for i in range(128):
        c_code = re.sub(rf'\by{i}\b', f'z[{i}]', c_code)

    return c_code

# Load the C code from a file
def process_c_file(input_file, output_file):
    with open(input_file, 'r') as file:
        c_code = file.read()

    # Perform the renaming
    updated_code = rename_variables(c_code)

    # Save the updated code to a new file
    with open(output_file, 'w') as file:
        file.write(updated_code)

    print(f"Updated code has been saved to {output_file}")

# Example usage
input_file = 'bs_multiply_32.c'  # Path to the input C program
output_file = 'updated_bs_multiply_32.c'  # Path to the output C program
process_c_file(input_file, output_file)