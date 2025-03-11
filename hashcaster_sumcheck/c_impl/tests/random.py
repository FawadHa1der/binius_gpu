#!/usr/bin/env python3
import re

# Regex to match a 128-bit decimal (or any decimal) number:
number_pattern = re.compile(r"\d+")

def split_128_to_64(num128: int):
    """
    Split a 128-bit unsigned integer into (low_64, high_64).
    """
    low = num128 & 0xFFFFFFFFFFFFFFFF
    high = (num128 >> 64) & 0xFFFFFFFFFFFFFFFF
    return (low, high)

def process_line(line: str) -> str:
    """
    Parse one line, replace:
      - '[' with '{'
      - ']' with '}'
      - each 128-bit decimal with "{low, high}"
    Keep commas, spaces, etc. in place.
    """
    i = 0
    length = len(line)
    out_tokens = []

    while i < length:
        ch = line[i]

        if ch == '[':
            out_tokens.append('{')
            i += 1
        elif ch == ']':
            out_tokens.append('}')
            i += 1
        else:
            # Maybe it's a number or something else
            # If it's a digit, parse the full number
            if ch.isdigit():
                # parse the full decimal
                match = number_pattern.match(line, i)
                if match:
                    num_str = match.group(0)
                    i_end = match.end()
                    # Convert to 128-bit int
                    val_128 = int(num_str)
                    low, high = split_128_to_64(val_128)
                    # produce a "C struct" snippet: "{low, high}"
                    out_tokens.append(f"{{0x{low:016x},0x{high:016x}}}")
                    i = i_end
                else:
                    # Shouldn't happen if ch.isdigit() matched
                    out_tokens.append(ch)
                    i += 1
            else:
                # Just copy the character
                out_tokens.append(ch)
                i += 1

    return "".join(out_tokens)

def main():
    input_file = "cobasis_frobenius_table.txt"
    output_file = "cobasis_frobenius_table_2_64_bits.txt"

    with open(input_file, "r") as fin, open(output_file, "w") as fout:
        for line in fin:
            line_out = process_line(line.rstrip("\n"))
            fout.write(line_out + "\n")

if __name__ == "__main__":
    main()