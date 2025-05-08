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
    # input_file = "cobasis_frobenius_table.txt"
    # output_file = "cobasis_frobenius_table_2_64_bits.txt"

    # with open(input_file, "r") as fin, open(output_file, "w") as fout:
    #     for line in fin:
    #         line_out = process_line(line.rstrip("\n"))
    #         fout.write(line_out + "\n")


    # convert the following list in 128 bit numbers to 2 64 bit numbers
    #        286199402842439698884600957666611691959ULL,
        # 92152884645276120751159248850998264247ULL,
        # 194088056572031856754469952786247188481ULL,
        # 299076299051606071403356588563077530026ULL,
        # 152881988702699464694451933917556507050ULL,
        # 194088056572031856754469952786247188481ULL,
        # 72193695521068243347088020961476214811ULL,
        # 72193695521068243347088020961476214811ULL,
        # 0ULL,
        # 18692268690725379462709786785192346071ULL,
        # 207463413279617572725564511330318156247ULL,
        # 194088056572031856754469952786247188481ULL,
        # 288774782084272973388352083845904859232ULL,
        # 100045175870249058746525603271412809824ULL,
        # 194088056572031856754469952786247188481ULL,
        # 286199402842439698884600957666611691945ULL,
        # 286199402842439698884600957666611691945ULL,
        # 0ULL,
        # 288774782084272973388352083845904859244ULL,
        # 288774782084272973388352083845904859244ULL,
        # 0ULL,
        # 74769074762901517850839147140769382878ULL,
        # 74769074762901517850839147140769382878ULL,
        # 0ULL,
        # 299076299051606071403356588563077530034ULL,
        # 299076299051606071403356588563077530034ULL,
        # 0ULL

    list = [ 286199402842439698884600957666611691959, 
             92152884645276120751159248850998264247,
             194088056572031856754469952786247188481,
             299076299051606071403356588563077530026,
             152881988702699464694451933917556507050,
             194088056572031856754469952786247188481,
             72193695521068243347088020961476214811,
             72193695521068243347088020961476214811,
             0,
             18692268690725379462709786785192346071,
             207463413279617572725564511330318156247,
             194088056572031856754469952786247188481,
             288774782084272973388352083845904859232,
             100045175870249058746525603271412809824,
             194088056572031856754469952786247188481,
             286199402842439698884600957666611691945,
             286199402842439698884600957666611691945,
             0,
             288774782084272973388352083845904859244,
             288774782084272973388352083845904859244,
             0,
             74769074762901517850839147140769382878,
             74769074762901517850839147140769382878,
             0,
             299076299051606071403356588563077530034,
             299076299051606071403356588563077530034,
             0]
    
    # convert the list to 2 64 bit numbers
    for i in list:
        # original = i
        print(f"original decimal  {i}")
        low, high = split_128_to_64(i)
        print(f"Hex 0x{low:016x}, 0x{high:016x}")
        # print in dec
        print(f"Decimal {low}, {high}")
    # low, high = split_128_to_64(9838263505978427528)
    # print(f"0x{low:016x}, 0x{high:016x}")
    # # print in dec
    # print(f"{low}, {high}")
    # Test the process_line function
if __name__ == "__main__":
    main()