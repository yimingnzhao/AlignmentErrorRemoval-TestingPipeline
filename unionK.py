import sys
import os


def getIntersection( small, large ):
    
    # Generates the intersection of the small and large counts
    result = ""
    for i in range(len(large)):
        if i < len(small):
            small_char = small[i]
        else:
            small_char = '.'
        large_char = large[i]
        if small_char == 'X' or large_char == 'X':
            result += ('X')
        else:
            result += (large_char)
    return result

    


USAGE = "python unionK.py [fasta file 1] [fasta file 2]"
DESCRIPTION = "Takes the union of the errors ('X') in the fasta files"

if not len(sys.argv) == 3:
    print()
    print("\tError: Incorrect number of parameters\n")
    print("\tUSAGE: " + USAGE)
    print()
    print("\tDESCRIPTION: " + DESCRIPTION)
    print()
    sys.exit()


file1 = sys.argv[1]
file2 = sys.argv[2]

f1 = open(file1, "r")
f2 = open(file2, "r")
output = open("OUTPUT", "a")

lines1 = f1.readlines()
lines2 = f2.readlines()

if len(lines1) != len(lines2):
    print("Error: Alignment length is not the same")
    sys.exit()

for i in range(len(lines1)):
    line1 = lines1[i]
    line2 = lines2[i]
    if ">" in line1:
        output.write(line1)
        continue
    if "X" in line1 or "X" in line2:
        output.write(getIntersection(line1, line2))
        continue
    output.write(line1)



f1.close()
f2.close()
output.close()
