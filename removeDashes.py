import sys
import os


USAGE = "python removeDashes.py [alignment file]"
DESCRIPTION = "Removes columns that only has '-' characters and outputs to file 'removed_dash_output'"


if not len(sys.argv) == 2:
    print();
    print("\tError: Incorrect number of parameters\n");
    print("\tUSAGE: " + USAGE );
    print()
    print("\tDESCRIPTION: " + DESCRIPTION)
    print()
    sys.exit();



fasta_file = sys.argv[1];

f = open( fasta_file, "r" );
lines = f.readlines()

all_dash_columns = []
print("Number of sequences: " +  str(len(lines)/2)) 
print("Chars in sequence: " + str(len(lines[1])))


for i in range(len(lines[1])):
    all_dash = True
    for j in range(1, len(lines), 2):
        if lines[j][i] != '-':
            all_dash = False
            break
    if all_dash:
        all_dash_columns.append(i)


f.close()


f = open("removed_dash_output", "a");
for line in lines:
    if ( len(line) < len(lines[1]) ):
        f.write(line);
        continue
    string = line;
    for column in reversed(all_dash_columns):
        string = string[:column] + string[column+1:];
    f.write(string)

f.close();
