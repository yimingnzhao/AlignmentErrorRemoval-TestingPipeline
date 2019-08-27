import treeswift
import sys
import random

"""
Checks if an input string is a float

Args:
    num (str): the string to check

Return:
    bool: whether the string is a float or not
"""
def isFloat( num ):
    try:
        float(num)
        return True
    except ValueError:
        return False


USAGE = "python chooseSequencesByDistance.py [newick tree file] [alignment file] [minimum diameter] [maximum diameter]"
DESCRIPTION = "Chooses a random tree that fits in the diameter range and outputs to file 'chosen_sequences.fasta'"
if not len(sys.argv) == 5:
    print()
    print("\tError: Incorrect number of parameters\n")
    print("\tUSAGE: " + USAGE)
    print()
    print("\tDESCRIPTION: " + DESCRIPTION)
    print()
    sys.exit()

if not isFloat(sys.argv[3]):
    print()
    print("\tError: [minimum diameter] is not a float")
    print("\tUSAGE: " + USAGE)
    sys.exit()

if not isFloat(sys.argv[4]):
    print()
    print("\tError: [maximum diameter] is not a float")
    print("\tUSAGE: " + USAGE)
    sys.exit()

newick_file = sys.argv[1]
align_file = sys.argv[2]
min_diameter = float(sys.argv[3])
max_diameter = float(sys.argv[4])

tree = treeswift.read_tree_newick(newick_file)
trees = []
max_alignments = 0
for node in tree.traverse_postorder(internal=True, leaves=False):
    subtree = tree.extract_subtree(node)
    diameter = subtree.diameter()
    if diameter <= max_diameter and diameter > min_diameter:
        if subtree.num_nodes(internal=False) > max_alignments:
            max_alignments = subtree.num_nodes(internal=False)
            trees = []
            trees.append(subtree)
        if subtree.num_nodes(internal=False) == max_alignments:
            trees.append(subtree)

if len(trees) == 0:
    print("Error: No clades found with the given diameter parameters");
    sys.exit()

chosen_tree = trees[random.randint(0, len(trees)-1)]

chosen_tree_labels = []
for node in chosen_tree.traverse_postorder(internal=False):
    chosen_tree_labels.append(">" + node.get_label())

output_f = open("chosen_sequences.fasta", "a")

count = 0;

f = open(align_file, "r")
line = f.readline().strip('\n')
while line:
    if line in chosen_tree_labels:
        output_f.write(line + '\n')
        output_f.write(f.readline())
    else:
        f.readline()
    line = f.readline().strip('\n')
output_f.close()
