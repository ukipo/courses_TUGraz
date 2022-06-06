#!/usr/local/bin/python3

from Bio import SeqIO
import sys

inargs = sys.argv
for i in inargs[1:]: #skip first since sys.argv[0] is the script name itself
    filename = i
    count = 0
    for record in SeqIO.parse(filename, "fasta"):
        count = count + 1
    print(str(count) + " sequences in " + filename)
