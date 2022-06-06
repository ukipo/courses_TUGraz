#!/usr/local/bin/python3

from Bio import SeqIO
import sys


with open(sys.argv[1], 'r') as filename:
    count = 0
    for record in SeqIO.parse(filename, "fasta"):
        count = count + 1
    print(str(count) + " sequences in " + sys.argv[1])
