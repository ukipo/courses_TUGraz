#!/usr/local/bin/python3

#6.3.1
from Bio import SeqIO
import sys

inargs = sys.argv
for i in inargs[1:]: #skip first since sys.argv[0] is the script name itself
    for record in SeqIO.parse(i, "fasta"): #for each protein sequence
        if not record.seq.startswith("M"): #not starting with M
            print(record.id + " startswith " + record.seq[0])

