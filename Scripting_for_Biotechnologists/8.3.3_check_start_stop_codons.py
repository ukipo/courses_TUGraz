#!/usr/local/bin/python3

#6.3.3
from Bio import SeqIO
import sys

inargs = sys.argv
for i in inargs[1:]: #skip first since sys.argv[0] is the script name itself
    bad = 0 #for counting the number of sequences not starting with M
    for record in SeqIO.parse(i, "fasta"): #for each protein sequence
        if not record.seq.startswith("M"): #not starting with M
            bad = bad + 1 #add one for eac instance
            print(
              ("id = %s, length %i, with %i features")
              % (record.id, len(record.seq), len(record.features))   #print sequence description
              )
            print(record.id + " startswith " + record.seq[0])
    print(str(bad) + " sequences in " + i + " do not start with M")

