#!/usr/local/bin/python3

#6.3.4
from Bio import SeqIO
import sys

inargs = sys.argv
for i in inargs[1:]: #skip first since sys.argv[0] is the script name itself
    nM = 0 #not starting with M
    nst = 0 #not ending with stop *
    nboth = 0 #not starting with M and ending with stop *
    for record in SeqIO.parse(i, "fasta"): #for each protein sequence
        if not record.seq.startswith("M"): #not starting with M
            nM = nM + 1
            if not record.seq.endswith("*"): #nested here, from those not starting with M, aslo not ending with *
              nboth = nboth + 1
        if not record.seq.endswith("*"): #not ending with *
          nst = nst + 1
    print(str(nM) + " sequences in " + i + " do not start with M")
    print(str(nst) + " sequences in " + i + " do not end with *")
    print(str(nboth) + " sequences in " + i + " do not start with M and end with *")

