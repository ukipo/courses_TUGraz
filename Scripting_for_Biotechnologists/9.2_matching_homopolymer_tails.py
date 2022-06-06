#!/usr/local/bin/python3

import re

#test with sequences:
seq1 = "ACTGGGTCTCTAAAGGGTCGTAAAAAACTCTCAAAAAAAAA"
seq2 = "GGTCTCTAAAGGGTCGCCCCCCCCTCTCGGGGGGGGG"
seq3 = "tgggtctctaaagggtcgtaaaaaactccccccccttttt"
seq4 = "CTGGGTCTCTAAAGGGTCGTAAAAAACTCTCCCCCC"
seq5 = "tgggtctctaaagggtcgtaaaaaactcccccccctttaaaaaaa"

seq = seq5

#case-insensitive regular expresson
#poly-A sequences whose length is at least 6

x = re.findall("[Aa]{6}[Aa]*", seq)
print(x)

if (x):
    print("Matches found")
else:
    print("No match")

#poly-A tails of 6 or more

x = re.findall("[Aa]{6}[Aa]*$", seq)
print(x)

if (x):
    print("Matches found")
else:
    print("No match")

#all homopolymer tails of minimum 6 bases, such
#the match happens regardless whether it is an
#A, C, G or T that is repeated

x = re.findall("([Aa]{6}[Aa]*$)|([Gg]{6}[Gg]*$)|([Tt]{6}[Tt]*$)|([Cc]{6}[Cc]*$)", seq)
print(x)

if (x):
    print("Matches found")
else:
    print("No match")