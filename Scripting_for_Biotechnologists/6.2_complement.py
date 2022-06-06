#!/usr/local/bin/python3

# This script complements a sequence.

seq = "GTATGGAGCTCTGTGAGTGTGTGCTTCTTTTATTATGGATCACGTTCAGCTTCTTTTTATTTGCCTCCTA"

# Replace A's with t's, assign to a new variable
comp_seq = seq.replace("A", "t")

# TO DO: Replace T's with a's in the complement sequence
comp_seq = comp_seq.replace("T", "a")

# TO DO: Replace C's with g's in the complement sequence
comp_seq = comp_seq.replace("C", "g")

# TO DO: Replace G's with c's in the complement sequence
comp_seq = comp_seq.replace("G", "c")

# TO DO: Convert all lowercase characters in complement sequence to uppercase
comp_seq.upper()

# Print the orignial sequence, with an informative messages
print("Original sequence: " + seq)

# TO DO: Print the complement sequence, with an informative message
print("Complement: " + comp_seq.upper())
