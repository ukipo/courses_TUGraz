#!/usr/local/bin/python3

# This script calculates the GC percentage of a sequence.

seq = "GTATGGAGCTCTGTGAGTGTGTGCTTCTTTTATTATGGATCACGTTCAGCTTCTTTTTATTTGCCTCCTA"

# Get the length of the sequence, assign to a variable
seq_length = len(seq)

# Count the number of G's, assign to a variable
g_count = seq.count("G")

# TO DO: Count the number of C's, assign to a variable
c_count = seq.count("C")

# TO DO: Calculate the GC percentage
gc_perc = (g_count + c_count)/seq_length

# Print informative messages, using a separate print statement for each item
# Print the sequence
print("The Sequence is " + seq)

# TO DO: Print the sequence length
print("The Sequence length is " + str(seq_length))

# Print the G count
print("The number of G bases is " + str(g_count))

# TO DO: Print the C count
print("The number of C bases is " + str(c_count))

# TO DO: Print the GC percentage
print("The GC content of the Sequence is " + str(gc_perc))
