#!/usr/local/bin/python3

# This script performs 3 separate tasks for a sequence consisting of 2 exons
# and 1 intron.

seq = "ATCGATCGATCGATCGACTGACTAGTCATAGCTATGCATGTAGCTACTCGATCGATCGATCGATCGATCGATCGATCGATCGATCATGCTATCATCGATCGATATCGATGCATCGACTACTAT"

###### Subtask 1: Extract the exonic region

# Extract the first exon (bases 1 to 63)
exon1 = seq[0:63]

# TO DO: Extract the second exon (bases 91 to end)
exon2 = seq[90:]

# Print a message to show the first exon
print("Exon 1: " + exon1)

# TO DO: Print a message to show the second exon
print("Exon 2: " + exon2)

# TO DO: Print a message to show the exonic sequence (exons put together)
print("Exonic sequence: " + exon1 + exon2)

###### Subtask 2 : Calculate the exonic region percentage

# Get the total sequence length
seq_length = len(seq)

# TO DO: Get the exonic region length
exon_length = len(exon1) + len(exon2)

# TO DO: Calculate the percentage 
exon_perc = exon_length/seq_length

# TO DO: Print a message to tell the exonic region percentage
print("Exonic region percentage: " + str(exon_perc))

###### Subtask 3 : Print the original sequence with introns masked in lowercase

# TO DO: Extract the intron (bases 64 to 90)
intron = seq[63:90]

# TO DO: Print a message to show the intron-masked sequence:
# first exon, intron in lower case, and second exon
print("Intron-masked sequence: " + exon1 + intron.lower() + exon2)
