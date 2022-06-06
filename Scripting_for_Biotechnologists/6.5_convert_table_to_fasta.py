#!/usr/local/bin/python3

# This script reads a tab-separated data file of exon information in a transcript,
# and creates a FASTA file of exon sequences, where information about each exon is
# included in the sequence header.

# There are 6 columns in the data file: exon ID, chromosome number, start position, 
# end position, length, and sequence.

# TO DO: Open the input table file
table_file = open("/export/home/up0289/lab06/BRCA2-202_exons.tsv", "r")
fasta_file = open("/export/home/up0289/lab06/BRCA2-202_fasta.tsv", 'w')
# TO DO: Start a for loop to read sequence lines one by one
for line in table_file:
  # Strip off the newline at the end of line
  line = line.rstrip()

  # Split the line using tab
  fields = line.split("\t")

  # Assign field 1 to ID
  exon_id = (fields[0])

  # TO DO: Assign the remaining fields to appropriate variables
  chromosome_no = str(fields[1])
  start = str(fields[2])
  end = str(fields[3])
  length = str(fields[4])
  sequence = str(fields[5])

  # TO DO: Print a very informative header for this exon. For example:
  # >ENSE00001930817_chr13:32379840_32379913_74bp
  # Remember to use str() function to convert a number to a string
  fasta_file.write('>'+ exon_id + "_chr" + chromosome_no + ":" + start + "_" + end + "_" + length + "bp" + "\n")

  # TO DO: Print the sequence
  fasta_file.write(sequence + "\n\n")

# TO DO: Close the input file
table_file.close()
fasta_file.close()
