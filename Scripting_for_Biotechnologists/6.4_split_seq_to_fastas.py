#!/usr/local/bin/python3

# This script reads a FASTA file of one sequence containing 3 exons and 2 introns,
# and generates two FASTA files. One FASTA file will have the 3 exons, and the other
# FASTA file will have the 2 introns.

# Open the input sequence file
seq_file = open(r"/export/home/up0289/lab06/AMY1A-202.fa", 'r')

# Read the first line (sequence header) and ignore
header = seq_file.readline()

# Initialize the sequence to be blank
seq = ""

# TO DO: Start a for loop to read sequence lines one by one
for line in seq_file:
  
   # TO DO: Append the line (with newline removed) to the growing sequence
   line = line.rstrip("\n")
   seq += line

# Close the input file
seq_file.close()
print("Sequence length: ", len(seq))

# TO DO: Extract exons and introns  
# exon 1: bases 1 to 42, intron 1: bases 43 to 267, exon 2: bases 268 to 481,
# intron 2: bases 482 to 826, exon 3: bases 827 to 973
exon_1 = seq[0:42]
intron_1 = seq[42:267]
exon_2 = seq[267:481]
intron_2 = seq[481:826]
exon_3 = seq[826:973]

# Open a file to write exon sequences
exons_file = open("AMY1A-202_exons.fa", "w")

# Open another file to write intron sequences
introns_file = open("AMY1A-202_introns.fa", "w")

# TO DO: Write each exon to the exon sequence file
# Write an informative header (e.g. >AMY1A-202_exon1), newline
# exon sequence, then newline
exons_file.write(">AMY1A-202_exon1\n")
exons_file.write(exon_1 + "\n")
exons_file.write(">AMY1A-202_exon2\n")
exons_file.write(exon_2 + "\n")
exons_file.write(">AMY1A-202_exon3\n")
exons_file.write(exon_3 + "\n")

print("Exon 2 length:", len(exon_1))
print("Exon 2 length:", len(exon_2))
print("Exon 3 length:", len(exon_3))

# Close the exon sequence output file
exons_file.close()

# Write each intron to the intron sequence file
# Write an informative header (e.g. >AMY1A-202_intron1), newline
# intron sequence, then newline
introns_file.write(">AMY1A-202_intron1\n")
introns_file.write(intron_1 + "\n")
introns_file.write(">AMY1A-202_intron2\n")
introns_file.write(intron_2 + "\n")

print("Intron 1 length:", len(intron_1))
print("Intron 3 length:", len(intron_2))

# Close the intron sequence output file
introns_file.close()
