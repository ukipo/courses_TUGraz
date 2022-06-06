#!/usr/local/bin/python3

# This script defines the genetic code, defines a function to translate 
# a codon to an amino acid code based on the genetic code, and tests
# translations of several DNA sequences using the function.

# Initialize a  dictionary to set the genetic code. There are 4x4x4=64 
# key:value pairs, # where each pair maps (translates) a 3-letter codon 
# to a 1-letter amino acid code.
genetic_code = {
  'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
  'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
  'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
  'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
  'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
  'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
  'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
  'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
  'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
  'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
  'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
  'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
  'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
  'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
  'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
  'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
}


# Define a function to translate a dna sequence to a protein sequence
def dna_to_protein(dna):

  # Get the last possible 1-based codon position 
  last_codon_start = len(dna) - 2

  # Intitialize the protein sequence as blank
  protein = ""

  # Start a for loop to iterate over possible (0-based) codon start positions
  for start in range(0, last_codon_start, 3):

    # TO DO: Extract the codon from the DNA, assign it to a variable
    end = start + 3
    codon = dna[start:end]

    # TO DO: Look up the amino acid translation of this codon. If a matching
    # translation is not found, set the amino acid code to X to mean unknown.
    AA = genetic_code.get(codon)
    if AA == None:
      AA = "X"

    # TO DO: Append the amino acid to the end of the protein sequence
    protein = protein + AA

  # Return the protein sequence built so far
  return protein


#import gzip
import gzip

#define function to open file with either gzip or the normal way
def opening_seq(selected_file):
  if selected_file.endswith(".gz"):
    open_file =  gzip.open(selected_file, "rt")
    print("File opened with gzip")
  else:
    open_file =  open(selected_file)
    print("File opened the good ol' way")
  return open_file

#open a file - first selected_file is the gzipped file, the second is a fasta file from previous labs to test if it works
#selected_file = ("/export/home/up0289/courses/scripting_for_biotechnologists/lab08/adh1.fasta.gz")
#selected_file = ("/export/home/up0289/courses/scripting_for_biotechnologists/lab07/exons.fasta")
sequence = opening_seq("/export/home/up0289/courses/scripting_for_biotechnologists/lab08/adh1.fasta.gz")

#open a new file in which it will write the translations
translated_sequence = open("/export/home/up0289/courses/scripting_for_biotechnologists/lab08/translated_sequences.fasta", "w")

#translate the sequences using dna_to_protein. do not translate the gene info. write gene info and translation into the new file
for line in sequence:
  line = line.rstrip("\n")
  if line.startswith('>'): #gene info should not be translated but just written down again
    title = line.rstrip()
    print(title)
    translated_sequence.write("\n"+ title + "\n")
  else:
    print(dna_to_protein(line))
    translated_sequence.write(dna_to_protein(line))

sequence.close
translated_sequence.close
