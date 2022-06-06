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


# List of dna sequences to test 
dna_list = [
  "ATGTTCGGTGAGCTGTAA", 
  "ATGTTCGGTGAGCTGTAAG",
  "ATGTTCGGTGAGCTGTAAGC",
  "ATGTNCGGTGAGCTGTAAGN",
  "ATGTNCGGTGAGCTGTAAGNT",
  "ATGTNCGGTGAGCTGTAAGTT"
]


# Start a for loop to iterate over the DNA sequences in the list
for dna in dna_list:

  # Print the DNA sequence
  print("DNA: " + dna)

  # TO DO: Print the translated protein sequence
  print("Protein: " + dna_to_protein(dna))

