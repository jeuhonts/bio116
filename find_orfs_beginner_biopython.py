# Find candidate genes (ORFs) in the lac operon and translate them, using Biopython.
# Same five steps and same output as find_orfs_beginner.py, for comparison.
#
# Needs Biopython:  pip install biopython
# Run from the repository root:  python3 find_orfs_beginner_biopython.py

import re
from Bio import SeqIO


# ------------------------------------------------------------------
# Step 1: Read the FASTA file and extract the sequence
# ------------------------------------------------------------------
# SeqIO.read opens the file, skips the ">" header line, joins the sequence
# lines and returns one record. It stops with an error if the file has
# no records or more than one.

record = SeqIO.read("data/J01636_lac_operon.fasta", "fasta")

print("Record ID:", record.id)
sequence = record.seq.upper()        # a Biopython Seq object, not a plain string

print("Sequence length:", len(sequence), "bases")


# ------------------------------------------------------------------
# Step 2: Derive the complementary strand
# ------------------------------------------------------------------
# Biopython knows the base pairs, so this is one method call.
# reverse_complement() complements AND reverses, giving the other strand 5' to 3'.

reverse_complement = sequence.reverse_complement()


# ------------------------------------------------------------------
# Step 3: A regular expression that finds ORFs
# ------------------------------------------------------------------
# Biopython has no ORF-finding function, so we use the same regex.
#   ATG                        start codon
#   (?:(?!TAA|TAG|TGA)...)     one codon (any 3 bases) that is NOT a stop codon
#   {100,}                     at least 100 of those codons
#   (?:TAA|TAG|TGA)            stop codon: the first one in the same frame
#   (?=( ... ))                lookahead: lets ORFs overlap, so none are missed;
#                              the ORF itself is saved in group 1

orf_pattern = r"(?=(ATG(?:(?!TAA|TAG|TGA)...){100,}(?:TAA|TAG|TGA)))"


# ------------------------------------------------------------------
# Step 4: Search both strands and keep the longest ORF for each stop codon
# ------------------------------------------------------------------
# re works on plain strings, so str() converts each Seq before searching.
# The first ORF we see for each stop is the longest (matches come left to right).

orfs = []        # the ORF sequences (as Seq objects, so step 5 can translate them)
strands = []     # "+" or "-" for each ORF
starts = []      # start position of each ORF on its strand (counting from 1)

for strand in ["+", "-"]:
    if strand == "+":
        dna = sequence
    else:
        dna = reverse_complement

    stops_seen = []                          # stop positions already used on this strand
    for match in re.finditer(orf_pattern, str(dna)):
        stop_position = match.end(1)
        if stop_position in stops_seen:      # a shorter ORF inside one we already have
            continue
        stops_seen.append(stop_position)

        orfs.append(dna[match.start(1):stop_position])   # slicing a Seq gives a Seq
        strands.append(strand)
        starts.append(match.start(1) + 1)

print("ORFs found:", len(orfs))


# ------------------------------------------------------------------
# Step 5: Translate each ORF with the genetic code
# ------------------------------------------------------------------
# translate() reads the ORF codon by codon using Biopython's built-in
# genetic code tables. Table 11 is the bacterial code (same codon meanings
# as the standard code). The stop codon is written as "*".

for i in range(len(orfs)):
    orf = orfs[i]
    protein = orf.translate(table=11)

    print()
    print("ORF", i + 1, "on strand", strands[i], "starting at", starts[i])
    print("  length:", len(orf), "bases,", len(protein) - 1, "amino acids")
    print("  protein:", protein)
