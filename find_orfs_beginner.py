# Find candidate genes (ORFs) in the lac operon and translate them.
# Run from the repository root:  python3 find_orfs_beginner.py

import re


# ------------------------------------------------------------------
# Step 1: Read the FASTA file and extract the sequence
# ------------------------------------------------------------------
# A FASTA file has one header line that starts with ">",
# followed by lines of sequence. We skip the header and join the rest.

# "with" opens the file and closes it again automatically at the end of the block.
# readlines() gives a list of strings, one per line, each ending in "\n".
with open("data/J01636_lac_operon.fasta") as fasta_file:
    lines = fasta_file.readlines()

print("Lines in the file:", len(lines))
print("First line:", lines[0].strip())

sequence = ""
for line in lines:
    line = line.strip()              # remove the "\n" at the end
    if line.startswith(">"):         # header line: skip it
        continue
    sequence = sequence + line.upper()
# (For a whole genome, collect the lines in a list and use "".join(list) instead: it is much faster.)

print("Sequence length:", len(sequence), "bases")


# ------------------------------------------------------------------
# Step 2: Derive the complementary strand
# ------------------------------------------------------------------
# Each base pairs with its partner: A-T and C-G.
# The other strand runs in the opposite direction, so after
# complementing we reverse it to read it 5' to 3' like the first strand.

pairs = {"A": "T", "T": "A", "C": "G", "G": "C"}

complement = ""
for base in sequence:
    complement = complement + pairs[base]

reverse_complement = complement[::-1]    # [::-1] reverses a string


# ------------------------------------------------------------------
# Step 3: A regular expression that finds ORFs
# ------------------------------------------------------------------
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
# The lookahead also finds shorter ORFs that start at an ATG inside a
# longer one and end at the same stop codon. Matches come in order from
# left to right, so the first ORF we see for each stop is the longest.

orfs = []        # the ORF sequences
strands = []     # "+" or "-" for each ORF
starts = []      # start position of each ORF on its strand (counting from 1)

for strand in ["+", "-"]:
    if strand == "+":
        dna = sequence
    else:
        dna = reverse_complement

    stops_seen = []                          # stop positions already used on this strand
    for match in re.finditer(orf_pattern, dna):
        orf = match.group(1)
        stop_position = match.end(1)
        if stop_position in stops_seen:      # a shorter ORF inside one we already have
            continue
        stops_seen.append(stop_position)

        orfs.append(orf)
        strands.append(strand)
        starts.append(match.start(1) + 1)

print("ORFs found:", len(orfs))


# ------------------------------------------------------------------
# Step 5: Translate each ORF with the genetic code
# ------------------------------------------------------------------
# Genetic code: codon -> amino acid (one-letter code, * = stop).
# Laid out like the textbook table.

genetic_code = {
    "TTT": "F",  "TCT": "S",  "TAT": "Y",  "TGT": "C",
    "TTC": "F",  "TCC": "S",  "TAC": "Y",  "TGC": "C",
    "TTA": "L",  "TCA": "S",  "TAA": "*",  "TGA": "*",
    "TTG": "L",  "TCG": "S",  "TAG": "*",  "TGG": "W",

    "CTT": "L",  "CCT": "P",  "CAT": "H",  "CGT": "R",
    "CTC": "L",  "CCC": "P",  "CAC": "H",  "CGC": "R",
    "CTA": "L",  "CCA": "P",  "CAA": "Q",  "CGA": "R",
    "CTG": "L",  "CCG": "P",  "CAG": "Q",  "CGG": "R",

    "ATT": "I",  "ACT": "T",  "AAT": "N",  "AGT": "S",
    "ATC": "I",  "ACC": "T",  "AAC": "N",  "AGC": "S",
    "ATA": "I",  "ACA": "T",  "AAA": "K",  "AGA": "R",
    "ATG": "M",  "ACG": "T",  "AAG": "K",  "AGG": "R",

    "GTT": "V",  "GCT": "A",  "GAT": "D",  "GGT": "G",
    "GTC": "V",  "GCC": "A",  "GAC": "D",  "GGC": "G",
    "GTA": "V",  "GCA": "A",  "GAA": "E",  "GGA": "G",
    "GTG": "V",  "GCG": "A",  "GAG": "E",  "GGG": "G",
}

for i in range(len(orfs)):
    orf = orfs[i]

    protein = ""
    for position in range(0, len(orf), 3):   # step through the ORF one codon at a time
        codon = orf[position:position + 3]
        protein = protein + genetic_code[codon]

    print()
    print("ORF", i + 1, "on strand", strands[i], "starting at", starts[i])
    print("  length:", len(orf), "bases,", len(protein) - 1, "amino acids")
    print("  protein:", protein)
