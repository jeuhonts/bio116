"""Find candidate bacterial genes (ORFs) with one regular expression, then translate them.

Run from the repository root:  python3 find_candidate_genes.py
"""
import re

MIN_CODONS = 100   # codons between start and stop; shorter ORFs are mostly chance

# Build the pattern from named pieces:  START  CODON{100,}  STOP
START = "ATG"                                  # try "(?:ATG|GTG|TTG)" for alternative starts
STOP = "(?:TAA|TAG|TGA)"
CODON = "(?:(?!" + STOP + ")...)"              # any 3 bases that are not a stop codon
# (?= ) is a lookahead: matches may overlap, so no ORF hides inside another; group 1 is the ORF
ORF = re.compile("(?=(" + START + CODON + "{%d,}" % MIN_CODONS + STOP + "))")

# Genetic code: codon -> amino acid (one-letter code, * = stop).
# Laid out like the textbook table: each row varies the second base (T, C, A, G).
CODE = {
    # first base T
    "TTT": "F",  "TCT": "S",  "TAT": "Y",  "TGT": "C",
    "TTC": "F",  "TCC": "S",  "TAC": "Y",  "TGC": "C",
    "TTA": "L",  "TCA": "S",  "TAA": "*",  "TGA": "*",
    "TTG": "L",  "TCG": "S",  "TAG": "*",  "TGG": "W",
    # first base C
    "CTT": "L",  "CCT": "P",  "CAT": "H",  "CGT": "R",
    "CTC": "L",  "CCC": "P",  "CAC": "H",  "CGC": "R",
    "CTA": "L",  "CCA": "P",  "CAA": "Q",  "CGA": "R",
    "CTG": "L",  "CCG": "P",  "CAG": "Q",  "CGG": "R",
    # first base A
    "ATT": "I",  "ACT": "T",  "AAT": "N",  "AGT": "S",
    "ATC": "I",  "ACC": "T",  "AAC": "N",  "AGC": "S",
    "ATA": "I",  "ACA": "T",  "AAA": "K",  "AGA": "R",
    "ATG": "M",  "ACG": "T",  "AAG": "K",  "AGG": "R",
    # first base G
    "GTT": "V",  "GCT": "A",  "GAT": "D",  "GGT": "G",
    "GTC": "V",  "GCC": "A",  "GAC": "D",  "GGC": "G",
    "GTA": "V",  "GCA": "A",  "GAA": "E",  "GGA": "G",
    "GTG": "V",  "GCG": "A",  "GAG": "E",  "GGG": "G",
}

# DNA from the FASTA file: skip the ">" header line, join the rest into one string
lines = open("data/J01636_lac_operon.fasta").read().splitlines()
seq = "".join(line.strip() for line in lines if not line.startswith(">")).upper()

# The other strand: complement every base, then read it backwards
other = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

for strand, dna in (("+", seq), ("-", other)):
    longest = {}                                   # stop position -> longest ORF ending there
    for m in ORF.finditer(dna):
        longest.setdefault(m.end(1), m.group(1))   # the first ORF found for a stop is the longest

    for end, orf in longest.items():
        protein = "".join(CODE.get(orf[i:i + 3], "X") for i in range(0, len(orf) - 3, 3))   # skip the stop; X = unknown (e.g. N)
        # positions on the + strand, counting from 1 (the same numbering GenBank uses)
        first, last = (end - len(orf) + 1, end) if strand == "+" else (len(seq) - end + 1, len(seq) - end + len(orf))
        print(f"{first}..{last} ({strand})  {len(protein)} aa  {protein[:40]}...")
