"""Find candidate bacterial genes (ORFs) with one regular expression, then translate them.

Run from the repository root:  python3 find_candidate_genes.py
"""
import re

MIN_CODONS = 100   # shorter ORFs are mostly chance

# ATG, then at least MIN_CODONS codons that are not stops, then a stop codon.
# The (?= ) lookahead lets matches overlap, so no ORF hides inside another.
ORF = re.compile(r"(?=(ATG(?:(?!TAA|TAG|TGA)[ACGT]{3}){%d,}(?:TAA|TAG|TGA)))" % MIN_CODONS)

# Genetic code: the 64 codons in TCAG order, paired with their amino acids in the same order
B = "TCAG"
CODE = dict(zip([a + b + c for a in B for b in B for c in B],
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))

# DNA from the GenBank file: the lines after ORIGIN, keeping only the bases
text = open("data/J01636_lac_operon.gb").read()
seq = re.sub("[^acgt]", "", text.split("\nORIGIN")[1].split("\n", 1)[1]).upper()

# The other strand: complement every base, then read it backwards
other = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

for strand, dna in (("+", seq), ("-", other)):
    longest = {}                                   # stop position -> longest ORF ending there
    for m in ORF.finditer(dna):
        longest.setdefault(m.end(1), m.group(1))   # the first ORF found for a stop is the longest

    for end, orf in longest.items():
        protein = "".join(CODE[orf[i:i + 3]] for i in range(0, len(orf) - 3, 3))   # skip the stop
        # GenBank-style positions on the + strand
        first, last = (end - len(orf) + 1, end) if strand == "+" else (len(seq) - end + 1, len(seq) - end + len(orf))
        print(f"{first}..{last} ({strand})  {len(protein)} aa  {protein[:40]}...")
