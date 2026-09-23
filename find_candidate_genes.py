"""Find candidate bacterial genes with regular expressions: the minimal goal of the Regex ORF Lab.

A candidate gene (open reading frame, ORF) is:
  a start codon, then at least MIN_CODONS codons with no stop, then a stop codon,
on either strand of the DNA. Each one is translated into protein.

Run from the repository root:  python3 find_candidate_genes.py
"""
import re

GENBANK_FILE = "data/J01636_lac_operon.gb"   # E. coli lac operon from GenBank
MIN_CODONS = 100                             # shorter ORFs are mostly chance

# The whole search in one regular expression:
#   (?= ... )                  lookahead, so ORFs that overlap are all found
#   ATG                        start codon
#   (?:(?!TAA|TAG|TGA)[ACGT]{3}){100,}   100+ codons, none of them a stop
#   (?:TAA|TAG|TGA)            the first in-frame stop
ORF_RE = re.compile(r"(?=(ATG(?:(?!TAA|TAG|TGA)[ACGT]{3}){%d,}(?:TAA|TAG|TGA)))" % MIN_CODONS)

# Genetic code: codons listed in TCAG order, first base changing slowest
BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {a + b + c: AMINO[16*i + 4*j + k]
               for i, a in enumerate(BASES)
               for j, b in enumerate(BASES)
               for k, c in enumerate(BASES)}


def read_sequence(path):
    """Return the DNA after the ORIGIN line of a GenBank file, as A/C/G/T only."""
    text = open(path).read()
    origin = text.split("\nORIGIN", 1)[1].split("\n", 1)[1]   # skip the ORIGIN line itself
    return re.sub(r"[^ACGT]", "", origin.split("//")[0].upper())


def reverse_complement(dna):
    """The other strand, read 5' to 3'."""
    return dna.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def translate(orf):
    """Codon by codon; drop the final stop (*)."""
    codons = re.findall(r"[ACGT]{3}", orf)
    return "".join(CODON_TABLE[c] for c in codons).rstrip("*")


def find_orfs(seq):
    """Yield (start, end, strand, orf) with 1-based + strand coordinates, like GenBank."""
    length = len(seq)
    for strand, dna in (("+", seq), ("-", reverse_complement(seq))):
        seen_stops = set()
        for m in ORF_RE.finditer(dna):
            orf = m.group(1)                 # the match itself is empty; the ORF is group 1
            start, end = m.start(1), m.start(1) + len(orf)
            if end in seen_stops:            # an inner ATG of an ORF we already have
                continue
            seen_stops.add(end)              # keep only the longest ORF for each stop
            if strand == "-":                # convert back to + strand coordinates
                start, end = length - end, length - start
            yield start + 1, end, strand, orf


if __name__ == "__main__":
    seq = read_sequence(GENBANK_FILE)
    print(f"{len(seq)} bp read from {GENBANK_FILE}")
    for start, end, strand, orf in sorted(find_orfs(seq)):
        protein = translate(orf)
        print(f"{start:>6}..{end:<6} {strand}  {len(protein):>4} aa  {protein[:40]}...")
