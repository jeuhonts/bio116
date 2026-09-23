"""Find candidate bacterial genes (ORFs) with Biopython instead of a regular expression.

Same result as find_candidate_genes.py. Idea: translate each of the six reading frames,
split the protein at every stop (*), and take each piece from its first M to the stop.

Needs Biopython:  pip install biopython
Run from the repository root:  python3 find_candidate_genes_biopython.py
"""
from Bio import SeqIO

MIN_CODONS = 100   # codons between start and stop; shorter ORFs are mostly chance

record = SeqIO.read("data/J01636_lac_operon.fasta", "fasta")   # one sequence in the file
seq = record.seq.upper()

orfs = []
for strand, dna in (("+", seq), ("-", seq.reverse_complement())):
    for frame in range(3):                                  # the three reading frames
        codons = dna[frame:]
        codons = codons[: len(codons) // 3 * 3]             # drop leftover bases at the end
        protein = str(codons.translate(table=11))           # bacterial code; * marks each stop

        codon_index = 0                                     # where the current piece starts
        pieces = protein.split("*")                         # stretches between stop codons
        for piece in pieces[:-1]:                           # the last piece has no stop after it
            m = piece.find("M")                             # first ATG = longest ORF for this stop
            if m != -1 and len(piece) - m - 1 >= MIN_CODONS:
                start = frame + 3 * (codon_index + m)       # 0-based, on this strand
                end = frame + 3 * (codon_index + len(piece) + 1)   # includes the stop codon
                if strand == "-":                           # convert to + strand positions
                    start, end = len(seq) - end, len(seq) - start
                orfs.append((start + 1, end, strand, piece[m:]))
            codon_index += len(piece) + 1                   # skip past this piece and its stop

for first, last, strand, protein in sorted(orfs):
    print(f"{first}..{last} ({strand})  {len(protein)} aa  {protein[:40]}...")
