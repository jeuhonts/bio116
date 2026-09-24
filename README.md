# bio116
Bioinformatics

## Five-Step ORF Finder

`orf-lab-5steps.html` is a shorter, beginner version of the lab built around
`find_orfs_beginner.py`: (1) read the FASTA file, (2) make the complementary strand, (3) build the
ORF regex, (4) search both strands and keep the longest ORFs in lists, (5) translate with an
explicit codon dictionary, then check the results against GenBank.

## Regex ORF Lab

`orf-regex-lab.html` is an interactive workshop page (open it in a browser). Students build Python
regular expressions step by step to find bacterial open reading frames in the *E. coli* lac operon
(GenBank J01636.1: `data/J01636_lac_operon.fasta` for the sequence students start from, `data/J01636_lac_operon.gb` for the annotation), translate them, and check every call against the
GenBank CDS annotation. Later steps search upstream for ribosome binding sites and sigma-70
promoters and measure how often those regexes match by chance. Students can also load their own
GenBank or FASTA file.

Python, run from the repository root:

- `find_orfs_beginner.py` does the same in five numbered steps with only loops, lists and an
  explicit codon dictionary: read FASTA, complement, regex, search both strands, translate.
- `find_candidate_genes.py` is the minimal goal students build, starting from the FASTA file: ORFs on both strands (with
  overlaps), a length cutoff, and translation, with a written-out codon table and about 20 commented lines of code, no functions.
- `find_candidate_genes_biopython.py` gives the same result with Biopython (`pip install biopython`):
  translate the six frames, split at stops, take each piece from its first M.
- `orf_finder.py` is the full reference solution: alternative start codons, start choice by
  ribosome binding site, and comparison with the GenBank annotation.

Step 9 (optional) scores ORFs by codon usage. Its genome-wide *E. coli* K-12 table comes from the
Kazusa Codon Usage Database via the `python_codon_tables` package (CC0).

`orf-map.html` maps the forward-strand ORFs that three versions of the regex find in the lac operon (the lookahead
version from the lab, and two common online versions without a lookahead) against the GenBank genes.
