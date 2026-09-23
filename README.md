# bio116
Bioinformatics

## Regex ORF Lab

`orf-regex-lab.html` is an interactive workshop page (open it in a browser). Students build Python
regular expressions step by step to find bacterial open reading frames in the *E. coli* lac operon
(GenBank J01636.1, `data/J01636_lac_operon.gb`), translate them, and check every call against the
GenBank CDS annotation. Later steps search upstream for ribosome binding sites and sigma-70
promoters and measure how often those regexes match by chance. Students can also load their own
GenBank or FASTA file.

`orf_finder.py` is the reference solution: run `python3 orf_finder.py` from the repository root.
