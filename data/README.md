# Data

## `dna_sequence.fasta`

A single ~10.3 kb DNA sequence used as the **offline fallback** for the
exercise notebook, so the whole lesson runs even when NCBI is unreachable
(exam rooms, no wifi, locked-down networks, CI sandboxes).

**Provenance.** It was reconstructed from the repository's original
`triplets.txt` file. That file contains the *same* DNA read three times, once in
each forward reading frame (each copy shifted by one base and printed as
whitespace-separated codons). This FASTA is the frame-0 read with all
whitespace and the `start N` markers removed:

```
length   : 10,309 bp
GC content: 65.9%  (GC-rich, typical of a human gene-dense region)
```

It is deliberately presented to students as an **uncharacterized sequence** —
the point of the exercise is to characterize it with code (GC landscape,
six-frame ORF search, codon usage, translation). Its longest open reading frame
is 330 aa on the + strand and translates cleanly (starts with `M`, ends with a
single stop, no internal stops) — a satisfying, verifiable result.

To regenerate it from scratch:

```bash
python - <<'PY'
import re, textwrap
text = open("../triplets.txt").read()
block1 = re.search(r'start\s+1(.*?)start\s+2', text, re.S).group(1)
seq = "".join(re.findall(r'[ACGT]', block1))
open("dna_sequence.fasta","w").write(
    ">local_mystery_seq reconstructed from triplets.txt | offline fallback dataset\n"
    + "\n".join(textwrap.wrap(seq, 70)) + "\n")
PY
```

## Live data

The notebook's *primary* data path fetches real records straight from NCBI with
Biopython (`Bio.Entrez` / `Bio.SeqIO`). When the network is available, students
work with a freshly downloaded gene (default: human β-globin, `HBB`); the file
here is only the safety net.
