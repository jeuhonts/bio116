# bio116
Bioinformatics

## Genome Browser Lab

Open [`genome-browser/index.html`](genome-browser/index.html) in a web browser: an interactive genome browser with guided missions on two real human regions (switch with the tabs above the browser):

- **HRAS**: RefSeqGene NG_007666.1, saved as [`NG_007666_HRAS.gb`](NG_007666_HRAS.gb) (same sequence as `triplets.txt`) (12 missions: navigation, a minus-strand neighbour gene (LRRC56), base-level zoom, splice sites, oncogenic and synonymous variants, strand orientation, CpG island, HRAS1 minisatellite, BED coordinates).
- **β-globin cluster**: GenBank U01317, saved as [`U01317_beta_globin_cluster.gb`](U01317_beta_globin_cluster.gb) (10 missions: gene order and development, the ψβ pseudogene, the sickle-cell mutation, splice-site and cryptic-splice thalassemia alleles, Gγ vs Aγ, HPFH promoter variants, strand orientation, Alu repeats).

### Graded exercise

The lab page ends with a graded exercise. Each student enters their student ID, which generates their own question set (a globin coding position and a hypothetical variant, an HRAS splice site, a codon in genome orientation, a BED line, a UCSC hg38 gene, and a written question). Pasting into the answer boxes is blocked, no answer key is shown, and **Submit** saves `BIO116-genome-browser-<ID>.json` for upload to the course site.

Grade with [`genome-browser/grader.html`](genome-browser/grader.html) (instructor only; do not share it with students, as it shows the answer keys): drop in the submission files to auto-score, flag edited files and paste attempts, add written scores, and export a CSV. Question generation and scoring live in `genome-browser/exercise.js`; the sequence data is in `genome-browser/loci.js`.
