# Tutorial: Protein Structure Prediction

**From a DNA sequence to a 3D model, and how to judge whether to trust it**

| | |
|---|---|
| **Audience** | Undergraduate bioinformatics (BIO116). Assumes you know the genetic code and the basics of protein structure. |
| **Time** | About 2–3 hours (Parts 1–4 offline, about 1 h; Part 5 needs a browser and about 30 min of compute) |
| **You need** | Python 3.8+ (standard library only), a web browser, and optionally [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) or the [Mol\* viewer](https://molstar.org/viewer/) |
| **Files** | `structure_tutorial.py` (companion script), `../../triplets.txt` (the DNA we start from) |

By the end you will be able to:

1. Find and translate the protein-coding sequence hidden in a stretch of genomic DNA.
2. Make and critique classical *sequence-only* predictions: hydropathy, secondary structure and disorder.
3. Run a modern deep-learning predictor (AlphaFold2 through ColabFold, ESMFold, or AlphaFold Server).
4. Read the confidence metrics (pLDDT, PAE, pTM) and compare a prediction against an experimental structure.
5. Explain what structure predictors **cannot** tell you.

---

## Part 0 — Background: why predicting structure is hard (and why it matters)

### Levels of protein structure

| Level | What it describes | Driving interactions |
|---|---|---|
| Primary | Amino-acid sequence | Covalent peptide bonds |
| Secondary | Local α-helices, β-strands, turns and loops | Backbone N–H···O=C hydrogen bonds |
| Tertiary | 3D fold of one chain | Hydrophobic burial, side-chain H-bonds, salt bridges, disulfides |
| Quaternary | Assembly of several chains | The same forces, acting across the interfaces between chains |

**Anfinsen's dogma (1961–73):** for many small globular proteins, the sequence alone determines the native
fold. So, in principle, structure is *computable* from sequence.

**Levinthal's paradox (1969):** a 100-residue chain with about 3 backbone conformations per residue has about
3¹⁰⁰ ≈ 10⁴⁸ possible conformations. Proteins can't search them all, but they still fold in milliseconds, and brute
force won't work for us either. We need *knowledge-based* shortcuts.

### Why predict at all?

Experimental structures (X-ray, cryo-EM, NMR) are slow and expensive. The PDB holds about 230,000 structures,
while UniProt holds about 250 million sequences. Predicted structures help with:

- designing mutants and interpreting disease variants,
- finding binding pockets for drug discovery,
- solving crystallography data by molecular replacement, and fitting cryo-EM maps,
- guessing the function of uncharacterised genes.

### A short history of methods

| Era | Approach | Idea | Typical accuracy |
|---|---|---|---|
| 1970s | **Statistical secondary structure** (Chou–Fasman, GOR) | Residue propensities for helix and strand | Q3 ≈ 50–65 % |
| 1990s–2000s | **Profile / neural-net secondary structure** (PSIPRED) | Use evolutionary profiles from multiple alignments | Q3 ≈ 80 % |
| 1990s– | **Homology (comparative) modelling** (MODELLER, SWISS-MODEL) | Copy the backbone of a related protein of known structure (>30 % identity) | Good when a template exists |
| 2000s | **Threading / fold recognition** (Phyre2, I-TASSER) | Fit the sequence onto known folds even at low identity | Variable |
| 2000s–2010s | **Ab initio / fragment assembly** (Rosetta) | Assemble 3–9-residue fragments while minimising an energy function | Small proteins only |
| 2010s | **Co-evolution contacts** (EVfold, RaptorX) | Residues that mutate together are usually in contact | Big jump for large families |
| 2020 → | **End-to-end deep learning** (AlphaFold2, RoseTTAFold) | A network reads the MSA and templates and outputs 3D coordinates directly | Often close to experimental accuracy |
| 2022 → | **Protein language models** (ESMFold) | Learns co-evolution implicitly from single sequences; no MSA, very fast | Slightly below AF2 |
| 2024 → | **Complexes and ligands** (AlphaFold3, RoseTTAFold All-Atom, Boltz, Chai) | Diffusion-based all-atom prediction of proteins with DNA, RNA, ligands and ions | Active research |

The yardstick is **CASP**, a blind competition held every two years. In CASP14 (2020), AlphaFold2 reached a median
GDT_TS of about 92 on free-modelling targets, which is roughly the error level of experimental structures. The
work was recognised with the 2024 Nobel Prize in Chemistry (Hassabis and Jumper, shared with David Baker for
protein design).

---

## Part 1 — Find the gene in the DNA

`triplets.txt` contains about 10 kb of genomic DNA, printed as codons in all three forward reading frames
(`start 1`, `start 2`, `start 3`). The first job is to find which protein it encodes.

```bash
cd tutorials/protein-structure-prediction
python3 structure_tutorial.py orfs
```

This translates each frame and lists open reading frames (ORFs): stretches from a Met (ATG) to a stop codon
that are at least 50 residues long.

```
frame  start    end length  first 40 aa
    2   3096   3425    330  MPLAGQGCRDTPPFVQGIPHSPFSPANTTLALLSRGRHTP
    1   2686   2884    199  MQEGGADGRRRKEGRKQGRKEGLLEPSHPGTVGRGDCRPS
    3    492    615    124  MQTSASAPGPGIMSSRRGVTSVLGASSVLLPDHLHVSSVW
  ...
    1   2077   2155     79  MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEVSL
    1   2327   2403     77  MVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQ
```

> **Question 1.1.** The longest ORF is the obvious candidate. Why is "longest ORF" a poor way to find genes
> in *eukaryotic genomic* DNA?

Copy the first 40 residues of each ORF into **BLASTP** (<https://blast.ncbi.nlm.nih.gov/>, database
*UniProtKB/Swiss-Prot*). Most ORFs give no convincing hit. The frame 1 ORF at codon 2077 gives a near-perfect match:

> `MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIE` → **GTPase HRas (UniProt P01112), human**

The residues right after `...DPTIE` don't match H-Ras (`VSL` instead of `DSYRK`). The gene continues in another
piece: this is an **intron** boundary. The ORF at 2327 matches a later part of H-Ras, and the C-terminus turns up
in frame 2.

### Splicing the exons

Take the H-Ras protein from UniProt and use **tBLASTn** (protein against your DNA) to map each part onto
`triplets.txt`. You'll find four coding exons. At every boundary the intron starts with **GT** and ends with
**AG** (the canonical splice sites). The script contains the answer:

```bash
python3 structure_tutorial.py splice
```

```
exon 1:  6229-6339  (111 nt)   intron  267 nt  GT...AG
exon 2:  6607-6785  (179 nt)   intron  153 nt  GT...AG
exon 3:  6939-7098  (160 nt)   intron  697 nt  GT...AG
exon 4:  7796-7915  (120 nt)

CDS 570 nt -> 189 aa, ends with stop: True

MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAG
QEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDL
AARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPG
CMSCKCVLS
```

That is the full 189-residue H-Ras protein: a small GTPase and molecular switch, and one of the most frequently
mutated oncogenes (G12, G13 and Q61 hotspots). We'll predict its structure and check the prediction against
the crystal structure **PDB 5P21**.

> **Question 1.2.** Exons 1–3 are all in reading frame 1, but exon 4 is in frame 2. Also, exon 2 ends in the
> *middle* of codon 97 (after its second base). Look at the intron lengths. Why does the reading frame
> change, and why does a "translate each frame separately" approach miss this?

Save the protein for later:

```bash
python3 structure_tutorial.py fasta --spliced > hras.fasta
```

---

## Part 2 — What the sequence alone tells you

Before any 3D modelling, simple sequence statistics already say a lot about a protein's likely structure.

```bash
python3 structure_tutorial.py analyze --spliced
```

### 2a. Global properties

```
Length            189 aa
Mol. weight       21.3 kDa
GRAVY             -0.42   (>0 hydrophobic, <0 hydrophilic)
Net charge ~pH7   -6
Disorder-promoting residues (PESQKAG): 41%
Order-promoting residues   (WCFIYVLN): 35%
```

- **GRAVY** (grand average of hydropathy) below 0 means a soluble, cytoplasmic protein.
- A balance of order- and disorder-promoting residues suggests a mostly folded protein.

### 2b. Hydropathy and transmembrane helices

The Kyte–Doolittle scale gives each residue a hydrophobicity score. Averaged over a 19-residue window (about
the length of a helix that spans a membrane), values above about 1.6 suggest a **transmembrane helix**.

```
No transmembrane helix predicted (no 19-residue window with KD > 1.6).
Max windowed hydropathy: +1.59
```

The peak is close to the cutoff but below it. That fits the biology: Ras is anchored to the membrane by a
**lipid** (farnesyl) attached to the C-terminal CAAX motif `CVLS`, not by a transmembrane helix.

### 2c. Secondary structure: Chou–Fasman

The script uses a simplified Chou–Fasman method. It averages each residue's statistical tendency to be in a
helix (P_α) or a strand (P_β) over a 7-residue window, then removes very short segments. You can read the whole
method in about 30 lines of `chou_fasman()` in the script.

For H-Ras, the script compares the prediction with the crystal structure. The reference string is built from
approximate helix and strand boundaries in PDB 5P21; you can recompute it exactly with DSSP in Part 4.

```
    1 pred CCCCEEEEEEECCCCCCCCEEEEEEEECCCCCCCCCCCCCCEEEEEECCCCCEEECCCHH
      5P21 CEEEEEEEECCCCCCHHHHHHHHHHCCCCCCCCCCCEEEEEEEEEECCEEEEEEEEEECC
   61 pred HHHHHHHHHHCCCCCCCEEEEEEEECCCCHHHHHHHHHHHHHHCCCCCCEEEEEECCCHH
      5P21 CCCCCHHHHHHHHHCCEEEEEEECCCHHHHHHHHHHHHHHHHHHCCCCCCEEEEEECCCC
  121 pred HHHHHHHHHHHHHCCCCEEEEECCCCCCCCCHHHHHEEEEEHHHHH
      5P21 CCCCCCHHHHHHHHHHHCCCEEECCCCCCCCHHHHHHHHHHHHHHH
Q3 accuracy (fraction of residues with the correct state): 57%
```

**Q3** is the fraction of residues assigned the correct state out of three (H, E, C). A random guess gives about
33–40 %. Chou–Fasman gets **57 %**. It finds β1, β5, α3 and α5 well, but it misses α1 (residues 16–25, part of the
P-loop region) and calls part of it a strand.

> **Question 2.1.** Chou–Fasman looks at each residue in isolation, within a small window. Name two kinds of
> information it ignores that modern predictors use. *(Hint: think about evolution, and about residues far
> apart in sequence.)*

> **Question 2.2.** Try `python3 structure_tutorial.py analyze --seq <sequence>` on a protein you know, such as
> a coiled-coil (for example the GCN4 leucine zipper `RMKQLEDKVEELLSKNYHLENEVARLKKLVGER`) or an all-β
> protein. Where does the simple method do well, and where does it fail?

### 2d. A negative control: a sequence that probably doesn't fold

Now look at the longest ORF, which BLAST didn't recognise:

```bash
python3 structure_tutorial.py analyze --frame 2
```

```
Disorder-promoting residues (PESQKAG): 50%
Most common       P 18%, L 14%, S 12%, G 11%, T 10%
Helix 5%   Strand 2%   Coil 93%
Best approximate period: 28 aa (75% of residues match the residue 28 positions later)
  -> strongly repetitive; expect low-confidence / disordered prediction
```

This ORF is 18 % proline and made of near-perfect 28-residue tandem repeats (`...PHSPFSPGDATLPLLSRGRHT...`). It
almost certainly comes from a repetitive stretch of noncoding DNA in this region, not a real protein. Keep it: in
Part 4 it shows what a predictor does with a sequence that has no real structure.

---

## Part 3 — Predicting the 3D structure

Choose **one** of these tools; all three are free. For a 189-residue protein, each takes a few minutes or less.

| Tool | How | Needs MSA? | Speed | Best for |
|---|---|---|---|---|
| **ESMFold** | [ESM Metagenomic Atlas "Fold sequence"](https://esmatlas.com/resources?action=fold), or `python3 structure_tutorial.py fold --spliced` | No | Seconds | Quick single chains of 400 residues or fewer |
| **ColabFold** (AlphaFold2) | [ColabFold AlphaFold2 notebook](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb): paste the sequence, then *Runtime → Run all* | Yes (MMseqs2 server) | Minutes | Best accuracy for single chains and simple complexes; gives the PAE plot |
| **AlphaFold Server** (AlphaFold3) | <https://alphafoldserver.com> (free Google account) | Yes (automatic) | Minutes | Proteins **with ligands, ions and nucleic acids** (for example Ras + GTP + Mg²⁺) |

You can also skip prediction and download the precomputed model from the **AlphaFold Protein Structure Database**:
<https://alphafold.ebi.ac.uk/entry/P01112>.

### Running ESMFold from the script

```bash
python3 structure_tutorial.py fold --spliced --out hras_esmfold.pdb
```

This sends the sequence to the public ESMFold API, saves the PDB file and prints a confidence summary. If your
network blocks the API, the script says so; use one of the web options above instead.

### What the network is doing (conceptually)

```
sequence ──► MSA search ──► Evoformer (48 blocks) ──► Structure module ──► 3D coordinates
            (homologs)      pair + MSA reasoning:       rotates/translates     + pLDDT
                            "which residues touch?"     each residue frame     + PAE
                                                        (recycled 3×)
```

- The **MSA** (multiple sequence alignment) provides co-evolution: if residue *i* mutates, residue *j*
  compensates, so they are probably in contact.
- **ESMFold** replaces the MSA with a protein language model trained on millions of sequences. It is faster,
  but it struggles more with proteins that have few relatives.
- **AlphaFold3** replaces the structure module with a **diffusion** model that generates atoms from noise, which
  lets it place ligands and nucleic acids too.

---

## Part 4 — Interpreting the prediction

A predicted structure is a **hypothesis**. Always check its confidence scores before believing it.

### 4a. pLDDT: per-residue confidence

Predictors store pLDDT (predicted local distance difference test, 0–100) in the **B-factor column** of the PDB
file. Summarise any predicted model with:

```bash
python3 structure_tutorial.py plddt hras_esmfold.pdb     # or the AlphaFold DB / ColabFold PDB
```

| pLDDT | Meaning | Typical interpretation |
|---|---|---|
| ≥ 90 | Very high | Backbone *and* side chains are reliable |
| 70–90 | Confident | Backbone is reliable; side chains less so |
| 50–70 | Low | Treat with caution; possibly flexible |
| < 50 | Very low | Probably **disordered**; the shape shown is meaningless |

In ChimeraX, colour by confidence with `color bfactor palette alphafold`. In Mol\*, choose *Coloring → pLDDT
Confidence*.

**What to expect for H-Ras:** residues 1–166 (the G domain) should be almost entirely above 90. The
**hypervariable region** (about 167–189, ending in `CVLS`) should drop below 50. This region is flexible in the
cell and isn't resolved in crystal structures either. Low pLDDT often correctly flags *real* disorder.

> **Question 4.1.** Fold the repetitive frame 2 ORF from Part 2d (`analyze --frame 2`, then `fasta --frame 2`)
> with ESMFold. What mean pLDDT do you get? Does the model look like a compact globule, or like "spaghetti"?
> What does that tell you about the protein?

### 4b. PAE: confidence in relative positions

The **predicted aligned error** (PAE) matrix (ColabFold and AlphaFold DB show it as a heat map) gives the
expected position error of residue *j* when the model is superimposed on residue *i*.

- Low PAE (dark) within a block means that domain is confidently folded.
- High PAE between two blocks means each domain may be right, but **their relative orientation is unknown**.
- For complexes, **ipTM** (interface pTM) above about 0.8 suggests a confident interface; below about 0.6, the
  interface is likely wrong.

> **Question 4.2.** A two-domain protein has pLDDT above 90 everywhere, but a high PAE between the domains. Can
> you use the model to measure the distance between two residues in different domains? Why or why not?

### 4c. Compare with experiment

Download the crystal structure **5P21** (H-Ras 1–166 bound to the GTP analogue GppNHp and Mg²⁺) from
<https://www.rcsb.org/structure/5P21>. Superimpose it on your model:

**ChimeraX**

```
open 5P21
open hras_esmfold.pdb
matchmaker #2 to #1            # prints Cα RMSD
color bfactor #2 palette alphafold
```

**Or online:** the RCSB *Pairwise Structure Alignment* tool (<https://www.rcsb.org/alignment>) reports RMSD and
**TM-score**.

| Metric | Range | Rule of thumb |
|---|---|---|
| Cα **RMSD** | 0 Å and up | < 2 Å is very close; sensitive to floppy tails, so restrict it to residues 1–166 |
| **TM-score** | 0–1 | > 0.5 means the same fold; > 0.9 is nearly identical; much less sensitive to protein length |
| **GDT_TS** | 0–100 | The CASP metric; > 90 is competitive with experiment |

> **Question 4.3.** Compute the RMSD twice: over all residues, and over residues 1–166 only. Explain the
> difference.

---

## Part 5 — What structure predictors can't tell you

H-Ras is a good case study for the limits, because its biology depends on details that a single static model
misses.

1. **Conformational states.** Ras switches between an inactive GDP-bound state and an active GTP-bound state.
   The difference is concentrated in two loops, *switch I* (residues about 30–38) and *switch II* (about
   59–76). AlphaFold2 gives one conformation, usually close to whichever state dominates the PDB.
2. **Point mutations.** Try predicting the oncogenic **G12V** mutant (change residue 12 from G to V). The model
   will look almost identical to wild type. The mutation is harmful because it blocks GTP hydrolysis (it sterically
   hinders the catalytic arginine of GAP), not because it changes the fold. Predictors are **not** reliable
   estimators of how mutations affect stability or function.
3. **Ligands, ions and modifications.** AlphaFold2 has no GTP, no Mg²⁺ and no farnesyl lipid. AlphaFold3 can
   add GTP and Mg²⁺, but most post-translational modifications are still out of reach.
4. **Membranes and partners.** Ras works at the membrane, bound to effectors (Raf, PI3K) and regulators (SOS,
   GAP). Single-chain predictions don't show these interactions.
5. **Disorder is not a failure.** Low pLDDT regions such as the hypervariable region may be *functionally*
   disordered. About 30 % of human protein residues fall in intrinsically disordered regions.
6. **Dynamics and folding pathways.** Predictors give a final structure, not how the protein moves or how
   it folds. For motion, use molecular dynamics or experiments (for example NMR or HDX).

> **Question 5.1.** Use AlphaFold Server to predict H-Ras 1–166 **with GTP and Mg²⁺**. Then predict it
> **with GDP**. Superimpose the two and colour by residue number. Where are the differences? Do they
> match switch I and switch II?

---

## Exercises (to hand in)

1. **Gene to protein.** Report the four exon coordinates and the splice-site dinucleotides you found with
   tBLASTn. Explain in two sentences why "longest ORF" failed to find this gene.
2. **Classical versus modern.** Make a table of Q3 for Chou–Fasman (from the script), PSIPRED
   (<http://bioinf.cs.ucl.ac.uk/psipred/>) and the secondary structure from your AlphaFold/ESMFold model (in
   ChimeraX: `dssp` then `info residues attribute ss_type`). Explain the ranking.
3. **Confidence.** Include a figure of your H-Ras model coloured by pLDDT. Mark the G domain, the switch
   regions and the hypervariable region. State the mean pLDDT for residues 1–166 and for 167–189.
4. **Validation.** Report Cα RMSD and TM-score between your model and 5P21 (residues 1–166).
5. **Negative control.** Report the mean pLDDT for the frame 2 repeat ORF, and argue whether it is a real protein.
6. **Limits.** In about 150 words, explain why an accurate AlphaFold model of H-Ras is still not enough to
   design a drug against KRAS-G12C. *(Hint: look up sotorasib and the switch II pocket.)*

---

## Hints for the in-text questions

<details>
<summary>Click to expand</summary>

- **1.1** Introns interrupt the coding sequence, and chance ORFs in noncoding DNA (especially GC-rich or
  repetitive DNA) can be longer than any single real exon. Gene finders use splice-site models, codon bias and
  homology.
- **1.2** Introns 1 and 2 are 267 and 153 nt long (multiples of 3), so exons 2 and 3 stay in frame 1, even
  though codon 97 (Arg, AGG) is split AG|G across intron 2. Intron 3 is 697 nt long (not a multiple of 3), so exon 4
  (residues 151–189) moves to frame 2. Only splicing puts the pieces back in one frame.
- **2.1** Evolutionary information (MSAs and profiles) and long-range contacts. β-sheets pair strands that
  can be far apart in sequence, which is exactly what a local window can't see.
- **4.1** Expect mean pLDDT well below 50 and an extended, non-compact model. That points to a disordered or
  non-real protein.
- **4.2** No: PAE says the relative placement of the domains is uncertain, even though each domain is
  confident on its own.
- **4.3** The unresolved tail (167–189) is placed arbitrarily and inflates the RMSD. A superposition restricted
  to the folded domain typically gives about 1 Å or less.

</details>

---

## Further reading

- Anfinsen, C. B. (1973). Principles that govern the folding of protein chains. *Science* 181, 223–230.
- Chou, P. Y. & Fasman, G. D. (1978). Prediction of the secondary structure of proteins from their amino acid
  sequence. *Adv. Enzymol.* 47, 45–148.
- Kyte, J. & Doolittle, R. F. (1982). A simple method for displaying the hydropathic character of a protein.
  *J. Mol. Biol.* 157, 105–132.
- Pai, E. F. et al. (1990). Refined crystal structure of the triphosphate conformation of H-ras p21 at 1.35 Å
  resolution. *EMBO J.* 9, 2351–2359. *(PDB 5P21)*
- Jumper, J. et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature* 596, 583–589.
- Lin, Z. et al. (2023). Evolutionary-scale prediction of atomic-level protein structure with a language model.
  *Science* 379, 1123–1130. *(ESMFold)*
- Mirdita, M. et al. (2022). ColabFold: making protein folding accessible to all. *Nat. Methods* 19, 679–682.
- Abramson, J. et al. (2024). Accurate structure prediction of biomolecular interactions with AlphaFold 3.
  *Nature* 630, 493–500.
