# Tutorial: Protein Structure Prediction

**From a protein sequence to a 3D model, with Ras and haemoglobin as worked examples**

| | |
|---|---|
| **Audience** | Undergraduate bioinformatics (BIO116). Assumes you know the amino acids and the basics of protein structure. |
| **Time** | About 2–3 hours (Parts 1–3 offline, about 1.5 h; Parts 4–5 need a browser and about 30 min of compute) |
| **You need** | Python 3.8+ (standard library only), a web browser, and optionally [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) or the [Mol\* viewer](https://molstar.org/viewer/) |
| **Files** | `protein-structure-tutorial.html` (interactive version: open it in any browser, works offline), `structure_tutorial.py` (companion script), `data/` (reference structures and secondary-structure datasets) |

The starting point is simple: you have a protein sequence and want to know its structure. We use two proteins
whose structures are known experimentally, so every prediction can be checked residue by residue:

- **H-Ras**, a GTPase and cancer driver with an α/β fold (a six-stranded β-sheet surrounded by helices);
- **haemoglobin**, the oxygen carrier, whose α and β chains are all-α globins (eight helices, no strands).

By the end you will be able to:

1. Make and critique sequence-only predictions: hydropathy and secondary structure (Chou–Fasman).
2. Train a small neural network for secondary structure and explain why it beats Chou–Fasman, and what beat it.
3. Map a prediction onto the experimental structure and see which parts were right.
4. Run a modern predictor (AlphaFold2 through ColabFold, ESMFold, or AlphaFold Server) and compare its model with
   the experimental structure.
5. Read the confidence metrics (pLDDT, PAE) and explain what predictors **cannot** tell you.

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
- guessing the function of uncharacterised proteins.

### A short history of methods

| Era | Approach | Idea | Typical accuracy |
|---|---|---|---|
| 1970s | **Statistical secondary structure** (Chou–Fasman, GOR) | Residue propensities for helix and strand | Q3 ≈ 50–65 % |
| 1988 | **Neural networks on single sequences** (Qian & Sejnowski) | Learn the rules from known structures instead of hand-made tables | Q3 ≈ 64 % |
| 1990s–2000s | **Profile neural networks** (PHD, PSIPRED) | Feed the network evolutionary profiles from multiple alignments | Q3 ≈ 71–80 % |
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

## Part 1 — The two proteins

```bash
cd tutorials/protein-structure-prediction
python3 structure_tutorial.py proteins
```

```
hras  H-Ras (GTPase HRas), human (UniProt P01112), 189 aa
      189 aa; G domain 1-166 (alpha/beta fold) + flexible C-terminal tail. Reference: PDB 1CRR, NMR model 1 (wild type, GDP-bound), residues 1-166.
hbb   Haemoglobin subunit beta, human (UniProt P68871), 146 aa
      146 aa; all-alpha globin fold, binds one haem; part of the alpha2beta2 tetramer. Reference: PDB 4HHB chain B (deoxyhaemoglobin, X-ray 1.74 A).
hba   Haemoglobin subunit alpha, human (UniProt P69905), 141 aa
      141 aa; globin fold, 43% identical to beta. Reference: PDB 4HHB chain A.
```

| | H-Ras | Haemoglobin |
|---|---|---|
| Job | Molecular switch: on with GTP, off with GDP | Carries O₂ in red blood cells |
| Fold | α/β: six-stranded β-sheet, five helices | All-α globin: eight helices (A–H), no strands |
| Chains | One chain, 189 residues (G domain 1–166 + membrane-anchoring tail) | Tetramer α₂β₂ (α 141, β 146 residues), one haem per chain |
| Famous mutations | G12, G13, Q61 (lock Ras "on" in cancer) | β E6V (sickle-cell disease) |
| Reference structure here | PDB 1CRR, NMR, wild type + GDP (Kraulis et al. 1994) | PDB 4HHB, X-ray 1.74 Å, deoxy form (Fermi et al. 1984) |

Save the sequences for later:

```bash
python3 structure_tutorial.py fasta --protein hras > hras.fasta
python3 structure_tutorial.py fasta --protein hbb  > hbb.fasta
```

Where the reference data come from: `data/make_reference.py` rebuilds `data/structures/` and
`data/reference_ss.tsv` from the public files, and documents every source. The secondary structure used as the
"right answer" is DSSP (Kabsch & Sander 1983) computed from these coordinates: H (α-helix), G (3₁₀ helix) and I (π
helix) count as helix; E (strand) and B (bridge) as strand; everything else as coil.

---

## Part 2 — What the sequence alone tells you

```bash
python3 structure_tutorial.py analyze --protein hras
python3 structure_tutorial.py analyze --protein hbb
```

### 2a. Global properties and hydropathy

```
== H-Ras (GTPase HRas), human (UniProt P01112) ==
Length            189 aa
Mol. weight       21.3 kDa
GRAVY             -0.42   (>0 hydrophobic, <0 hydrophilic)
Net charge ~pH7   -6
...
No transmembrane helix predicted (no 19-residue window with KD > 1.6).
Max windowed hydropathy: +1.59
```

- **GRAVY** (grand average of hydropathy) below 0 means a soluble protein. Both examples are soluble.
- The Kyte–Doolittle scale gives each residue a hydrophobicity score. Averaged over a 19-residue window (about
  the length of a helix that crosses a membrane), values above about 1.6 suggest a **transmembrane helix**.
  H-Ras peaks just below the cutoff (+1.59). It is anchored to the membrane by a **lipid** attached to its
  C-terminal `CVLS` motif, not by a transmembrane helix. Haemoglobin β peaks at +0.99.

### 2b. Secondary structure: Chou–Fasman

The script uses a simplified Chou–Fasman method. It averages each residue's statistical tendency to be in a
helix (P_α) or a strand (P_β) over a 7-residue window, then removes very short segments. You can read the whole
method in about 30 lines of `chou_fasman()` in the script. The output compares it with the experimental
structure; **lower-case letters are wrong**.

```
    1 seq  MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAG
      CF   CCccEEEEEeeCCCCcccceeeeeeeeCCCCCCCCCCCcccEEeeeeCCCCcEEEcCChh
      DSSP CCEEEEEEECCCCCCHHHHHHHHHCCCCCCCCCCCCCCEEEEECCCCCCCCEEEEECCCC
...
Chou-Fasman: helix 29%, strand 27%, coil 44%
Experiment:  helix 35%, strand 23%, coil 42%  (residues 1-166)
Q3 (fraction of residues in the correct state): 60%
```

For H-Ras, Chou–Fasman gets **60 %** of residues right (random guessing gives about 33–40 %). It finds β1 but
calls most of the first helix (residues 16–24) coil or strand. Now haemoglobin β:

```
    1 seq  VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV
      CF   CCChHHHHHHHeeeeeeeCcccccccccceeeeeeeeeeHHhHHHCCCCCcccccCCccc
      DSSP CCCCHHHHHHHHHHHCCCCHHHHHHHHHHHHHHHCHHHHHHCHHHCCCCCHHHHHCCHHH
...
Chou-Fasman: helix 25%, strand 26%, coil 49%
Experiment:  helix 78%, strand 0%, coil 22%  (residues 1-146)
Q3 (fraction of residues in the correct state): 39%
```

Only **39 %**. Chou–Fasman calls a quarter of haemoglobin strand, although the protein has no strands at all.

> **Question 2.1.** Why does Chou–Fasman mistake haemoglobin's helices for strands? *(Hint: look at which
> residues fill the buried face of a globin helix, and at their P_β values in the script.)*

> **Question 2.2.** Chou–Fasman looks at each residue in isolation, within a small window. Name two kinds of
> information it ignores that modern predictors use.

### 2c. Neural networks: learning the rules from data

Chou and Fasman tabulated their propensities by hand. In 1988, Qian and Sejnowski instead trained a **neural
network** on known structures. The companion script builds the same kind of network:

```
 window of 13 residues        one-hot input           hidden layer      output (softmax)
 ... L T I Q [L] I Q N ...  ->  13 x 22 units   ->   10 tanh units  ->   P(helix), P(strand), P(coil)
```

- Each residue in the window becomes 22 inputs (20 amino acids, unknown, "past the chain end"), exactly one of
  them on. This is **one-hot encoding**.
- The network is trained by **gradient descent**: for each residue it nudges every weight to make the true
  state (from DSSP) a little more likely.
- **Data:** the CB6133-filtered set for training and the standard **CB513** benchmark (514 proteins, 84,765
  residues) for testing; the two share no proteins. Seventeen training proteins related to Ras or haemoglobin
  (≥ 25 % identity, or a small-GTPase motif pair) are removed, so the network has never seen either family. See
  `data/README.md` and `data/excluded_training.txt`.

```bash
python3 structure_tutorial.py nn        # 300 training proteins, about 15 seconds
```

```
Training set: 300 proteins, 68263 residues (17 relatives of Ras and haemoglobin excluded)
Test set (CB513): 514 proteins, 84765 residues
Network: window 13 x 22 inputs -> 10 hidden units -> 3 outputs

epoch  train Q3  test Q3   time
    1     60.7%    61.4%     3s
    2     62.4%    61.5%     4s
    3     62.9%    62.3%     4s

Chou-Fasman on the same test set: 53.9%
...
Haemoglobin subunit beta, human (no relatives in the training set); lower case = disagrees with the experimental structure
    1 seq  VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV
      CF   CCChHHHHHHHeeeeeeeCcccccccccceeeeeeeeeeHHhHHHCCCCCcccccCCccc
      NN   CCCCccHHHHHHHHHCCCCccHcccHHHHHeeeeecHHHHehHccCCCCCccHecCCccH
      DSSP CCCCHHHHHHHHHHHCCCCHHHHHHHHHHHHHHHCHHHHHHCHHHCCCCCHHHHHCCHHH
...
Q3: Chou-Fasman 39%, neural network 58%
```

On the same 84,765 test residues, the network beats Chou–Fasman by about 8 points after training on only 300
proteins. On H-Ras it gets 63 % (Chou–Fasman 60 %), and on haemoglobin β 58 % (Chou–Fasman 39 %). Try other
settings:

| Command | CB513 Q3 (84,765 residues) | H-Ras | Haemoglobin β |
|---|---|---|---|
| `nn --window 1 --hidden 0` (one residue, linear) | 48.7 % | 49 % | 53 % |
| `nn` (window 13, 10 hidden, 300 proteins, 3 epochs) | 62.3 % | 63 % | 58 % |
| `nn --window 17 --hidden 20 --proteins 2000 --epochs 4` (about 2 minutes) | 65.2 % | 72 % | 54 % |
| Chou–Fasman, for comparison | 53.9 % | 60 % | 39 % |

Accuracy on a single protein moves by several points from run to run and between settings, because one protein
is only 146–166 residues. The CB513 column, averaged over 514 proteins, is the reliable measure.

Three lessons:

1. **Context matters.** A single residue gives under 50 %, barely better than always guessing coil (43 % of CB513).
   A window of 13–17 residues adds about 15 points.
2. **Model size isn't the bottleneck.** A bigger network with 7 times more data gains only about 3 points on
   CB513. A single sequence doesn't carry much more information about local structure, and single-sequence
   methods level off around 65 %.
3. **Evolution breaks the barrier.** In 1993, **PHD** (Rost and Sander) fed the network a *profile* from a
   multiple sequence alignment instead of a single sequence and reached 70.8 %. **PSIPRED** (1999) reached about
   76.5 % with PSI-BLAST profiles, and deep networks with profiles reach about 85 % on CB513 today. The practical
   ceiling is roughly 88–90 %, because DSSP assigns slightly different states to different structures of the same
   protein. (Published values come from different test sets, so compare them loosely.)

> **Question 2.3.** Train with `--proteins 50 --epochs 8`. Compare the training and test Q3. What is happening, and
> why must a predictor always be tested on proteins it has never seen?

> **Question 2.4.** A position in a profile says "this column is hydrophobic in 95 % of homologues". Why is that
> more useful for predicting a buried β-strand than the single residue at that position?

---

## Part 3 — Which parts were predicted correctly?

Numbers such as "Q3 = 58 %" hide *where* a method goes wrong. Open the interactive page
(`protein-structure-tutorial.html`, step 4). It shows the experimental structure of the selected protein and
colours every residue green (predicted correctly) or red (predicted wrongly), for Chou–Fasman or for the network
you trained in step 3. For haemoglobin, the other three chains of the tetramer are drawn faintly and the haem irons
as red spheres. Hover over a residue to see its experimental and predicted state.

What to look for:

- **Haemoglobin β, Chou–Fasman:** whole helices turn red, because it calls them strands (38 residues wrongly
  called strand).
- **Haemoglobin β, neural network:** much of the red disappears. The remaining errors are short false strands and
  helices that end too early or start too late.
- **H-Ras:** both methods get only part of the β2–β3 hairpin (residues 37–58). What makes a strand a strand is its
  partner strand, which can be far away in the sequence, so local methods find strands harder than helices.

> **Question 3.1.** In the viewer, are the network's errors on H-Ras mostly inside secondary-structure elements or at
> their ends? What does that suggest about how to score a predictor fairly?

---

## Part 4 — Predicting the 3D structure

Choose **one** of these tools; all are free. For proteins of this size, each takes a few minutes or less.

| Tool | How | Needs MSA? | Speed | Best for |
|---|---|---|---|---|
| **ESMFold** | [ESM Metagenomic Atlas "Fold sequence"](https://esmatlas.com/resources?action=fold), or `python3 structure_tutorial.py fold --protein hbb` | No | Seconds | Quick single chains of 400 residues or fewer |
| **ColabFold** (AlphaFold2) | [ColabFold AlphaFold2 notebook](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb): paste the sequence, then *Runtime → Run all*. Join chains with `:` for a complex | Yes (MMseqs2 server) | Minutes | Best accuracy for single chains and simple complexes; gives the PAE plot |
| **AlphaFold Server** (AlphaFold3) | <https://alphafoldserver.com> (free Google account) | Yes (automatic) | Minutes | Proteins **with ligands, ions and nucleic acids**: Ras + GDP + Mg²⁺, or the haemoglobin tetramer (2 α + 2 β) + 4 haem |

You can also download precomputed models from the **AlphaFold Protein Structure Database**:
[P01112](https://alphafold.ebi.ac.uk/entry/P01112) (H-Ras), [P68871](https://alphafold.ebi.ac.uk/entry/P68871)
(haemoglobin β), [P69905](https://alphafold.ebi.ac.uk/entry/P69905) (haemoglobin α). The database models of
haemoglobin include the initiator methionine (147 and 142 residues); the comparison below handles that by aligning
the sequences first.

### What the network is doing (conceptually)

```
sequence ──► MSA search ──► Evoformer (48 blocks) ──► Structure module ──► 3D coordinates
            (homologs)      pair + MSA reasoning:       rotates/translates     + pLDDT
                            "which residues touch?"     each residue frame     + PAE
                                                        (recycled 3×)
```

- The **MSA** (multiple sequence alignment) provides co-evolution: if residue *i* mutates, residue *j*
  compensates, so they are probably in contact. Globins and Ras-family GTPases both have thousands of relatives.
- **ESMFold** replaces the MSA with a protein language model trained on millions of sequences. It is faster,
  but it struggles more with proteins that have few relatives.
- **AlphaFold3** replaces the structure module with a **diffusion** model that generates atoms from noise, which
  lets it place ligands such as haem and GTP too.

---

## Part 5 — Judging the model

### 5a. Compare the model with the experiment

```bash
python3 structure_tutorial.py compare my_hbb_model.pdb --protein hbb
```

`compare` aligns the model's sequence with the reference, superimposes the two structures on their Cα atoms,
re-fits on the residues within 4 Å (so a floppy loop can't drag the whole superposition), and reports how far each
residue is from where it should be. Without a model of your own yet, try a stand-in: the experimental **α chain**
used as a "model" of the β chain. The two are 43 % identical, which is the situation homology modelling exploits.

```bash
python3 structure_tutorial.py compare data/structures/hb_4hhb.pdb --chain A --protein hbb
```

```
Model: data/structures/hb_4hhb.pdb, chain A, 141 residues
Reference: Haemoglobin subunit beta, human, hb_4hhb.pdb chain B, 146 residues
Matched residues: 136 (67 identical in sequence, 49%)
Cα RMSD over all matched residues: 2.08 A
Cα RMSD over the 125 residues within 4 A: 1.56 A
Fraction of reference residues within  1 A: 38%   2 A: 73%   4 A: 86%   8 A: 93%
GDT-style score (mean of the four fractions): 72.4

Per residue:  # within 2 A (correct)   + 2-4 A   . more than 4 A   - not in model
    1 ref  VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV
      fit  +-+++++#######################..-############-.+.-----+#####

   61 ref  KAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGK
      fit  #################+.++..++.-#####+###+#########..++.-#######+

  121 ref  EFTPPVQAAYQKVVAGVANALAHKYH
      fit  +#########################
```

Three quarters of β is within 2 Å of the α-chain template. The unmatched stretch (`-----`) around residues 45–55 is
where α chains are shorter: they lack the D helix. The interactive page (step 6) does the same comparison in 3D:
load a PDB file, or press the demo button, and the experimental chain is coloured green (≤ 2 Å), amber (2–4 Å) or
red (> 4 Å), with the model drawn on top as a thin line.

| Metric | Range | Rule of thumb |
|---|---|---|
| Fraction within 2 Å | 0–100 % | The simplest "how much is right" measure |
| Cα **RMSD** | 0 Å and up | < 2 Å is very close; one floppy tail can inflate it, so read it together with the fractions |
| **GDT**-style score | 0–100 | Mean of the fractions within 1, 2, 4 and 8 Å; CASP's GDT_TS takes the best superposition for each cutoff |
| **TM-score** | 0–1 | > 0.5 means the same fold; use the [RCSB pairwise alignment tool](https://www.rcsb.org/alignment) |

> **Question 5.1.** Compare your AlphaFold or ESMFold model of H-Ras with 1CRR. Which residues fall outside 2 Å?
> 1CRR is an NMR structure: why might a good model still differ from it in a few loops?

### 5b. pLDDT: per-residue confidence

Predictors store pLDDT (predicted local distance difference test, 0–100) in the **B-factor column** of the PDB
file:

```bash
python3 structure_tutorial.py plddt my_model.pdb
```

| pLDDT | Meaning | Typical interpretation |
|---|---|---|
| ≥ 90 | Very high | Backbone *and* side chains are reliable |
| 70–90 | Confident | Backbone is reliable; side chains less so |
| 50–70 | Low | Treat with caution; possibly flexible |
| < 50 | Very low | Probably **disordered**; the shape shown is meaningless |

**What to expect:** for H-Ras, residues 1–166 should be almost entirely above 90 and the C-terminal tail
(167–189) below 50: it is flexible in the cell and missing from experimental structures. For haemoglobin chains,
expect high pLDDT across all eight helices. In the interactive page, switch the viewer to "Model confidence" to
colour the structure by pLDDT; in ChimeraX use `color bfactor palette alphafold`.

> **Question 5.2.** Is pLDDT high where your model is close to the experiment (step 6, "Distance from experiment")?
> Find a residue where the two disagree and suggest why.

### 5c. PAE: confidence in relative positions

The **predicted aligned error** (PAE) matrix (ColabFold and AlphaFold DB show it as a heat map) gives the expected
position error of residue *j* when the model is superimposed on residue *i*. Low PAE within a block means that part
is confidently folded; high PAE between two blocks means their relative orientation is uncertain. For complexes
such as the haemoglobin tetramer, **ipTM** above about 0.8 suggests confident interfaces.

> **Question 5.3.** You predict the haemoglobin tetramer. All four chains have pLDDT above 90, but PAE between the
> α and β chains is high. Can you trust the α/β interface? Why or why not?

---

## Part 6 — What structure predictors can't tell you

1. **Conformational states.** Ras switches between an inactive GDP-bound and an active GTP-bound state, with the
   changes concentrated in switch I (residues about 30–38) and switch II (about 59–76). Haemoglobin switches
   between the T state (deoxy, like 4HHB) and the R state (oxy). A predictor gives one conformation, usually close
   to whichever state dominates the PDB.
2. **Point mutations.** Predict H-Ras **G12V** or sickle-cell haemoglobin **β E6V**. The models will look almost
   identical to wild type. G12V is oncogenic because it blocks GTP hydrolysis; E6V causes disease because deoxy HbS
   tetramers stick together into long fibres. Predictors are **not** reliable estimators of how mutations change
   function, stability or aggregation.
3. **Ligands and cofactors.** AlphaFold2 has no GDP, no Mg²⁺ and no haem. AlphaFold3 can add them, but most
   post-translational modifications (such as Ras's farnesyl lipid) are still out of reach.
4. **Assemblies and partners.** Haemoglobin works only as an α₂β₂ tetramer; Ras works at the membrane, bound to
   effectors (Raf, PI3K) and regulators (SOS, GAP). Single-chain predictions don't show these.
5. **Disorder is not a failure.** Low-pLDDT regions such as the Ras C-terminal tail may be *functionally*
   disordered.
6. **Dynamics and folding pathways.** Predictors give a final structure, not how the protein moves or how it folds.

> **Question 6.1.** Use AlphaFold Server to predict the haemoglobin tetramer with and without four haem groups.
> Does adding haem change pLDDT or ipTM? Where do the haem groups sit compared with the iron positions in 4HHB?

---

## Exercises (to hand in)

1. **Classical versus learned.** Make a table of Q3 for H-Ras and haemoglobin β with Chou–Fasman, your best neural
   network from `nn` (give its settings), PSIPRED (<http://bioinf.cs.ucl.ac.uk/psipred/>), and the DSSP secondary
   structure of your AlphaFold or ESMFold model (in ChimeraX: `dssp` then `info residues attribute ss_type`).
   Explain the ranking.
2. **Where were they wrong?** Include screenshots from step 4 of the interactive page showing haemoglobin β
   coloured by Chou–Fasman and by your network. Describe where each method fails and why.
3. **3D accuracy.** For your H-Ras and haemoglobin β models, report the fraction of residues within 2 Å, the Cα
   RMSD from `compare` (or step 6) and the mean pLDDT. Is the model worse where pLDDT is low?
4. **Assemblies.** Predict the haemoglobin tetramer with AlphaFold Server, with and without haem. Report ipTM and
   compare the model with 4HHB.
5. **Limits.** In about 150 words, explain why an accurate AlphaFold model of Ras wasn't enough to design sotorasib
   against KRAS-G12C. *(Hint: look up the switch II pocket.)*

---

## Hints for the in-text questions

<details>
<summary>Click to expand</summary>

- **2.1** The buried face of a globin helix is lined with Val, Leu, Phe and Ile, which all have high strand
  propensities (V 1.70, I 1.60, F 1.38, L 1.30). A 7-residue window full of them looks like a strand to
  Chou–Fasman, which can't see the helical spacing of hydrophobic residues every 3–4 positions.
- **2.2** Evolutionary information (MSAs and profiles) and long-range contacts. β-sheets pair strands that can be
  far apart in sequence, which is exactly what a local window can't see.
- **2.3** Training Q3 keeps climbing while test Q3 stalls, below what 300 proteins give. The growing gap is
  **overfitting**: the network starts memorising the 50 training proteins instead of learning general rules.
  Only accuracy on unseen proteins tells you how the method will do on a new protein.
- **2.4** A single residue is a noisy signal: many amino acids occur in strands, helices and loops alike. A column
  that stays hydrophobic across many homologues shows that the *position* is buried in all of them. The strand's
  alternating buried/exposed pattern is also much clearer when averaged over many sequences.
- **3.1** Many errors sit at the ends of helices and strands, where DSSP itself is sensitive to small coordinate
  changes. This is why some benchmarks also report segment-overlap scores (SOV) that forgive small boundary shifts.
- **5.1** NMR structures are ensembles fitted to distance restraints; loops and the switch regions are often
  poorly defined and vary between models. The model may match a crystal structure of Ras better than NMR model 1.
- **5.2** pLDDT and real error usually agree, but pLDDT measures the model's *local* self-consistency. A region
  can be confidently predicted in a different conformation (for example the other switch state of Ras), or placed
  differently relative to the rest of the chain.
- **5.3** No: PAE says the relative placement of the chains is uncertain, even though each chain is confident on
  its own.
- **6.1** Answers vary. Check whether the haem irons land close to the positions in 4HHB (listed as `HETATM FE`
  lines in `data/structures/hb_4hhb.pdb`).

</details>

---

## Further reading

- Anfinsen, C. B. (1973). Principles that govern the folding of protein chains. *Science* 181, 223–230.
- Chou, P. Y. & Fasman, G. D. (1978). Prediction of the secondary structure of proteins from their amino acid
  sequence. *Adv. Enzymol.* 47, 45–148.
- Kyte, J. & Doolittle, R. F. (1982). A simple method for displaying the hydropathic character of a protein.
  *J. Mol. Biol.* 157, 105–132.
- Kabsch, W. & Sander, C. (1983). Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded
  and geometrical features. *Biopolymers* 22, 2577–2637. *(DSSP)*
- Qian, N. & Sejnowski, T. J. (1988). Predicting the secondary structure of globular proteins using neural
  network models. *J. Mol. Biol.* 202, 865–884.
- Rost, B. & Sander, C. (1993). Prediction of protein secondary structure at better than 70% accuracy.
  *J. Mol. Biol.* 232, 584–599. *(PHD)*
- Jones, D. T. (1999). Protein secondary structure prediction based on position-specific scoring matrices.
  *J. Mol. Biol.* 292, 195–202. *(PSIPRED)*
- Klausen, M. S. et al. (2019). NetSurfP-2.0: improved prediction of protein structural features by integrated
  deep learning. *Proteins* 87, 520–527.
- Zhou, J. & Troyanskaya, O. G. (2014). Deep supervised and convolutional generative stochastic network for
  protein secondary structure prediction. *ICML*. *(CB6133 and CB513 data; cleaned version from Drori et al. 2018,
  github.com/idrori/cu-ssp)*
- Fermi, G., Perutz, M. F., Shaanan, B. & Fourme, R. (1984). The crystal structure of human deoxyhaemoglobin at
  1.74 Å resolution. *J. Mol. Biol.* 175, 159–174. *(PDB 4HHB)*
- Kraulis, P. J., Domaille, P. J., Campbell-Burk, S. L., Van Aken, T. & Laue, E. D. (1994). Solution structure and
  dynamics of ras p21·GDP determined by heteronuclear three- and four-dimensional NMR spectroscopy.
  *Biochemistry* 33, 3515. *(PDB 1CRR)*
- Jumper, J. et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature* 596, 583–589.
- Lin, Z. et al. (2023). Evolutionary-scale prediction of atomic-level protein structure with a language model.
  *Science* 379, 1123–1130. *(ESMFold)*
- Mirdita, M. et al. (2022). ColabFold: making protein folding accessible to all. *Nat. Methods* 19, 679–682.
- Abramson, J. et al. (2024). Accurate structure prediction of biomolecular interactions with AlphaFold 3.
  *Nature* 630, 493–500.
