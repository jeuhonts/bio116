# Tutorial: Structure prediction for special protein classes

**Membrane proteins and coiled coils: where generic predictors need your help**

| | |
|---|---|
| **Audience** | Undergraduate bioinformatics (BIO116). Assumes you have done the [Ras and haemoglobin structure prediction tutorial](../protein-structure-prediction/README.md): you know hydropathy plots, pLDDT, PAE and ipTM. |
| **Time** | About 2–3 hours (Parts A1–A3 and B1–B3 offline, about 1.5 h; the AlphaFold and web-server exercises need a browser and about 1 h) |
| **You need** | Python 3.8+ (standard library only), a web browser, and optionally [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) |
| **Files** | `special-classes-tutorial.html` (interactive version: open it in any browser, works offline), `special_classes.py` (companion script) |

By the end you will be able to:

1. Explain why AlphaFold-style predictors need extra help with membrane proteins and coiled coils.
2. Find transmembrane helices with a hydropathy plot, and name the two classic ways this fails.
3. Predict membrane topology with the positive-inside rule.
4. Explain why β-barrel membrane proteins are invisible to hydropathy plots.
5. Place a predicted model in a membrane and check that the result makes sense.
6. Read a heptad repeat, assign a coiled-coil register, and use the rules of thumb for oligomer state.
7. Use AlphaFold-Multimer to test, not assume, how many chains a coiled coil has.

---

## Part 0 — Why these proteins are special

In the H-Ras tutorial a single soluble chain folded well, and pLDDT told you where to trust the model. Two
features of the way predictors work break that simple picture.

| What the predictor does | Why it matters |
|---|---|
| Predicts the chains you give it, **in a vacuum** | There is **no lipid bilayer**. AlphaFold learned what membrane proteins look like from the PDB, but it never places a membrane. You must work out where the membrane is and which side each loop faces. |
| Models **exactly the number of copies you submit** | The oligomeric state is an **input**, not an output. Submit one copy of a coiled coil and you get one helix; submit three and you get a trimer, whether or not that is the real state. |
| Learns from the PDB | Membrane proteins are hard to crystallise and are under-represented, although cryo-EM is closing the gap. Short, idealised coiled coils are over-represented in designed-protein entries. |

**Membrane proteins** are encoded by roughly a quarter of human genes, and many drugs act on them. They come in
two architectures: **α-helical bundles** (most plasma-membrane and inner-membrane proteins) and **β-barrels** (the
outer membranes of Gram-negative bacteria, mitochondria and chloroplasts).

**Coiled coils** are two to seven α-helices wound around each other. The sequence pattern is easy to detect, but two
things are hard: the **register** (which residues sit in the core) and the **oligomer state** (how many helices, and
whether they are parallel or antiparallel).

The theme of this tutorial: a predictor answers the question you ask it. For these classes you have to supply part of
the answer yourself.

### The example proteins

The script and the HTML page contain these proteins as presets. List them with:

```bash
cd tutorials/special-protein-classes
python3 special_classes.py list
```

| Preset | Protein | UniProt | Length | Why it's here |
|---|---|---|---|---|
| `bR` | Bacteriorhodopsin (precursor) | P02945 | 262 | Classic 7-helix membrane protein |
| `gpa` | Glycophorin A (precursor) | P02724 | 150 | Single TM helix plus a signal peptide |
| `rho` | Rhodopsin, bovine | P02699 | 348 | G protein-coupled receptor (GPCR), 7 helices |
| `lacy` | Lactose permease LacY | P02920 | 417 | 12 helices, strong positive-inside signal |
| `aqp1` | Aquaporin-1 | P29972 | 269 | 6 helices plus 2 half-helices; a hard case |
| `b2ar` | β2-adrenergic receptor | P07550 | 413 | GPCR, 7 helices |
| `ompa` | Outer membrane protein A (precursor) | P0A910 | 346 | 8-stranded β-barrel |
| `tpm1` | Tropomyosin α-1 | P09493 | 284 | Continuous two-stranded coiled coil |
| `hras` | H-Ras | P01112 | 189 | Soluble negative control, from the H-Ras tutorial |
| `gcn4-p1`, `-pIL`, `-pII`, `-pLI` | GCN4 leucine zipper and core variants | P03069 (249–281) | 33 | Dimer / trimer / tetramer switch |

Where each sequence and reference topology comes from is listed at the end of this README. None of the membrane
protein sequences was typed by hand.

---

## Part A1 — Finding transmembrane helices

The hydrocarbon core of a lipid bilayer is about **30 Å** thick. An α-helix rises 1.5 Å per residue, so about
**20 hydrophobic residues** span it. Kyte and Doolittle (1982) found that averaging hydropathy over a
**19-residue window** and calling windows above **1.6** picks out many transmembrane (TM) helices.

The script calls helices one at a time: it takes the highest window above the cutoff, marks those 19 residues as a
helix, then takes the next-highest window that doesn't overlap a helix it has already called.

```bash
python3 special_classes.py tm --preset bR
```

```
== bR: Bacteriorhodopsin, Halobacterium salinarum (precursor; 7 TM helices) ==
Length 262 aa
Kyte-Doolittle, 19-residue window: max +2.12 at residue 228; cutoff 1.6

Predicted TM helices (5):
   #   start   end   peak KD
   1      24    42   +1.84
   2      57    75   +1.82
   3     121   139   +1.81
   4     148   166   +1.88
   5     219   237   +2.12

Positive-inside rule: K + R within 15 residues of each helix end
  loop        N    1    2    3    4    C
  side        A    B    A    B    A    B
  K+R         1    3    0    2    2    2
  side A (N-terminal side) 3   side B 7   ->  N-terminus predicted OUT, C-terminus in

Reference: set160 (TMHMM benchmark, experimental topology; github.com/tjs23/python_ml_course set160.labels)
  7 TM helices: 23-42, 57-76, 95-114, 121-140, 148-167, 191-210, 217-236
  N-terminus out
  prediction vs reference: 5 matched, 0 extra, 2 missed (match = overlap of 5 or more residues)

Topology track: i inside, o outside, M membrane (reference letters as in its source)
    1 seq  MLELLPTAVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITT
      pred oooooooooooooooooooooooMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiMMMM
      ref  ooooooooooooooooooooooMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiMMMM

   61 seq  LVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGT
      pred MMMMMMMMMMMMMMMooooooooooooooooooooooooooooooooooooooooooooo
      ref  MMMMMMMMMMMMMMMMooooooooooooooooooMMMMMMMMMMMMMMMMMMMMiiiiii
  ...
```

Five of seven helices are found, and they line up well with the reference. Helices C (95–114) and F (191–210) are
missed.

> **Question A1.1.** Look at the sequence of helix C (`RYADWLFTTPLLLLDLALLV`, residues 95–114 of the
> precursor). Which residues pull its average down? Numbering note: the mature protein starts at residue 14 of
> the precursor, so precursor Asp98 and Asp109 are the well-known **Asp85** and **Asp96** of the proton pump.

Hydropathy analysis has two classic failure modes.

**1. Signal peptides look like TM helices.** A signal peptide has a hydrophobic core of 7–15 residues that sends the
protein into the membrane. It's then cut off, but in the sequence it looks like a short TM helix.

```bash
python3 special_classes.py tm --preset gpa
```

```
Predicted TM helices (2):
   #   start   end   peak KD
   1       5    23   +2.04
   2      92   110   +2.67
  ! helix 1 starts at residue 5: it could be a cleaved signal peptide, not a TM helix. Hydropathy cannot tell; use Phobius, SignalP 6.0 or DeepTMHMM.

Positive-inside rule: K + R within 15 residues of each helix end
  loop        N    1    C
  side        A    B    A
  K+R         1    2    4
  side A (N-terminal side) 5   side B 2   ->  N-terminus predicted IN, C-terminus in

Reference: set160 (TMHMM benchmark, experimental topology; github.com/tjs23/python_ml_course set160.labels)
  1 TM helix: 92-114
  N-terminus out
  prediction vs reference: 1 matched, 1 extra, 0 missed (match = overlap of 5 or more residues)
```

Glycophorin A has one TM helix. The extra "helix" at 5–23 is its signal peptide. Worse, the fake helix flips the
topology: the N-terminus is now predicted inside, which is wrong. This is exactly the problem **Phobius** (Käll et
al., 2004) was built to solve: it models signal peptides and TM helices together.

The warning is a crude rule (any first helix starting in the first 20 residues). It also fires for LacY (helix 1 at
14–32) and aquaporin-1 (13–31), where the first helix is real. Hydropathy alone can't separate the two cases.

**2. Amphipathic helices.** A helix with one hydrophobic face and one polar face (for example a helix that lies
flat on the membrane surface) can have a high average hydropathy without crossing the membrane. And, as helix C of
bacteriorhodopsin shows, a real TM helix with functional charged residues can have a low average.

> **Question A1.2.** Use the interactive page (step 2) to change the window and cutoff for bacteriorhodopsin. Can
> you find a setting that finds all seven helices without also calling extra ones? What does that tell you about
> tuning a cutoff on one protein?

### Modern tools

| Tool | Method | Notes |
|---|---|---|
| **TMHMM 2.0** (Krogh et al., 2001) | Hidden Markov model with states for helix core, caps and loops | Fast and classic. Has no signal peptide model, so it can make the glycophorin A mistake. |
| **Phobius** (Käll et al., 2004) | HMM that combines signal peptides and TM helices | Use it whenever the N-terminus is hydrophobic. <https://phobius.sbc.su.se/> |
| **TOPCONS** | Consensus of several predictors | <https://topcons.cbr.su.se/> |
| **DeepTMHMM** (Hallgren et al., 2022) | Deep learning on protein language-model embeddings | Predicts α-helical *and* β-barrel topologies and signal peptides. <https://dtu.biolib.com/DeepTMHMM> |

---

## Part A2 — Topology: the positive-inside rule

Finding the helices is only half the job. The **topology** says which loops are inside the cell and which are outside.
Von Heijne (1992) showed that loops on the cytoplasmic side are rich in **lysine and arginine**: the
**positive-inside rule**. Because loops alternate sides, you can add up K + R on "side A" (the N-terminal side and
every second loop after it) and on "side B". The side with more positive charge is probably inside.

The script counts K + R within 15 residues of each helix end (for the termini, only the end next to the membrane).
Use `--use-reference` to apply the rule to the reference helices, which tests the rule separately from helix finding:

```bash
python3 special_classes.py tm --preset lacy --use-reference
```

```
Positive-inside rule: K + R within 15 residues of each helix end
  loop        N    1    2    3    4    5    6    7    8    9   10   11    C
  side        A    B    A    B    A    B    A    B    A    B    A    B    A
  K+R         1    1    3    0    5    0    2    1    2    0    2    0    2
  side A (N-terminal side) 17   side B 2   ->  N-terminus predicted IN, C-terminus in
```

LacY is a textbook case: 17 against 2. Both reference topologies (set160 and the DeepTMHMM set) agree that both ends
are in the cytoplasm.

Results for all the presets with reference topologies (the positive-inside step, real output of `tm` and
`tm --use-reference`):

| Protein | Helices found (19 / 1.6) | Reference helices | K+R side A / B, predicted helices | Prediction | K+R side A / B, reference helices | Prediction | Reference N-terminus |
|---|---|---|---|---|---|---|---|
| Bacteriorhodopsin | 5 | 7 | 3 / 7 | N-out | 3 / 9 | N-out | out |
| Glycophorin A | 2 (one is the signal peptide) | 1 | 5 / 2 | **N-in (wrong)** | 1 / 5 | N-out | out |
| Rhodopsin | 6 | 7 | 1 / 10 | N-out, **C-out (wrong)** | 1 / 12 | N-out, C-in | out |
| LacY | 10 | 12 | 16 / 1 | N-in | 17 / 2 | N-in | in |

And for the two presets without a reference topology on the page:

```bash
python3 special_classes.py tm --preset aqp1
python3 special_classes.py tm --preset b2ar
```

| Protein | Helices found | Prediction | What is known |
|---|---|---|---|
| Aquaporin-1 | 5 (13–31, 94–112, 138–156, 166–184, 214–232) | K+R 6 / 9, N-out, C-in | 6 TM helices, both termini in the cytoplasm (see UniProt P29972, *Topology*) |
| β2-adrenergic receptor | 6 (31–49, 71–89, 111–129, 152–170, 200–218, 277–295) | K+R 5 / 15, N-out, C-out | GPCR: 7 TM helices, N-terminus outside, C-terminus inside |

> **Question A2.1.** When the rule fails above, is it the rule itself or the helix finding that went wrong? What
> single mistake flips every loop after it?

> **Question A2.2.** Aquaporin-1 contains two short "half-helices" that dip into the membrane from opposite sides
> and meet in the middle, each carrying an NPA motif (Asn-Pro-Ala at residues 76–78 and 192–194). Why would a model
> that only knows "in, out or across" struggle with them?

---

## Part A3 — β-barrels hide from hydropathy plots

In a β-barrel each membrane-crossing strand is only 7–10 residues long, and its side chains point alternately into
the lipid and into the barrel's interior. So the sequence **alternates** hydrophobic and polar residues, and a
19-residue average stays low.

Instead, measure the **hydrophobic moment at 180°** (Eisenberg et al., 1982): give each residue's hydropathy a
direction 180° from its neighbour's, add the vectors over a short window, and divide by the window length. Strict
alternation gives a large moment. (At 100° the same calculation detects amphipathic α-helices.)

```bash
python3 special_classes.py barrel --preset ompa
```

```
== ompa: Outer membrane protein A, E. coli (precursor; 8-strand beta-barrel + periplasmic domain) ==
TM helices called by hydropathy (19-window > 1.6): 1 at 4-22
Max 19-window hydropathy: +1.84
Hydrophobic moment at 180 deg, 9-residue window (high = hydrophobic and polar residues alternate)
  residues with moment >= 1.2: 90 of 338
  mean moment in reference strands  1.31
  mean moment everywhere else       0.79
  reference strands with 3 or more residues above 1.2: 5 of 8
  mean 19-window hydropathy over strand residues -0.32 (far below the 1.6 helix cutoff)
  '#' moment >= 1.2   '-' below   '.' window does not fit

    1 seq  MKKTAIAIAVALAGFATVAQAAPKDNTWYTGAKLGWSQYHDTGFINNNGPTHENQLGAGA
      mu   ....----------------------------#---------------------------
      ref  SSSSSSSSSSSSSSSSSSSSSPPPPPPBBBBBBBBOOOOOOOOOOOOOOOOOOOOOOBBB

   61 seq  FGGYQVNPYVGFEMGYDWLGRMPYKGSVENGAYKAQGVQLTAKLGYPITDDLDIYTRLGG
      mu   -----#########---------##---#---------#####################-
      ref  BBBBBPPPPBBBBBBBOOOOOOOOOOOOOOOOOOOOOOBBBBBBBBPPPPPPBBBBBBBB
  ...
```

The reference letters come from the DeepTMHMM training set: `S` signal peptide, `B` membrane β-strand, `O` outside,
`P` periplasm. The only "TM helix" is the signal peptide again. The strands average −0.32 on the 19-residue scale, but
their mean moment (1.31) is clearly higher than the rest of the protein (0.79).

> **Question A3.1.** The moment is not a clean detector: only 5 of 8 strands have three or more residues above 1.2,
> and parts of the soluble periplasmic domain (residues 191–346) score high too. Run
> `python3 special_classes.py barrel --preset hras`. Why would a soluble protein like H-Ras also show high-moment
> stretches?

Dedicated β-barrel predictors (DeepTMHMM, and earlier tools such as PRED-TMBB and BOMP) combine this kind of
signal with evolutionary information and with the fact that outer-membrane proteins are made with a signal peptide.

---

## Part A4 — From sequence to a model in a membrane

AlphaFold models of membrane proteins are often good for the fold, because the PDB contains many membrane protein
structures to learn from. What the model can't give you is the membrane. A workflow:

1. **Topology first.** Run DeepTMHMM (or Phobius if the N-terminus is hydrophobic). Compare with the UniProt
   *Topology* section, which is often based on experiments.
2. **Get a model.** Download it from the AlphaFold DB (<https://alphafold.ebi.ac.uk>), or run ColabFold or AlphaFold
   Server. Check pLDDT as usual.
3. **Place the membrane.** Upload the model to the **PPM server** (<https://opm.phar.umich.edu/ppm_server>). It finds
   the orientation and hydrophobic thickness that minimise the energy of transferring the protein from water into
   the membrane (Lomize et al., 2012). For experimental structures the **OPM database**
   (<https://opm.phar.umich.edu>) already has the answer, and **MemProtMD** (<https://memprotmd.bioch.ox.ac.uk/>)
   shows PDB structures self-assembled into a bilayer by coarse-grained simulation.
4. **Check it.**
   - The hydrophobic belt should be about **30 Å** thick.
   - Trp and Tyr should cluster at the two interfaces ("aromatic belts").
   - K and R should be mostly on the cytoplasmic side (Part A2).
   - The predicted TM segments should cross the slab and the loops should stay outside it.
5. **Multimers.** Many membrane proteins are oligomers (aquaporins are tetramers). A monomer model leaves a
   hydrophobic protein–protein interface exposed to the lipid. Submit the right number of copies.

In ChimeraX you can open the PPM output directly: it contains dummy atoms (residue name `DUM`) that mark the two
membrane planes.

> **Question A4.1.** AlphaFold 3 can place some lipids and detergents as ligands. Does that mean it models the
> membrane? What would you still check?

---

## Part B1 — The heptad repeat

An α-helix has 3.6 residues per turn (100° per residue). When helices wind around each other in a coiled coil, the
repeat becomes **3.5 residues per turn**: seven residues in two turns, the **heptad**. The positions are labelled
**a b c d e f g**:

| Position | Where it points | Typical residues |
|---|---|---|
| **a**, **d** | Into the core, packing against the other helices | Leu, Ile, Val, Met (sometimes Asn at *a*) |
| **e**, **g** | The edges of the core | Glu, Lys, Arg: salt bridges between *g* of one helix and *e′* of the next |
| **b**, **c**, **f** | Out into the solvent | Polar and charged residues, Ala |

The classic example is the **GCN4 leucine zipper**, GCN4-p1 (PDB 2ZTA; O'Shea et al., 1991). Its register, as found
by the script's heptad scan and shown on a text helical wheel:

```bash
python3 special_classes.py wheel --preset gcn4-p1 --end 33
```

```
== gcn4-p1: GCN4-p1 leucine zipper (wild type, dimer, PDB 2ZTA): residues 1-33, residue 1 = 'g' (best register from `coils`) ==

Helical wheel (100 degrees per residue, 3.6 residues per turn).
Spokes in the order they appear around the wheel; innermost residue first.

  a   0 deg  core     M2 V9 N16 V23 V30
  e  40 deg  edge     E6 L13 E20 K27
  b 100 deg  surface  K3 E10 Y17 A24 G31
  f 140 deg  surface  D7 S14 N21 K28
  c 200 deg  surface  Q4 E11 H18 R25 E32
  g 240 deg  edge     R1 K8 K15 E22 L29
  d 300 deg  core     L5 L12 L19 L26 R33

Heptad table
   a  b  c  d  e  f  g
                     R
   M  K  Q  L  E  D  K
   V  E  E  L  L  S  K
   N  Y  H  L  E  N  E
   V  A  R  L  K  K  L
   V  G  E  R

Hydrophobic (L I V M F Y) at a/d: 8 of 10
Charged (E K R D) at e/g: 7 of 9
```

The leucines at every *d* position (5, 12, 19, 26) gave the "leucine zipper" its name. The *a* positions are
Met2, Val9, **Asn16**, Val23 and Val30. The buried Asn16 is a key determinant of the parallel dimer.

Force a wrong register to see the difference (`--register c` makes residue 1 a *c* position):

```bash
python3 special_classes.py wheel --preset gcn4-p1 --end 33 --register c
```

```
  a   0 deg  core     E6 L13 E20 K27
  ...
  d 300 deg  core     M2 V9 N16 V23 V30
```

> **Question B1.1.** In the wrong register, what ends up in the "core"? Try the interactive wheel (step 6 of the
> HTML page) and find the register by eye before pressing "Show the best register".

---

## Part B2 — Scanning for coiled coils

**COILS** (Lupas et al., 1991) slides a 28-residue window (four heptads) along a sequence, scores all seven
registers with a position-specific scoring matrix built from known coiled coils, and keeps the best score for each
residue.

The script uses the same idea with a **simplified, transparent scoring table that we made up for teaching**. It is
**not** the published COILS (MTK/MTIDK) matrices:

| Position | Score |
|---|---|
| *a*, *d* | 2 for L, I, V, M, F, Y; otherwise 0 |
| *e*, *g* | 1 for E, K, R, D; 0.5 for Q, N, S, T, H, A; otherwise 0 |
| *b*, *c*, *f* | 1 for E, K, R, D, Q, N, S, T, H, A, G; otherwise 0 |
| Pro, anywhere | −1 (proline breaks helices) |

A perfect heptad scores 9; each window score is divided by the maximum, giving 0–1. A residue is called coiled coil
when its best score is at least **0.75** and it is part of a run of 14 or more such residues. You can read the whole
method in `heptad_score()` and `coils_scan()` in the script.

```bash
python3 special_classes.py coils --preset tpm1
```

```
== tpm1: Tropomyosin alpha-1, human (continuous two-stranded coiled coil) ==
Simplified heptad scoring, 28-residue window (teaching method, not COILS)
Score range 0.65 to 0.88; 220 of 284 residues >= 0.75
Coiled-coil segments (14 or more residues >= threshold): 30-119, 155-284

score: first decimal of the best score (8 = 0.80-0.89)   reg: heptad letter where score >= 0.75
    1 seq  MDAIKKKMQMLKLDKENALDRAEQAEADKKAAEDRSKQLEDELVSLQKKLKGTEDELDKY
      scr  666666666666666666666677777777778888888888888888888888888888
      reg  .............................bcdefgabcdefgabcdefgabcdefgabcd
  ...
```

Summary over the presets (real output of `coils` and `analyze`):

| Preset | Max score | Coiled-coil segments (≥ 0.75) |
|---|---|---|
| GCN4-p1 | 0.89 | 1–33 |
| GCN4-pIL, -pII, -pLI | 0.94 | 1–33 |
| Tropomyosin | 0.88 | 30–119, 155–284 |
| H-Ras | 0.72 | none |
| LacY | 0.74 | none |
| β2-adrenergic receptor | 0.71 | none |
| Bacteriorhodopsin | 0.67 | none |
| OmpA | 0.65 | none |

The toy method separates these examples, but only just: H-Ras (0.72) and LacY (0.74) come close to the threshold.
It also misses residues 1–29 and 120–154 of tropomyosin, which is a coiled coil along its whole length but has
alanines and other small residues in parts of its core.

> **Question B2.1.** Transmembrane helices are hydrophobic at every position, including *b*, *c* and *f*. How does
> the scoring table keep them from scoring as coiled coils? Check with `coils --preset bR`.

### Modern tools

| Tool | Method |
|---|---|
| **Marcoil** | HMM; in the MPI Bioinformatics Toolkit (<https://toolkit.tuebingen.mpg.de/tools/marcoil>) |
| **DeepCoil** (Ludwiczak et al., 2019) | Deep learning; DeepCoil2 is in the MPI Toolkit (<https://toolkit.tuebingen.mpg.de/tools/deepcoil>) |
| **CoCoNat** | Predicts segments, register and oligomer state from protein language-model embeddings (<https://coconat.biocomp.unibo.it>) |
| **Waggawagga** | Coiled coils and single α-helices, with helical wheels (<https://waggawagga.motorprotein.de>) |

---

## Part B3 — Dimer, trimer or tetramer?

Harbury, Zhang, Kim and Alber (1993) changed only the core of GCN4-p1. They replaced the four *a* residues at
positions 9, 16, 23 and 30 and the four *d* residues at 5, 12, 19 and 26 (Met2 and Arg33 were not changed), and named
each peptide after its *a* and *d* residues. The script builds the variants by exactly this substitution:

| Peptide | *a* | *d* | Sequence | State (Harbury et al., 1993) | PDB |
|---|---|---|---|---|---|
| GCN4-p1 | Val (+ Asn16) | Leu | `RMKQLEDKVEELLSKNYHLENEVARLKKLVGER` | dimer | 2ZTA |
| pIL | Ile | Leu | `RMKQLEDKIEELLSKIYHLENEIARLKKLIGER` | dimer | — |
| pII | Ile | Ile | `RMKQIEDKIEEILSKIYHIENEIARIKKLIGER` | trimer | 1GCM |
| pLI | Leu | Ile | `RMKQIEDKLEEILSKLYHIENELARIKKLLGER` | tetramer | 1GCL |

Changing eight core residues switched the oligomer state, because β-branched side chains (Ile, Val) and Leu pack
differently at *a* and *d*. The script turns this into rules of thumb:

```bash
python3 special_classes.py coils --preset gcn4-pLI
```

```
a positions: MLLLL
d positions: IIIIR
Oligomer rule of thumb: TETRAMER-like (Leu a, Ile d: like GCN4-pLI)
  (Only for short, GCN4-like zippers. Real oligomer state also depends on e/g residues, length and context: test 2, 3 and 4 copies in AlphaFold-Multimer.)
```

For the other peptides it prints `DIMER-like ...; Asn at an a position favours a dimer` (p1), `DIMER-like` (pIL) and
`TRIMER-like` (pII). For tropomyosin it prints `no clear rule applies`: natural coiled coils are more varied than the
Harbury peptides.

### Testing oligomer state with AlphaFold

AlphaFold never refuses the number of copies you give it; it builds the best bundle it can for each. So you have to
**compare** runs:

1. In AlphaFold Server (<https://alphafoldserver.com>) or ColabFold (join copies with `:` in the sequence box),
   predict GCN4-pLI with **2, 3 and 4 copies**.
2. For each run record **ipTM** and **pTM**, and look at the PAE between chains.
3. Repeat for GCN4-p1 and pII.

> **Question B3.1.** Does the run with the highest ipTM match the experimental state for all three peptides? If
> two runs are almost equally confident, what experiment would you use to decide? *(Hint: size-exclusion
> chromatography with multi-angle light scattering, or analytical ultracentrifugation.)*

### Coiled coils can be built, not only predicted

Crick (1953) showed that a coiled coil can be described with a few parameters: the supercoil radius, the pitch
and the phase of each helix. Because the geometry is so regular, you can build coiled coils from these equations.
**CCBuilder 2.0** (Wood & Woolfson, 2018; <http://coiledcoils.chm.bris.ac.uk/ccbuilder2>) and the Python library
**ISAMBARD** (<https://github.com/woolfson-group/isambard>) do this, and they are widely used to design new coiled coils.

> **Question B3.2.** Build a GCN4-p1 dimer in CCBuilder 2.0 and compare it with 2ZTA and with your AlphaFold
> dimer (Cα RMSD in ChimeraX with `matchmaker`). Which model is closer to the crystal structure?

---

## Part C — Why it matters (briefly)

- **Drug targets.** GPCRs such as the β2-adrenergic receptor, and ion channels and transporters, are membrane
  proteins. To model a drug-binding pocket you need the right topology and the right membrane position.
- **Disease variants.** A mutation that adds a charge in the middle of a TM helix, or that changes the K/R balance
  of a loop, can stop the protein from inserting correctly.
- **Assembly.** Coiled coils hold together transcription factors (GCN4, Fos/Jun), muscle proteins (tropomyosin,
  myosin), membrane-fusion machinery (SNAREs) and viral fusion proteins. The oligomer state is central to how they
  work.

---

## Exercises (to hand in)

1. **Hydropathy versus deep learning.** Run DeepTMHMM on bacteriorhodopsin, glycophorin A, aquaporin-1 and OmpA
   (`python3 special_classes.py fasta --preset bR`, and so on). Make a table of the number of TM segments from the
   19 / 1.6 scan, from DeepTMHMM and from the reference. Explain each disagreement in one sentence.
2. **Topology.** Take the six TM helices of aquaporin-1 from UniProt P29972 (*Topology* section), count K + R by hand
   within 15 residues of each helix end, and apply the positive-inside rule. Does it now give the
   right orientation? Where do the two NPA half-helices end up?
3. **Membrane placement.** Download the AlphaFold DB model of aquaporin-1 and run it through the PPM server. Report
   the hydrophobic thickness and the first and last residue of each TM segment that lies inside the slab. Include a
   figure with the membrane planes shown.
4. **Coiled coils.** Choose a protein you're interested in (for example a SNARE or a myosin tail from UniProt).
   Compare the segments from `coils --seq ...` with DeepCoil or CoCoNat.
5. **Oligomer state.** Do the 2 / 3 / 4-copy AlphaFold experiment from Part B3 for GCN4-p1, pII and pLI. Tabulate
   ipTM and discuss whether it recovers the Harbury results.

---

## Hints for the in-text questions

<details>
<summary>Click to expand</summary>

- **A1.1** Arg95, Asp98 and Asp109 are charged (−4.5, −3.5 and −3.5 on the Kyte–Doolittle scale); Tyr96 (−1.3),
  Trp99 (−0.9) and Pro104 (−1.6) also pull the average down. In mature numbering these are Arg82, Asp85 and Asp96,
  residues that move protons across the membrane, so the helix is "less hydrophobic" for a functional reason.
- **A1.2** Lowering the cutoff finds more of the real helices, but at some point loops or signal peptides start
  being called too. A cutoff tuned on one protein doesn't transfer to others.
- **A2.1** Every error above comes from the helix list. With the reference helices the rule is right for all four
  proteins. A single missing or extra helix shifts the alternation of every loop after it (glycophorin A's signal
  peptide; the missing helix 7 of rhodopsin; the missing helices of aquaporin-1).
- **A2.2** A half-helix enters and leaves the membrane on the same side, so it doesn't flip the side of the next
  loop. A model with only "in / across / out" states has to call it either a full TM helix (which flips the topology)
  or a loop.
- **A3.1** β-strands in soluble proteins also alternate: one face packs into the hydrophobic core and the other faces
  the solvent. H-Ras has a six-stranded β-sheet.
- **A4.1** No. Ligand lipids are placed where they bind, not as a continuous bilayer. You still need to find the
  membrane plane and thickness (PPM) and check the aromatic belts and K/R distribution.
- **B1.1** In the *c* register the core positions hold E6, L13, E20, K27 at *a* and M2, V9, N16, V23, V30 at *d*:
  the leucines at 5, 12, 19, 26 move to *g*. Only half of the core is hydrophobic.
- **B2.1** A TM helix has hydrophobic residues at *b*, *c* and *f* too, which score 0 there. Bacteriorhodopsin
  never gets above 0.67.
- **B3.1** Record what you get; don't assume the answer. If the ipTMs are similar, the model isn't telling you the
  oligomer state.
- **B3.2** Report your own RMSD values. Expect both models to be close to 2ZTA; a parametric model is perfectly
  regular, so it can't reproduce local irregularities of the real structure.

</details>

---

## Where the sequences come from

UniProt, RCSB and EBI were not reachable when this tutorial was written, so every sequence was downloaded from a
public GitHub repository and checked against at least one other, independent copy. Only the 33-residue GCN4-p1
peptide was typed by hand, and it matches the copy listed below. `python3 special_classes.py list` prints the same
information.

| Preset | Sequence source (repository and path) | Cross-checked against |
|---|---|---|
| `bR` | `tjs23/python_ml_course` `set160.labels`, entry `BACR_HALHA` (the old Swiss-Prot name of P02945) | UniProt FASTA of P02945 in `jomimc/AF2_Stability_PRL_2024` `fasta/P02945.fasta` (identical); mature chain from residue 14 in `NENUBioCompute/DMCTOP` `Dataset.fasta` |
| `gpa` | `set160.labels`, entry `GLPA_HUMAN` | Residues 20–150 (mature chain) identical in `DMCTOP` `Dataset.fasta` |
| `rho` | `set160.labels`, entry `OPSD_BOVIN` | Identical in `DMCTOP` `Dataset.fasta` (P02699) |
| `lacy` | `set160.labels`, entry `LACY_ECOLI` | Identical to the P02920 entry of `andregodinhodtu/TMHMM` `data/raw_data/DeepTMHMM.3line`. (The `DMCTOP` copy has Gly154, the C154G mutant used for crystallisation.) |
| `aqp1` | `SBRG/Recon3D` `.../358__46__1_protein/sequences/P29972-1.fasta` | `korcsmarosgroup/HMIpipeline` `.../protein_sequences/P29972.fasta` and `DMCTOP` `Dataset.fasta` (all identical) |
| `b2ar` | `pablogainza/gprots` `uniref/ADRB2/P07550.fasta` | `andrewcboardman/pyGEMME` `examples/ADRB2_HUMAN/results/ADRB2_HUMAN.fasta` and ProteinGym `ADRB2_HUMAN_Jones_2020.fasta` in `ai4protein/Pro-Prime` (all identical) |
| `ompa` | `andregodinhodtu/TMHMM` `data/raw_data/DeepTMHMM.3line`, entry `P0A910` | UniProt FASTA in `javieriserte/hack-a-ton-dp` `data/fastas/P0A910.fasta`, `AlessioDelConte/pytorch_disprot` and `YaoYinYing/FoldDock` (all identical) |
| `tpm1` | `Thivby/Master_dissertation` `Fastas/pet_own_method/lmcd1_pet_P09493.fasta` | `dmx2/myocarditis` `data/myocarditis_antigens.fasta`, `pixelatedbus/klasifikasi-protein` `sitoplasma.fasta` (identical) |
| `hras` | Translated from the four H-Ras coding exons in `../../triplets.txt` (nt 6229–6339, 6607–6785, 6939–7098, 7796–7915) | UniProt P01112 (see the H-Ras tutorial) |
| `gcn4-p1` | Typed (well-established peptide) | `neeleshsoni21/COCONUT` `coconut/example/pipeline1/2zta/2zta.fasta` (identical) |
| `gcn4-pII` | Built from p1 by the substitution rule in Part B3 | `Indicator/RaptorX-SS8` `examples/1gcm.seq` (identical) |
| `gcn4-pLI` | Built from p1 by the substitution rule | Residues 1–31 identical to SCOP/ASTRAL 2.06 domain `d1gcla_` (PDB 1GCL), in `wiwie/clustevalDockerRepository` `.../astral_class_h/...` |
| `gcn4-pIL` | Built from p1 by the substitution rule | Not independently checked |

Sanity checks on the sequences: bacteriorhodopsin is 262 residues with Lys229 (= Lys216 of the mature protein,
which carries the retinal) and Asp98/Asp109 (= Asp85/Asp96); glycophorin A is 150 residues; bovine rhodopsin 348;
LacY 417; aquaporin-1 269, with NPA motifs at 76–78 and 192–194; the β2-adrenergic receptor 413, with the DRY motif
at 130–132; OmpA 346; tropomyosin 284.

**Reference topologies.**

- **set160** (bacteriorhodopsin, glycophorin A, rhodopsin, LacY) is the 160-protein set with experimentally
  determined topologies used to develop and test TMHMM (Sonnhammer et al., 1998; Krogh et al., 2001), as
  redistributed in `tjs23/python_ml_course`. A second copy, `simonrozsival/mff-transmembrane-hmm`
  `set160.labels.txt`, has identical sequences. Its TM boundaries come from the older topology literature and differ
  by a few residues from modern structure-based boundaries: compare the two LacY references.
- **DeepTMHMM training set** (OmpA, LacY) has topologies derived from 3D structures (Hallgren et al., 2022). The
  copy used here is `DeepTMHMM.3line` in the student project `andregodinhodtu/TMHMM`, not the official DTU
  download, which was not reachable.

No reference topology is shipped for aquaporin-1 or the β2-adrenergic receptor: look them up in UniProt (P29972,
P07550) or OPM.

---

## Further reading

- Kyte, J. & Doolittle, R. F. (1982). A simple method for displaying the hydropathic character of a protein.
  *J. Mol. Biol.* 157, 105–132.
- Eisenberg, D., Weiss, R. M. & Terwilliger, T. C. (1982). The helical hydrophobic moment: a measure of the
  amphiphilicity of a helix. *Nature* 299, 371–374.
- von Heijne, G. (1992). Membrane protein structure prediction: hydrophobicity analysis and the positive-inside
  rule. *J. Mol. Biol.* 225, 487–494.
- Sonnhammer, E. L. L., von Heijne, G. & Krogh, A. (1998). A hidden Markov model for predicting transmembrane
  helices in protein sequences. *Proc. Int. Conf. Intell. Syst. Mol. Biol.* 6, 175–182.
- Krogh, A., Larsson, B., von Heijne, G. & Sonnhammer, E. L. L. (2001). Predicting transmembrane protein topology
  with a hidden Markov model: application to complete genomes. *J. Mol. Biol.* 305, 567–580. *(TMHMM)*
- Käll, L., Krogh, A. & Sonnhammer, E. L. L. (2004). A combined transmembrane topology and signal peptide
  prediction method. *J. Mol. Biol.* 338, 1027–1036. *(Phobius)*
- Hallgren, J. et al. (2022). DeepTMHMM predicts alpha and beta transmembrane proteins using deep neural networks.
  *bioRxiv* 2022.04.08.487609.
- Lomize, M. A., Pogozheva, I. D., Joo, H., Mosberg, H. I. & Lomize, A. L. (2012). OPM database and PPM web server:
  resources for positioning of proteins in membranes. *Nucleic Acids Res.* 40, D370–D376.
- Crick, F. H. C. (1953). The packing of α-helices: simple coiled-coils. *Acta Crystallogr.* 6, 689–697.
- Lupas, A., Van Dyke, M. & Stock, J. (1991). Predicting coiled coils from protein sequences. *Science* 252,
  1162–1164. *(COILS)*
- O'Shea, E. K., Klemm, J. D., Kim, P. S. & Alber, T. (1991). X-ray structure of the GCN4 leucine zipper, a
  two-stranded, parallel coiled coil. *Science* 254, 539–544. *(PDB 2ZTA)*
- Harbury, P. B., Zhang, T., Kim, P. S. & Alber, T. (1993). A switch between two-, three-, and four-stranded
  coiled coils in GCN4 leucine zipper mutants. *Science* 262, 1401–1407.
- Wood, C. W. & Woolfson, D. N. (2018). CCBuilder 2.0: powerful and accessible coiled-coil modeling. *Protein Sci.*
  27, 103–111.
- Ludwiczak, J., Winski, A., Szczepaniak, K., Alva, V. & Dunin-Horkawicz, S. (2019). DeepCoil—a fast and accurate
  prediction of coiled-coil domains in protein sequences. *Bioinformatics* 35, 2790–2795.
- Jumper, J. et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature* 596, 583–589.
- Evans, R. et al. (2021). Protein complex prediction with AlphaFold-Multimer. *bioRxiv* 2021.10.04.463034.
