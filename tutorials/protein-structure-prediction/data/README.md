# Secondary-structure datasets

| File | Proteins | Residues | Use |
|---|---|---|---|
| `cb6133filtered.tsv` | 5,534 | 1,183,318 | Training |
| `cb513.tsv` | 514 | 84,765 | Test (the standard CB513 benchmark) |

Each line holds a protein sequence and its DSSP secondary structure in 8 states, separated by a tab. The DSSP
states are H (α-helix), G (3₁₀ helix), I (π helix), E (strand), B (bridge), T (turn), S (bend) and L (coil).
The tutorial reduces them to 3 states in the usual way: H/G/I → H, E/B → E, everything else → C.

The sets come from Zhou & Troyanskaya (2014), *Deep supervised and convolutional generative stochastic network
for protein secondary structure prediction*, ICML. The cleaned versions used here (duplicates removed, training
and test sets disjoint) come from Drori et al. (2018), <https://github.com/idrori/cu-ssp>, `model_3/`, and are
redistributed under that repository's MIT License:

> MIT License. Copyright (c) 2018 Iddo Drori. Permission is hereby granted, free of charge, to any person
> obtaining a copy of this software and associated documentation files (the "Software"), to deal in the
> Software without restriction, including without limitation the rights to use, copy, modify, merge, publish,
> distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
> furnished to do so, subject to the following conditions: The above copyright notice and this permission
> notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED
> "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED.

## Which training proteins are excluded

`excluded_training.txt` lists 17 training proteins that the `nn` command and the interactive page leave out,
so that the predictions for the example proteins are a fair test: every protein with at least 25 % identity over at
least half of H-Ras, haemoglobin α or haemoglobin β (including a haemoglobin β variant), plus every protein with
both small-GTPase motifs (P-loop `GxxxxGK[ST]` and G3 `DxxG[QH]`). `make_reference.py` builds the list.
CB513 itself contains an H-Ras structure (a G12P mutant, residues 1–166), which stays in the test set.

# Reference structures

| File | Content | Source |
|---|---|---|
| `structures/hras_1crr.pdb` | H-Ras 1–166, wild type with GDP, NMR model 1 of 20, backbone atoms | PDB 1CRR (Kraulis et al. 1994), via `github.com/biotite-dev/biotite` `tests/structure/data/pdb/1crr.pdb` |
| `structures/hb_4hhb.pdb` | Human deoxyhaemoglobin, chains A/C = α, B/D = β, backbone atoms and the four haem irons | PDB 4HHB (Fermi et al. 1984, X-ray 1.74 Å), via `github.com/3dmol/3Dmol.js` `tests/auto/data/4hhb.bcif.gz` |
| `reference_ss.tsv` | 8-state secondary structure of H-Ras (1CRR) and haemoglobin α and β (4HHB), computed from the coordinates above | `make_reference.py` |

The secondary structure is computed with pydssp, a re-implementation of the Kabsch & Sander DSSP hydrogen-bond
rules, plus 3₁₀ (G) and π (I) helices from consecutive 3- and 5-turns. As a check, it agrees with the independent
DSSP labels in the datasets above on 93 % of haemoglobin β residues (training-set entry, a variant structure) and
on 87 % of H-Ras residues (CB513, a crystal structure of the G12P mutant; 1CRR is an NMR structure, which is
expected to differ a little more).

`make_reference.py` rebuilds all of these files from the sources (needs `pip install biotite pydssp biopython`).
You don't need to run it for the tutorial.
