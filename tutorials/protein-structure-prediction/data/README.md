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

Note: CB513 contains an H-Ras structure (a G12P mutant, residues 1–166). The tutorial's `nn` command removes
small-GTPase-like proteins (both a P-loop `GxxxxGK[ST]` and a G3 `DxxG[QH]` motif) from the training set, so the
network never sees a Ras relative before it predicts H-Ras.
