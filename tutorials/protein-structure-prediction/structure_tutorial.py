#!/usr/bin/env python3
"""
Protein structure prediction tutorial -- companion script.

Start from a protein sequence and ask what can be predicted about its
structure. Two worked examples: H-Ras and human haemoglobin (alpha and beta).

  proteins   list the example proteins and their reference structures
  fasta      print an example sequence in FASTA format
  analyze    sequence-only predictions: composition, Kyte-Doolittle hydropathy,
             Chou-Fasman secondary structure, scored against the experimental
             structure (DSSP)
  nn         train a small neural network (Qian & Sejnowski style) on real
             DSSP data and compare it with Chou-Fasman
  compare    superimpose a predicted model on the experimental structure and
             report which residues are correct (within 1, 2, 4, 8 A)
  fold       (needs internet) submit a sequence to the ESMFold API
  plddt      summarise per-residue confidence in a predicted PDB file

Only the Python standard library is needed (fold uses urllib).

Examples:
  python3 structure_tutorial.py proteins
  python3 structure_tutorial.py analyze --protein hras
  python3 structure_tutorial.py analyze --protein hbb
  python3 structure_tutorial.py nn                    # about 15 seconds
  python3 structure_tutorial.py fasta --protein hbb > hbb.fasta
  python3 structure_tutorial.py compare my_model.pdb --protein hbb
  python3 structure_tutorial.py compare data/structures/hb_4hhb.pdb --chain A --protein hbb
  python3 structure_tutorial.py plddt my_model.pdb
"""

import argparse
import math
import os
import random
import re
import sys
import time
import urllib.request
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")

# --------------------------------------------------------------------------
# The example proteins
# --------------------------------------------------------------------------

PROTEINS = {
    "hras": {
        "name": "H-Ras (GTPase HRas), human",
        "uniprot": "P01112",
        # Full-length protein. Residues 1-166 are identical to the SEQRES of PDB 1CRR;
        # the whole sequence also matches the translation of the H-Ras gene in ../../triplets.txt.
        "seq": "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINN"
               "TKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLN"
               "PPDESGPGCMSCKCVLS",
        "structure": "hras_1crr.pdb", "chain": "A",
        "about": "189 aa; G domain 1-166 (alpha/beta fold) + flexible C-terminal tail. "
                 "Reference: PDB 1CRR, NMR model 1 (wild type, GDP-bound), residues 1-166.",
    },
    "hbb": {
        "name": "Haemoglobin subunit beta, human",
        "uniprot": "P68871",
        # Mature chain as in PDB 4HHB chain B (the initiator Met of the 147-aa precursor is removed).
        "seq": "VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFAT"
               "LSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH",
        "structure": "hb_4hhb.pdb", "chain": "B",
        "about": "146 aa; all-alpha globin fold, binds one haem; part of the alpha2beta2 tetramer. "
                 "Reference: PDB 4HHB chain B (deoxyhaemoglobin, X-ray 1.74 A).",
    },
    "hba": {
        "name": "Haemoglobin subunit alpha, human",
        "uniprot": "P69905",
        "seq": "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDL"
               "HAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
        "structure": "hb_4hhb.pdb", "chain": "A",
        "about": "141 aa; globin fold, 43% identical to beta. Reference: PDB 4HHB chain A.",
    },
}


def load_reference_ss():
    """{protein: (first_residue, sequence, dssp8)} from data/reference_ss.tsv."""
    out = {}
    with open(os.path.join(DATA_DIR, "reference_ss.tsv")) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            name, _pdb, _chain, first, seq, ss = line.rstrip("\n").split("\t")
            out[name] = (int(first), seq, ss)
    return out


DSSP_TO_Q3 = {"H": "H", "G": "H", "I": "H", "E": "E", "B": "E"}  # everything else -> C


def to_q3(dssp):
    return "".join(DSSP_TO_Q3.get(c, "C") for c in dssp)


def reference_q3(protein):
    """3-state reference aligned to the protein sequence (covers residues 1..len(reference))."""
    first, seq, ss = load_reference_ss()[protein]
    assert PROTEINS[protein]["seq"][first - 1:first - 1 + len(seq)] == seq, "reference does not match sequence"
    return to_q3(ss)


def q3(pred, ref):
    n = min(len(pred), len(ref))
    return sum(1 for a, b in zip(pred[:n], ref[:n]) if a == b) / n


def wrap(s, width=60):
    return [s[i:i + width] for i in range(0, len(s), width)]


# --------------------------------------------------------------------------
# Sequence-only predictions
# --------------------------------------------------------------------------

# Kyte & Doolittle (1982) hydropathy scale.
KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Chou & Fasman (1978) propensities: (P_alpha, P_beta).
CHOU_FASMAN = {
    "E": (1.51, 0.37), "M": (1.45, 1.05), "A": (1.42, 0.83), "L": (1.21, 1.30),
    "K": (1.16, 0.74), "F": (1.13, 1.38), "Q": (1.11, 1.10), "W": (1.08, 1.37),
    "I": (1.08, 1.60), "V": (1.06, 1.70), "D": (1.01, 0.54), "H": (1.00, 0.87),
    "R": (0.98, 0.93), "T": (0.83, 1.19), "S": (0.77, 0.75), "C": (0.70, 1.19),
    "Y": (0.69, 1.47), "N": (0.67, 0.89), "P": (0.57, 0.55), "G": (0.57, 0.75),
}

DISORDER_PROMOTING = set("PESQKAG")
ORDER_PROMOTING = set("WCFIYVLN")

MASS = {
    "A": 71.08, "R": 156.19, "N": 114.10, "D": 115.09, "C": 103.14, "E": 129.12,
    "Q": 128.13, "G": 57.05, "H": 137.14, "I": 113.16, "L": 113.16, "K": 128.17,
    "M": 131.19, "F": 147.18, "P": 97.12, "S": 87.08, "T": 101.10, "W": 186.21,
    "Y": 163.18, "V": 99.13,
}


def sliding_mean(seq, scale, window):
    half = window // 2
    vals = [scale.get(a, 0.0) for a in seq]
    out = [None] * len(seq)
    for i in range(half, len(seq) - half):
        out[i] = sum(vals[i - half:i + half + 1]) / window
    return out


def chou_fasman(seq, window=6):
    """
    Simplified Chou-Fasman: average helix/strand propensity over a window,
    call H if <P_alpha> >= 1.03 and beats <P_beta>, E if <P_beta> >= 1.05,
    otherwise C (coil). Real Chou-Fasman adds nucleation/extension rules;
    this version is deliberately short so students can read it.
    """
    pa = {a: v[0] for a, v in CHOU_FASMAN.items()}
    pb = {a: v[1] for a, v in CHOU_FASMAN.items()}
    ha = sliding_mean(seq, pa, window + 1)
    hb = sliding_mean(seq, pb, window + 1)
    ss = []
    for a, b in zip(ha, hb):
        if a is None:
            ss.append("C")
        elif a >= 1.03 and a > b:
            ss.append("H")
        elif b >= 1.05 and b >= a:
            ss.append("E")
        else:
            ss.append("C")
    # Remove isolated calls shorter than 4 (helix) / 3 (strand).
    s = "".join(ss)
    s = re.sub(r"(?<!H)H{1,3}(?!H)", lambda m: "C" * len(m.group()), s)
    s = re.sub(r"(?<!E)E{1,2}(?!E)", lambda m: "C" * len(m.group()), s)
    return s


def print_tracks(seq, tracks, ref=None):
    """tracks: [(label, string)]; residues that disagree with ref are shown in lower case."""
    for i in range(0, len(seq), 60):
        print(f"{i + 1:>5} seq  {seq[i:i + 60]}")
        for lab, s in tracks:
            part = s[i:i + 60]
            if ref is not None and lab != "DSSP":
                part = "".join(c if i + j >= len(ref) or c == ref[i + j] else c.lower()
                               for j, c in enumerate(part))
            if part:
                print(f"      {lab:<4} {part}")
        print()


def analyze(seq, label, ref=None):
    n = len(seq)
    comp = Counter(seq)
    mw = sum(MASS.get(a, 110) for a in seq) + 18.02
    gravy = sum(KD.get(a, 0) for a in seq) / n
    dis = sum(comp[a] for a in DISORDER_PROMOTING) / n
    order = sum(comp[a] for a in ORDER_PROMOTING) / n
    charge = comp["K"] + comp["R"] - comp["D"] - comp["E"]

    print(f"== {label} ==")
    print(f"Length            {n} aa")
    print(f"Mol. weight       {mw / 1000:.1f} kDa")
    print(f"GRAVY             {gravy:+.2f}   (>0 hydrophobic, <0 hydrophilic)")
    print(f"Net charge ~pH7   {charge:+d}")
    print(f"Disorder-promoting residues (PESQKAG): {dis:.0%}")
    print(f"Order-promoting residues   (WCFIYVLN): {order:.0%}")
    top = ", ".join(f"{a} {c / n:.0%}" for a, c in comp.most_common(5))
    print(f"Most common       {top}")
    print()

    hyd = sliding_mean(seq, KD, 19)
    tm = [i for i, v in enumerate(hyd) if v is not None and v > 1.6]
    peak = max(v for v in hyd if v is not None) if n >= 19 else float("nan")
    if tm:
        print(f"Possible transmembrane segment(s): windows centred on residues {tm[0] + 1}-{tm[-1] + 1} "
              f"(KD window 19 > 1.6)")
    else:
        print("No transmembrane helix predicted (no 19-residue window with KD > 1.6).")
    print(f"Max windowed hydropathy: {peak:+.2f}")
    print()

    ss = chou_fasman(seq)
    print("Secondary structure: H helix, E strand, C coil.")
    if ref:
        print("DSSP = experimental structure; lower-case letters disagree with it.")
    print()
    print_tracks(seq, [("CF", ss)] + ([("DSSP", ref)] if ref else []), ref)
    print(f"Chou-Fasman: helix {ss.count('H') / n:.0%}, strand {ss.count('E') / n:.0%}, "
          f"coil {ss.count('C') / n:.0%}")
    if ref:
        print(f"Experiment:  helix {ref.count('H') / len(ref):.0%}, strand {ref.count('E') / len(ref):.0%}, "
              f"coil {ref.count('C') / len(ref):.0%}  (residues 1-{len(ref)})")
        print(f"Q3 (fraction of residues in the correct state): {q3(ss, ref):.0%}")


# --------------------------------------------------------------------------
# A neural network for secondary structure (Qian & Sejnowski 1988)
# --------------------------------------------------------------------------

NN_AA = "ACDEFGHIKLMNPQRSTVWY"   # index 20 = unknown (X), 21 = past the chain end
NN_SS = "HEC"


def load_ss_dataset(name):
    """Return [(sequence, 3-state labels)] from data/<name>.tsv."""
    out = []
    with open(os.path.join(DATA_DIR, name + ".tsv")) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            seq, dssp = line.rstrip("\n").split("\t")
            out.append((seq, to_q3(dssp)))
    return out


def excluded_training():
    """Indices of training proteins related to the example proteins (see data/excluded_training.txt)."""
    with open(os.path.join(DATA_DIR, "excluded_training.txt")) as fh:
        return {int(l.split("\t")[0]) for l in fh if l.strip() and not l.startswith("#")}


def encode_windows(seq, window):
    """
    One-hot encoding, stored sparsely: for each residue, the list of active
    input units (one per window position). Input unit = position * 22 + letter.
    """
    half = window // 2
    idx = [NN_AA.find(a) if a in NN_AA else 20 for a in seq]
    rows = []
    for i in range(len(seq)):
        row = []
        for w in range(window):
            j = i - half + w
            row.append(w * 22 + (idx[j] if 0 <= j < len(seq) else 21))
        rows.append(row)
    return rows


class WindowNet:
    """
    Sequence window -> one-hot -> tanh hidden layer -> softmax over H/E/C.
    With hidden=0 it is a single-layer (linear) network.
    Trained with plain stochastic gradient descent on cross-entropy.
    """

    def __init__(self, window=13, hidden=10, seed=1):
        import random
        self.window, self.hidden = window, hidden
        self.rng = random.Random(seed)
        n_in = 22 * window
        r = self.rng.uniform
        if hidden:
            self.w1 = [[r(-0.1, 0.1) for _ in range(hidden)] for _ in range(n_in)]
            self.b1 = [0.0] * hidden
            self.w2 = [[r(-0.1, 0.1) for _ in range(3)] for _ in range(hidden)]
        else:
            self.w2 = [[0.0] * 3 for _ in range(n_in)]
        self.b2 = [0.0] * 3

    def _forward(self, row):
        import math
        if self.hidden:
            a = self.b1[:]
            for u in row:
                wu = self.w1[u]
                for k in range(self.hidden):
                    a[k] += wu[k]
            a = [math.tanh(v) for v in a]
            z = [self.b2[c] + sum(a[k] * self.w2[k][c] for k in range(self.hidden)) for c in range(3)]
        else:
            a = None
            z = [self.b2[c] + sum(self.w2[u][c] for u in row) for c in range(3)]
        m = max(z)
        e = [math.exp(v - m) for v in z]
        t = sum(e)
        return a, [v / t for v in e]

    def predict(self, seq):
        return "".join(NN_SS[max(range(3), key=p.__getitem__)]
                       for _, p in map(self._forward, encode_windows(seq, self.window)))

    def train_epoch(self, examples, lr):
        self.rng.shuffle(examples)
        correct = 0
        for row, y in examples:
            a, p = self._forward(row)
            correct += max(range(3), key=p.__getitem__) == y
            g = p[:]
            g[y] -= 1.0                      # dLoss/dz for softmax + cross-entropy
            if self.hidden:
                da = [sum(self.w2[k][c] * g[c] for c in range(3)) * (1 - a[k] * a[k])
                      for k in range(self.hidden)]
                for k in range(self.hidden):
                    for c in range(3):
                        self.w2[k][c] -= lr * a[k] * g[c]
                for u in row:
                    wu = self.w1[u]
                    for k in range(self.hidden):
                        wu[k] -= lr * da[k]
                for k in range(self.hidden):
                    self.b1[k] -= lr * da[k]
            else:
                for u in row:
                    for c in range(3):
                        self.w2[u][c] -= lr * g[c]
            for c in range(3):
                self.b2[c] -= lr * g[c]
        return correct / len(examples)


def score_q3(predict, data):
    ok = n = 0
    for seq, ss in data:
        pred = predict(seq)
        ok += sum(a == b for a, b in zip(pred, ss))
        n += len(ss)
    return ok / n


def run_nn(window, hidden, n_proteins, epochs, lr, seed):
    train = load_ss_dataset("cb6133filtered")
    test = load_ss_dataset("cb513")
    skip = excluded_training()
    kept = [t for i, t in enumerate(train) if i not in skip]
    random.Random(seed).shuffle(kept)
    subset = kept[:n_proteins]
    examples = [(row, NN_SS.index(y)) for seq, ss in subset
                for row, y in zip(encode_windows(seq, window), ss)]
    print(f"Training set: {len(subset)} proteins, {len(examples)} residues "
          f"({len(skip)} relatives of Ras and haemoglobin excluded)")
    print(f"Test set (CB513): {len(test)} proteins, {sum(len(s) for s, _ in test)} residues")
    print(f"Network: window {window} x 22 inputs -> {hidden or 'no'} hidden units -> 3 outputs\n")

    cf = score_q3(chou_fasman, test)
    net = WindowNet(window, hidden, seed)
    print(f"{'epoch':>5} {'train Q3':>9} {'test Q3':>8} {'time':>6}")
    for ep in range(1, epochs + 1):
        t0 = time.time()
        tr = net.train_epoch(examples, lr / ep ** 0.5)
        te = score_q3(net.predict, test)
        print(f"{ep:>5} {tr:>9.1%} {te:>8.1%} {time.time() - t0:>5.0f}s")
    print(f"\nChou-Fasman on the same test set: {cf:.1%}")

    for name in ("hras", "hbb"):
        seq, ref = PROTEINS[name]["seq"], reference_q3(name)
        cf_ss, nn_ss = chou_fasman(seq), net.predict(seq)
        print(f"\n{PROTEINS[name]['name']} (no relatives in the training set); "
              f"lower case = disagrees with the experimental structure")
        print_tracks(seq[:len(ref)], [("CF", cf_ss[:len(ref)]), ("NN", nn_ss[:len(ref)]), ("DSSP", ref)], ref)
        print(f"Q3: Chou-Fasman {q3(cf_ss, ref):.0%}, neural network {q3(nn_ss, ref):.0%}")


# --------------------------------------------------------------------------
# Comparing a 3D model with the experimental structure
# --------------------------------------------------------------------------

THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
         "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
         "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M"}


def read_ca(pdb_text):
    """{chain: [(resnum, aa, (x, y, z), bfactor)]} for CA atoms of the first model."""
    chains = {}
    for line in pdb_text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA" and line[16] in " A":
            aa = THREE.get(line[17:20].strip())
            if aa:
                chains.setdefault(line[21], []).append(
                    (int(line[22:26]), aa, (float(line[30:38]), float(line[38:46]), float(line[46:54])),
                     float(line[60:66] or 0)))
    return chains


def align_sequences(a, b, gap=-4):
    """Needleman-Wunsch with a simple identity score; returns [(i, j)] of aligned positions."""
    n, m = len(a), len(b)
    S = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        S[i][0] = i * gap
    for j in range(1, m + 1):
        S[0][j] = j * gap
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            S[i][j] = max(S[i - 1][j - 1] + (5 if a[i - 1] == b[j - 1] else -2),
                          S[i - 1][j] + gap, S[i][j - 1] + gap)
    pairs, i, j = [], n, m
    while i and j:
        if S[i][j] == S[i - 1][j - 1] + (5 if a[i - 1] == b[j - 1] else -2):
            pairs.append((i - 1, j - 1))
            i, j = i - 1, j - 1
        elif S[i][j] == S[i - 1][j] + gap:
            i -= 1
        else:
            j -= 1
    return pairs[::-1]


def jacobi_eigen(A, sweeps=60):
    """Eigen-decomposition of a small symmetric matrix. Returns (eigenvalues, eigenvector columns)."""
    n = len(A)
    a = [row[:] for row in A]
    v = [[float(i == j) for j in range(n)] for i in range(n)]
    for _ in range(sweeps):
        off = sum(a[p][q] ** 2 for p in range(n) for q in range(p + 1, n))
        if off < 1e-20:
            break
        for p in range(n):
            for q in range(p + 1, n):
                if abs(a[p][q]) < 1e-15:
                    continue
                th = (a[q][q] - a[p][p]) / (2 * a[p][q])
                t = (1 if th >= 0 else -1) / (abs(th) + math.sqrt(th * th + 1))
                c = 1 / math.sqrt(t * t + 1)
                s = t * c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p], a[k][q] = c * akp - s * akq, s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k], a[q][k] = c * apk - s * aqk, s * apk + c * aqk
                for k in range(n):
                    vkp, vkq = v[k][p], v[k][q]
                    v[k][p], v[k][q] = c * vkp - s * vkq, s * vkp + c * vkq
    return [a[i][i] for i in range(n)], v


def superpose(X, Y):
    """
    Optimal rotation + translation moving points X onto Y (Horn's quaternion method).
    Returns a function that transforms a point.
    """
    n = len(X)
    cx = [sum(p[k] for p in X) / n for k in range(3)]
    cy = [sum(p[k] for p in Y) / n for k in range(3)]
    S = [[0.0] * 3 for _ in range(3)]
    for p, q in zip(X, Y):
        u = [p[k] - cx[k] for k in range(3)]
        w = [q[k] - cy[k] for k in range(3)]
        for i in range(3):
            for j in range(3):
                S[i][j] += u[i] * w[j]
    (xx, xy, xz), (yx, yy, yz), (zx, zy, zz) = S
    N = [[xx + yy + zz, yz - zy, zx - xz, xy - yx],
         [yz - zy, xx - yy - zz, xy + yx, zx + xz],
         [zx - xz, xy + yx, -xx + yy - zz, yz + zy],
         [xy - yx, zx + xz, yz + zy, -xx - yy + zz]]
    vals, vecs = jacobi_eigen(N)
    k = max(range(4), key=vals.__getitem__)
    q0, q1, q2, q3_ = (vecs[i][k] for i in range(4))
    R = [[q0 * q0 + q1 * q1 - q2 * q2 - q3_ * q3_, 2 * (q1 * q2 - q0 * q3_), 2 * (q1 * q3_ + q0 * q2)],
         [2 * (q1 * q2 + q0 * q3_), q0 * q0 - q1 * q1 + q2 * q2 - q3_ * q3_, 2 * (q2 * q3_ - q0 * q1)],
         [2 * (q1 * q3_ - q0 * q2), 2 * (q2 * q3_ + q0 * q1), q0 * q0 - q1 * q1 - q2 * q2 + q3_ * q3_]]

    def move(p):
        u = [p[k] - cx[k] for k in range(3)]
        return tuple(sum(R[i][j] * u[j] for j in range(3)) + cy[i] for i in range(3))
    return move


def compare_structures(model, ref, core_cutoff=4.0, rounds=5):
    """
    model, ref: lists of (resnum, aa, xyz, b) for one chain each.
    Aligns the sequences, superimposes on all matched CA atoms, then re-fits on
    the residues within core_cutoff A (a few rounds) so that one floppy loop
    cannot drag the whole superposition. Returns per-residue distances.
    """
    pairs = align_sequences("".join(r[1] for r in model), "".join(r[1] for r in ref))
    use = list(range(len(pairs)))
    for _ in range(rounds):
        move = superpose([model[pairs[k][0]][2] for k in use], [ref[pairs[k][1]][2] for k in use])
        dist = [math.dist(move(model[i][2]), ref[j][2]) for i, j in pairs]
        new = [k for k, d in enumerate(dist) if d < core_cutoff]
        if len(new) < 3 or new == use:
            break
        use = new
    return pairs, dist, move


def report_comparison(model_chain, ref_chain, ref_label):
    pairs, dist, _ = compare_structures(model_chain, ref_chain)
    n_ref = len(ref_chain)
    ident = sum(model_chain[i][1] == ref_chain[j][1] for i, j in pairs)
    rmsd = math.sqrt(sum(d * d for d in dist) / len(dist))
    within = {t: sum(d <= t for d in dist) / n_ref for t in (1, 2, 4, 8)}
    print(f"Reference: {ref_label}, {n_ref} residues")
    print(f"Matched residues: {len(pairs)} ({ident} identical in sequence, {ident / len(pairs):.0%})")
    print(f"Cα RMSD over all matched residues: {rmsd:.2f} A")
    core = [d for d in dist if d <= 4]
    if core:
        print(f"Cα RMSD over the {len(core)} residues within 4 A: {math.sqrt(sum(d * d for d in core) / len(core)):.2f} A")
    print("Fraction of reference residues within  " + "   ".join(f"{t} A: {v:.0%}" for t, v in within.items()))
    print(f"GDT-style score (mean of the four fractions): {100 * sum(within.values()) / 4:.1f}")
    marks = ["-"] * n_ref
    for (i, j), d in zip(pairs, dist):
        marks[j] = "#" if d <= 2 else "+" if d <= 4 else "."
    ref_seq = "".join(r[1] for r in ref_chain)
    print("\nPer residue:  # within 2 A (correct)   + 2-4 A   . more than 4 A   - not in model")
    for k in range(0, n_ref, 60):
        print(f"{ref_chain[k][0]:>5} ref  {ref_seq[k:k + 60]}")
        print(f"      fit  {''.join(marks[k:k + 60])}\n")


def pick_chain(chains, target_seq, wanted=None):
    if wanted:
        if wanted not in chains:
            sys.exit(f"Chain {wanted!r} not found; chains in file: {', '.join(chains)}")
        return wanted
    def score(ch):
        pairs = align_sequences("".join(r[1] for r in chains[ch]), target_seq)
        return sum(chains[ch][i][1] == target_seq[j] for i, j in pairs)
    return max(chains, key=score)


# --------------------------------------------------------------------------
# ESMFold and pLDDT
# --------------------------------------------------------------------------

ESMFOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"


def esmfold(seq, timeout=300):
    """Send a sequence (<= ~400 aa) to the public ESMFold API; returns PDB text."""
    req = urllib.request.Request(ESMFOLD_URL, data=seq.encode(), method="POST")
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode()


def plddt_from_pdb(pdb_text):
    """
    AlphaFold, ColabFold and ESMFold store pLDDT in the B-factor column.
    Returns [(residue_number, plddt)] using the CA atom of each residue.
    ESMFold's API reports 0-1; AlphaFold reports 0-100 -- normalised to 0-100.
    """
    vals = [(r, b) for chain in read_ca(pdb_text).values() for r, _, _, b in chain]
    if vals and max(v for _, v in vals) <= 1.0:
        vals = [(r, v * 100) for r, v in vals]
    return vals


def summarise_plddt(vals):
    bands = [(90, "very high"), (70, "confident"), (50, "low"), (0, "very low")]
    n = len(vals)
    mean = sum(v for _, v in vals) / n
    print(f"Residues: {n}   mean pLDDT: {mean:.1f}")
    lo = 101
    for cut, name in bands:
        k = sum(1 for _, v in vals if cut <= v < lo)
        print(f"  {name:<10} ({cut:>2}-{lo - 1:<3}) {k:>5} residues  {k / n:.0%}")
        lo = cut
    track = "".join("9" if v >= 90 else "7" if v >= 70 else "5" if v >= 50 else "."
                    for _, v in vals)
    print("\nPer-residue confidence (9 >=90, 7 >=70, 5 >=50, . <50):")
    for i, line in enumerate(wrap(track)):
        print(f"{i * 60 + 1:>5} {line}")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("proteins", help="list the example proteins")

    for name in ("analyze", "fasta", "fold"):
        p = sub.add_parser(name)
        p.add_argument("--protein", choices=sorted(PROTEINS), default="hras")
        p.add_argument("--seq", help="use this protein sequence instead")
        if name == "fold":
            p.add_argument("--out", default="prediction.pdb")

    p = sub.add_parser("nn", help="train a window neural network for secondary structure")
    p.add_argument("--window", type=int, default=13)
    p.add_argument("--hidden", type=int, default=10)
    p.add_argument("--proteins", type=int, default=300, help="training proteins to use")
    p.add_argument("--epochs", type=int, default=3)
    p.add_argument("--lr", type=float, default=0.02)
    p.add_argument("--seed", type=int, default=1)

    p = sub.add_parser("compare", help="compare a predicted model with the experimental structure")
    p.add_argument("model", help="predicted structure (PDB format)")
    p.add_argument("--protein", choices=sorted(PROTEINS), default="hras")
    p.add_argument("--chain", help="chain of the model to use (default: best sequence match)")

    p = sub.add_parser("plddt", help="summarise confidence in a predicted PDB")
    p.add_argument("pdb")

    args = ap.parse_args()

    if args.cmd == "proteins":
        for key, p in PROTEINS.items():
            print(f"{key:<5} {p['name']} (UniProt {p['uniprot']}), {len(p['seq'])} aa")
            print(f"      {p['about']}")
        return
    if args.cmd == "nn":
        run_nn(args.window, args.hidden, args.proteins, args.epochs, args.lr, args.seed)
        return
    if args.cmd == "plddt":
        with open(args.pdb) as fh:
            summarise_plddt(plddt_from_pdb(fh.read()))
        return
    if args.cmd == "compare":
        prot = PROTEINS[args.protein]
        with open(os.path.join(DATA_DIR, "structures", prot["structure"])) as fh:
            ref = read_ca(fh.read())[prot["chain"]]
        with open(args.model) as fh:
            chains = read_ca(fh.read())
        if not chains:
            sys.exit("No CA atoms found in the model file.")
        ch = pick_chain(chains, "".join(r[1] for r in ref), args.chain)
        print(f"Model: {args.model}, chain {ch}, {len(chains[ch])} residues")
        report_comparison(chains[ch], ref, f"{prot['name']}, {prot['structure']} chain {prot['chain']}")
        return

    if args.seq:
        seq, label, ref = re.sub(r"[^A-Z]", "", args.seq.upper()), "your sequence", None
    else:
        prot = PROTEINS[args.protein]
        seq, label, ref = prot["seq"], f"{prot['name']} (UniProt {prot['uniprot']})", reference_q3(args.protein)

    if args.cmd == "analyze":
        analyze(seq, label, ref)
    elif args.cmd == "fasta":
        print(f">{args.protein if not args.seq else 'query'} {label}")
        print("\n".join(wrap(seq)))
    elif args.cmd == "fold":
        print(f"Submitting {len(seq)} aa to ESMFold ...", file=sys.stderr)
        try:
            pdb = esmfold(seq)
        except Exception as e:  # network blocked, rate limited, sequence too long
            sys.exit(f"ESMFold request failed: {e}\n"
                     "Use the FASTA with ColabFold or the ESMFold web page instead (see README, Part 4).")
        with open(args.out, "w") as fh:
            fh.write(pdb)
        print(f"Wrote {args.out}", file=sys.stderr)
        summarise_plddt(plddt_from_pdb(pdb))


if __name__ == "__main__":
    main()
