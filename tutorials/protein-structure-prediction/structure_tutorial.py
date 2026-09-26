#!/usr/bin/env python3
"""
Protein structure prediction tutorial -- companion script.

Walks from the codons in ../../triplets.txt to a predicted protein structure:

  1. orfs        translate each reading frame and find open reading frames
  2. splice      join the coding exons and translate the full protein
  3. analyze     sequence-based predictions for one ORF:
                   - composition and physicochemical properties
                   - Kyte-Doolittle hydropathy (buried vs exposed, TM helices)
                   - Chou-Fasman secondary structure (helix / strand / coil)
                   - low-complexity and tandem-repeat detection (disorder hints)
  3b. nn         train a small neural network (Qian & Sejnowski style) on
                 real DSSP data and compare it with Chou-Fasman on CB513
  4. fold        (optional, needs internet) submit the ORF to the ESMFold API
                 and read per-residue confidence (pLDDT) from the returned PDB
  5. plddt       summarise pLDDT from any AlphaFold/ESMFold/ColabFold PDB file

Only the Python standard library is needed, except `fold`, which uses urllib.

Examples:
  python3 structure_tutorial.py orfs
  python3 structure_tutorial.py analyze            # longest ORF overall
  python3 structure_tutorial.py analyze --frame 3  # longest ORF in frame 3
  python3 structure_tutorial.py splice             # join exons -> full protein
  python3 structure_tutorial.py analyze --spliced
  python3 structure_tutorial.py fasta --spliced > protein.fasta
  python3 structure_tutorial.py nn                 # about 15 seconds
  python3 structure_tutorial.py fold --out orf.pdb
  python3 structure_tutorial.py plddt orf.pdb
"""

import argparse
import os
import re
import sys
import urllib.request
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_TRIPLETS = os.path.join(HERE, "..", "..", "triplets.txt")

# --------------------------------------------------------------------------
# Step 1-2: codons -> protein
# --------------------------------------------------------------------------

BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {
    a + b + c: AMINO[16 * i + 4 * j + k]
    for i, a in enumerate(BASES)
    for j, b in enumerate(BASES)
    for k, c in enumerate(BASES)
}


def read_frames(path):
    """Return {frame_number: [codon, ...]} from a triplets.txt-style file."""
    frames, current = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.lower().startswith("start"):
                current = int(line.split()[1])
                frames[current] = []
            elif current is not None:
                frames[current].extend(line.split())
    return frames


def translate(codons):
    # Incomplete trailing codons (e.g. "CC") become X.
    return "".join(CODON_TABLE.get(c.upper(), "X") for c in codons)


def find_orfs(protein, min_len=50):
    """ORFs = Met ... stop. Returns (start, end, seq) with 1-based residue start."""
    orfs = []
    for m in re.finditer(r"M[^*]*\*", protein):
        # re.finditer does not overlap, so each stop gets its first upstream Met,
        # which is the longest ORF ending at that stop.
        seq = m.group()[:-1]
        if len(seq) >= min_len:
            orfs.append((m.start() + 1, m.end() - 1, seq))
    return orfs


def all_orfs(path, min_len=50):
    out = []
    for frame, codons in sorted(read_frames(path).items()):
        for start, end, seq in find_orfs(translate(codons), min_len):
            out.append({"frame": frame, "start": start, "end": end, "seq": seq})
    return out


# Coding exons found in Part 2 of the tutorial (1-based, inclusive, on the
# forward strand = frame 1 read as one continuous DNA string). Each intron
# starts with GT and ends with AG; the last exon includes the stop codon.
HRAS_EXONS = [(6229, 6339), (6607, 6785), (6939, 7098), (7796, 7915)]


def read_dna(path):
    return "".join(read_frames(path)[1]).upper()


def splice(path, exons=HRAS_EXONS):
    dna = read_dna(path)
    cds = "".join(dna[a - 1:b] for a, b in exons)
    return cds, translate(cds[i:i + 3] for i in range(0, len(cds) - 2, 3))


def pick_orf(path, frame=None):
    orfs = all_orfs(path, min_len=1)
    if frame is not None:
        orfs = [o for o in orfs if o["frame"] == frame]
    if not orfs:
        sys.exit("No ORF found.")
    return max(orfs, key=lambda o: len(o["seq"]))


# --------------------------------------------------------------------------
# Step 3: sequence-based structure predictions
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

# Residues enriched in intrinsically disordered regions vs. ordered cores.
DISORDER_PROMOTING = set("PESQKAG")
ORDER_PROMOTING = set("WCFIYVLN")

# Average residue masses (Da) for a rough molecular weight.
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


def tandem_repeats(seq, min_unit=5, max_unit=40, min_copies=3):
    """Find the period with the most exact back-to-back repeats (a simple scan)."""
    best = None
    for unit in range(min_unit, max_unit + 1):
        for i in range(len(seq) - unit * min_copies + 1):
            motif = seq[i:i + unit]
            copies, j = 1, i + unit
            while seq[j:j + unit] == motif:
                copies += 1
                j += unit
            if copies >= min_copies and (best is None or copies * unit > best[2] * best[1]):
                best = (i + 1, unit, copies, motif)
    return best


def approximate_period(seq, max_unit=60):
    """Period with the highest self-identity when the sequence is shifted by it."""
    scores = []
    for p in range(3, min(max_unit, len(seq) // 2) + 1):
        same = sum(1 for i in range(len(seq) - p) if seq[i] == seq[i + p])
        scores.append((same / (len(seq) - p), p))
    return max(scores)


def shannon_entropy(s):
    import math
    counts = Counter(s)
    n = len(s)
    return -sum(c / n * math.log2(c / n) for c in counts.values())


def low_complexity(seq, window=12, cutoff=2.2):
    """Mask windows whose Shannon entropy (bits) is low (SEG-like, simplified)."""
    mask = [False] * len(seq)
    for i in range(len(seq) - window + 1):
        if shannon_entropy(seq[i:i + window]) < cutoff:
            for j in range(i, i + window):
                mask[j] = True
    return mask


def wrap(s, width=60):
    return [s[i:i + width] for i in range(0, len(s), width)]


# Approximate secondary structure of H-Ras in crystal structure PDB 5P21
# (residues 1-166; 167-189 are not resolved). Used to score Chou-Fasman.
HRAS_SS_5P21 = {
    "E": [(2, 9), (37, 46), (49, 58), (77, 83), (111, 116), (141, 143)],
    "H": [(16, 25), (66, 74), (87, 104), (127, 137), (152, 166)],
}


def reference_ss(length=166, ranges=HRAS_SS_5P21):
    ss = ["C"] * length
    for state, segs in ranges.items():
        for a, b in segs:
            for i in range(a - 1, b):
                ss[i] = state
    return "".join(ss)


def q3(pred, ref):
    n = min(len(pred), len(ref))
    return sum(1 for a, b in zip(pred[:n], ref[:n]) if a == b) / n


def analyze(seq, label, reference=None):
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

    # Hydropathy: 19-residue window average > 1.6 suggests a TM helix.
    hyd = sliding_mean(seq, KD, 19)
    tm = [i for i, v in enumerate(hyd) if v is not None and v > 1.6]
    if tm:
        segs, s = [], tm[0]
        for a, b in zip(tm, tm[1:] + [None]):
            if b != a + 1:
                segs.append((s + 1, a + 1))
                s = b
        print("Possible transmembrane segments (KD window 19 > 1.6):",
              ", ".join(f"{a}-{b}" for a, b in segs))
    else:
        print("No transmembrane helix predicted (no 19-residue window with KD > 1.6).")
    peak = max(v for v in hyd if v is not None) if n >= 19 else float("nan")
    print(f"Max windowed hydropathy: {peak:+.2f}")
    print()

    ss = chou_fasman(seq)
    lc = low_complexity(seq)
    lc_str = "".join("x" if m else "." for m in lc)
    print("Secondary structure (Chou-Fasman, simplified):  H helix  E strand  C coil")
    print("Low complexity mask:                             x low-entropy window")
    for i, (a, b, c) in enumerate(zip(wrap(seq), wrap(ss), wrap(lc_str))):
        print(f"{i * 60 + 1:>5} seq {a}")
        print(f"      ss  {b}")
        print(f"      lc  {c}")
    print()
    print(f"Helix {ss.count('H') / n:.0%}   Strand {ss.count('E') / n:.0%}   "
          f"Coil {ss.count('C') / n:.0%}   Low-complexity {sum(lc) / n:.0%}")

    rep = tandem_repeats(seq)
    ident, period = approximate_period(seq)
    print()
    if rep:
        start, unit, copies, motif = rep
        print(f"Exact tandem repeat: '{motif}' x{copies} (unit {unit} aa, starts at {start})")
    print(f"Best approximate period: {period} aa "
          f"({ident:.0%} of residues match the residue {period} positions later)")
    if ident > 0.4:
        print("  -> strongly repetitive; expect low-confidence / disordered prediction")

    if reference:
        print()
        print("Comparison with the crystal structure (PDB 5P21, residues 1-166):")
        for i, (a, b) in enumerate(zip(wrap(ss[:len(reference)]), wrap(reference))):
            print(f"{i * 60 + 1:>5} pred {a}")
            print(f"      5P21 {b}")
        print(f"Q3 accuracy (fraction of residues with the correct state): "
              f"{q3(ss, reference):.0%}")


# --------------------------------------------------------------------------
# Step 3b: a neural network for secondary structure (Qian & Sejnowski 1988)
# --------------------------------------------------------------------------

DATA_DIR = os.path.join(HERE, "data")
DSSP_TO_Q3 = {"H": "H", "G": "H", "I": "H", "E": "E", "B": "E"}  # rest -> C
NN_AA = "ACDEFGHIKLMNPQRSTVWY"   # index 20 = unknown (X), 21 = past the chain end
NN_SS = "HEC"
GTPASE_PLOOP = re.compile(r"G.{4}GK[ST]")
GTPASE_G3 = re.compile(r"D.{2}G[QH]")


def load_ss_dataset(name):
    """Return [(sequence, 3-state labels)] from data/<name>.tsv."""
    out = []
    with open(os.path.join(DATA_DIR, name + ".tsv")) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            seq, dssp = line.rstrip("\n").split("\t")
            out.append((seq, "".join(DSSP_TO_Q3.get(c, "C") for c in dssp)))
    return out


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


def run_nn(window, hidden, n_proteins, epochs, lr, seed, hras_seq):
    import random
    import time
    train = load_ss_dataset("cb6133filtered")
    test = load_ss_dataset("cb513")
    kept = [t for t in train if not (GTPASE_PLOOP.search(t[0]) and GTPASE_G3.search(t[0]))]
    random.Random(seed).shuffle(kept)
    subset = kept[:n_proteins]
    examples = [(row, NN_SS.index(y)) for seq, ss in subset
                for row, y in zip(encode_windows(seq, window), ss)]
    print(f"Training set: {len(subset)} proteins, {len(examples)} residues "
          f"({len(train) - len(kept)} small-GTPase-like proteins excluded)")
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

    ref = reference_ss()
    nn_ss, cf_ss = net.predict(hras_seq), chou_fasman(hras_seq)
    print("\nH-Ras, residues 1-166 (not in the training set):")
    for i in range(0, 166, 60):
        print(f"{i + 1:>5} seq  {hras_seq[i:min(i + 60, 166)]}")
        print(f"      CF   {cf_ss[i:min(i + 60, 166)]}")
        print(f"      NN   {nn_ss[i:min(i + 60, 166)]}")
        print(f"      5P21 {ref[i:i + 60]}")
    print(f"Q3 vs 5P21: Chou-Fasman {q3(cf_ss, ref):.0%}, neural network {q3(nn_ss, ref):.0%}")


# --------------------------------------------------------------------------
# Step 4-5: 3D prediction with ESMFold and reading pLDDT
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
    vals = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            vals.append((int(line[22:26]), float(line[60:66])))
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
    # One character per residue: 9 = 90+, 7 = 70-89, 5 = 50-69, . = <50
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
    ap.add_argument("--triplets", default=DEFAULT_TRIPLETS)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("orfs", help="list ORFs in all frames")
    p.add_argument("--min-len", type=int, default=50)

    for name in ("analyze", "fasta", "fold"):
        p = sub.add_parser(name)
        p.add_argument("--frame", type=int, help="restrict to one reading frame")
        p.add_argument("--seq", help="analyse this protein sequence instead")
        p.add_argument("--spliced", action="store_true",
                       help="use the protein from the spliced exons (see `splice`)")
        if name == "fold":
            p.add_argument("--out", default="prediction.pdb")

    sub.add_parser("splice", help="join the exons and translate the full protein")

    p = sub.add_parser("nn", help="train a window neural network for secondary structure")
    p.add_argument("--window", type=int, default=13)
    p.add_argument("--hidden", type=int, default=10)
    p.add_argument("--proteins", type=int, default=300, help="training proteins to use")
    p.add_argument("--epochs", type=int, default=3)
    p.add_argument("--lr", type=float, default=0.02)
    p.add_argument("--seed", type=int, default=1)

    p = sub.add_parser("plddt", help="summarise confidence in a predicted PDB")
    p.add_argument("pdb")

    args = ap.parse_args()

    if args.cmd == "orfs":
        orfs = all_orfs(args.triplets, args.min_len)
        print(f"{'frame':>5} {'start':>6} {'end':>6} {'length':>6}  first 40 aa")
        for o in sorted(orfs, key=lambda o: -len(o["seq"])):
            print(f"{o['frame']:>5} {o['start']:>6} {o['end']:>6} "
                  f"{len(o['seq']):>6}  {o['seq'][:40]}")
        return

    if args.cmd == "plddt":
        with open(args.pdb) as fh:
            summarise_plddt(plddt_from_pdb(fh.read()))
        return

    if args.cmd == "nn":
        run_nn(args.window, args.hidden, args.proteins, args.epochs, args.lr, args.seed,
               splice(args.triplets)[1].rstrip("*"))
        return

    if args.cmd == "splice":
        dna = read_dna(args.triplets)
        cds, prot = splice(args.triplets)
        for n, ((a, b), nxt) in enumerate(zip(HRAS_EXONS, HRAS_EXONS[1:] + [None]), 1):
            line = f"exon {n}: {a:>5}-{b:<5} ({b - a + 1:>3} nt)"
            if nxt:
                intron = dna[b:nxt[0] - 1]
                line += f"   intron {len(intron):>4} nt  {intron[:2]}...{intron[-2:]}"
            print(line)
        print(f"\nCDS {len(cds)} nt -> {len(prot.rstrip('*'))} aa, "
              f"ends with stop: {prot.endswith('*')}\n")
        print("\n".join(wrap(prot.rstrip("*"))))
        return

    if args.seq:
        seq, label = args.seq.upper(), "user sequence"
    elif args.spliced:
        seq, label = splice(args.triplets)[1].rstrip("*"), "spliced protein (exons 1-4)"
    else:
        o = pick_orf(args.triplets, args.frame)
        seq = o["seq"]
        label = f"frame {o['frame']} ORF, residues {o['start']}-{o['end']}"

    if args.cmd == "analyze":
        analyze(seq, label, reference_ss() if args.spliced else None)
    elif args.cmd == "fasta":
        print(f">{label.replace(' ', '_').replace(',', '')}")
        print("\n".join(wrap(seq)))
    elif args.cmd == "fold":
        print(f"Submitting {len(seq)} aa to ESMFold ...", file=sys.stderr)
        try:
            pdb = esmfold(seq)
        except Exception as e:  # network blocked, rate limited, sequence too long
            sys.exit(f"ESMFold request failed: {e}\n"
                     "Use the FASTA with ColabFold or the ESMFold web page instead "
                     "(see README, Part 5).")
        with open(args.out, "w") as fh:
            fh.write(pdb)
        print(f"Wrote {args.out}", file=sys.stderr)
        summarise_plddt(plddt_from_pdb(pdb))


if __name__ == "__main__":
    main()
