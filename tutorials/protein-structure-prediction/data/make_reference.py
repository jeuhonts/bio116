#!/usr/bin/env python3
"""
Rebuild the reference files in this folder from their public sources.

You do NOT need to run this for the tutorial; the outputs are committed. It
documents exactly where every coordinate and secondary-structure label comes
from. Needs: pip install biotite pydssp biopython numpy

Outputs
  structures/hras_1crr.pdb    H-Ras 1-166 (wild type, GDP), NMR model 1 of PDB 1CRR,
                              backbone atoms only
  structures/hb_4hhb.pdb      human deoxyhaemoglobin, PDB 4HHB, chains A-D backbone
                              atoms + the four haem iron atoms
  reference_ss.tsv            DSSP-style secondary structure (8 states) computed
                              from those coordinates
  excluded_training.txt       training proteins removed because they are relatives
                              of the example proteins

Sources (fetched from GitHub because RCSB was not reachable when this was built)
  1CRR: github.com/biotite-dev/biotite  tests/structure/data/pdb/1crr.pdb
  4HHB: github.com/3dmol/3Dmol.js       tests/auto/data/4hhb.bcif.gz
"""
import gzip
import io
import os
import re
import urllib.request
import warnings

import numpy as np

warnings.filterwarnings("ignore")
import biotite.structure as struc  # noqa: E402
import biotite.structure.io.pdb as pdbio  # noqa: E402
import biotite.structure.io.pdbx as pdbx  # noqa: E402
import pydssp  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
URL_1CRR = "https://raw.githubusercontent.com/biotite-dev/biotite/main/tests/structure/data/pdb/1crr.pdb"
URL_4HHB = "https://raw.githubusercontent.com/3dmol/3Dmol.js/master/tests/auto/data/4hhb.bcif.gz"


def fetch(url):
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read()


def backbone(arr, chain):
    a = arr[(arr.chain_id == chain) & struc.filter_amino_acids(arr) & (arr.element != "H")]
    coords, seq, ids = [], "", []
    for r in np.unique(a.res_id):
        res = a[a.res_id == r]
        coords.append([res[res.atom_name == n].coord[0] for n in ("N", "CA", "C", "O")])
        seq += struc.info.one_letter_code(res.res_name[0])
        ids.append(int(r))
    return np.array(coords, float), seq, ids


def dssp8(coords):
    """
    Helix (H) and strand (E) from pydssp, a re-implementation of the Kabsch &
    Sander DSSP hydrogen-bond rules. 3-10 (G) and pi (I) helices are added the
    DSSP way: two consecutive 3-turns (or 5-turns) where nothing else is assigned.
    """
    ss = ["C" if c == "-" else c for c in pydssp.assign(coords, out_type="c3")]
    hb = np.asarray(pydssp.get_hbond_map(coords)) > 0.5
    n = len(ss)
    # pydssp's map is [acceptor-side?]: pick the orientation that helices use.
    fwd = sum(hb[i, i + 4] for i in range(n - 4))
    bwd = sum(hb[i + 4, i] for i in range(n - 4))
    bond = (lambda i, j: hb[j, i]) if bwd > fwd else (lambda i, j: hb[i, j])
    for k, lab in ((3, "G"), (5, "I")):
        turn = [i + k < n and bond(i, i + k) for i in range(n)]
        for i in range(1, n):
            span = range(i, min(i + k, n))
            if turn[i - 1] and turn[i] and all(ss[j] not in "HE" for j in span):
                for j in span:
                    ss[j] = lab
    return "".join(ss)


def write_pdb(path, header, chains, hetero=()):
    """chains: [(chain_id, seq3names, ids, coords[L,4,3])]; hetero: [(name, chain, resid, xyz)]"""
    lines = [f"REMARK   1 {h}" for h in header]
    serial = 1
    for ch, names, ids, xyz in chains:
        for name, rid, atoms in zip(names, ids, xyz):
            for an, (x, y, z) in zip(("N", "CA", "C", "O"), atoms):
                el = an[0]
                lines.append(f"ATOM  {serial:5d}  {an:<3s} {name:3s} {ch}{rid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}")
                serial += 1
        lines.append(f"TER   {serial:5d}      {names[-1]:3s} {ch}{ids[-1]:4d}")
        serial += 1
    for name, ch, rid, (x, y, z) in hetero:
        lines.append(f"HETATM{serial:5d} FE   {name:3s} {ch}{rid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          FE")
        serial += 1
    lines.append("END")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def three(seq):
    inv = {v: k for k, v in {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
        "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
        "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}.items()}
    return [inv[a] for a in seq]


def excluded_training(queries):
    """Training proteins with >=25% identity over >=50% of an example protein, or a small-GTPase motif pair."""
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score, al.mode = -11, -1, "local"
    rows = [l.split("\t")[0] for l in open(os.path.join(HERE, "cb6133filtered.tsv")) if not l.startswith("#")]
    out = {}
    for i, s in enumerate(rows):
        if re.search(r"G.{4}GK[ST]", s) and re.search(r"D.{2}G[QH]", s):
            out.setdefault(i, []).append("small-GTPase motifs (P-loop + G3)")
        for name, q in queries.items():
            t = s.replace("X", "A")
            if al.score(q, t) < 60:
                continue
            a = al.align(q, t)[0]
            ident = sum(q[qs + k] == t[ss + k] for (qs, qe), (ss, _) in zip(*a.aligned) for k in range(qe - qs))
            alen = sum(qe - qs for (qs, qe), _ in zip(*a.aligned))
            if alen >= 0.5 * len(q) and ident / alen >= 0.25:
                out.setdefault(i, []).append(f"{ident / alen:.0%} identical to {name} over {alen} residues")
    return out


def main():
    os.makedirs(os.path.join(HERE, "structures"), exist_ok=True)
    ss_rows = []

    ras = pdbio.PDBFile.read(io.StringIO(fetch(URL_1CRR).decode())).get_structure(model=1)
    c, seq, ids = backbone(ras, "A")
    write_pdb(os.path.join(HERE, "structures", "hras_1crr.pdb"),
              ["H-Ras 1-166 (wild type, GDP-bound), PDB 1CRR, NMR model 1 of 20 (Kraulis et al. 1994).",
               "Backbone atoms only. Source file: " + URL_1CRR],
              [("A", three(seq), ids, c)])
    ss_rows.append(("hras", "1CRR", "A", ids[0], seq, dssp8(c)))

    hb = pdbx.get_structure(pdbx.BinaryCIFFile.read(io.BytesIO(gzip.decompress(fetch(URL_4HHB)))), model=1)
    chains = []
    for ch, name in (("A", "hba"), ("B", "hbb"), ("C", None), ("D", None)):
        c, seq, ids = backbone(hb, ch)
        chains.append((ch, three(seq), ids, c))
        if name:
            ss_rows.append((name, "4HHB", ch, ids[0], seq, dssp8(c)))
    fe = hb[(hb.res_name == "HEM") & (hb.atom_name == "FE")]
    write_pdb(os.path.join(HERE, "structures", "hb_4hhb.pdb"),
              ["Human deoxyhaemoglobin, PDB 4HHB (Fermi et al. 1984), chains A,C = alpha, B,D = beta.",
               "Backbone atoms and haem iron atoms only. Source file: " + URL_4HHB],
              chains, [("HEM", a.chain_id, int(a.res_id), a.coord) for a in fe])

    with open(os.path.join(HERE, "reference_ss.tsv"), "w") as fh:
        fh.write("# protein\tpdb\tchain\tfirst_residue\tsequence\tdssp8 (computed by make_reference.py)\n")
        for r in ss_rows:
            fh.write("\t".join(map(str, r)) + "\n")

    queries = {name: seq for name, _, _, _, seq, _ in ss_rows}
    ex = excluded_training(queries)
    with open(os.path.join(HERE, "excluded_training.txt"), "w") as fh:
        fh.write("# 0-based line index in cb6133filtered.tsv (comment line not counted)\treason\n")
        for i in sorted(ex):
            fh.write(f"{i}\t{'; '.join(ex[i])}\n")
    print("wrote structures/, reference_ss.tsv, excluded_training.txt;", len(ex), "training proteins excluded")


if __name__ == "__main__":
    main()
