#!/usr/bin/env python3
"""
Special protein classes tutorial -- companion script.

Sequence-only tools for two classes of protein that generic structure
predictors handle badly: membrane proteins and coiled coils.

  list      show the built-in example proteins (presets) and their sources
  fasta     print a preset (or --seq) as FASTA, e.g. to paste into a web server
  tm        Kyte-Doolittle hydropathy, transmembrane (TM) helix calls and a
            positive-inside topology guess; compares with a reference topology
            when the preset has one
  barrel    why hydropathy misses beta-barrels: hydrophobic moment at 180 deg
            (strict hydrophobic/polar alternation) versus the 19-residue average
  coils     simplified heptad scoring (a teaching method, NOT the published
            COILS matrices): best register per residue, coiled-coil segments,
            and the a/d core used for the oligomer-state rules of thumb
  wheel     text helical wheel and heptad table for a stretch of sequence
  analyze   run tm, barrel and coils together and print a short summary

Only the Python standard library is needed.

Examples:
  python3 special_classes.py list
  python3 special_classes.py tm --preset bR
  python3 special_classes.py tm --preset bR --use-reference
  python3 special_classes.py barrel --preset ompa
  python3 special_classes.py coils --preset gcn4-p1
  python3 special_classes.py wheel --preset gcn4-pLI
  python3 special_classes.py analyze --seq MKTAYIAKQRQISFVKSHFSRQ...
"""

import argparse
import math
import re
import sys

# --------------------------------------------------------------------------
# Example proteins
# --------------------------------------------------------------------------
# Nothing below was typed from memory except the GCN4-p1 leucine zipper.
# Every other sequence was downloaded and cross-checked (identical in at
# least two independent files) -- see the "source" field of each preset and
# the "Where the sequences come from" section of README.md.
#
# Reference topologies are stored run-length encoded, one letter per residue:
#   set160 (TMHMM benchmark set; experimentally determined topologies):
#       i inside (cytoplasm)   o outside   M membrane helix
#   DeepTMHMM training set (topologies derived from 3D structures):
#       I inside   O outside   M membrane helix   B membrane beta-strand
#       P periplasm (the inside of an outer-membrane protein)   S signal peptide

PRESETS = {
    "bR": dict(
        name="Bacteriorhodopsin, Halobacterium salinarum (precursor; 7 TM helices)", acc="P02945",
        source="set160.labels entry BACR_HALHA (old name of this entry); identical to UniProt FASTA in github.com/jomimc/AF2_Stability_PRL_2024 fasta/P02945.fasta",
        seq=(
            "MLELLPTAVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITT"
            "LVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGT"
            "ILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEV"
            "ASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSR"
            "AIFGEAEAPEPSAGDGAAATSD"),
        refs=[
            dict(name="set160 (TMHMM benchmark, experimental topology; github.com/tjs23/python_ml_course set160.labels)",
                 rle=(
                    "o22 M20 i14 M20 o18 M20 i6 M20 o7 M20 i23 M20 o6 M20 i26")),
        ],
    ),
    "gpa": dict(
        name="Glycophorin A, human (precursor with signal peptide; 1 TM helix)", acc="P02724",
        source="set160.labels entry GLPA_HUMAN; residues 20-150 identical to github.com/NENUBioCompute/DMCTOP Dataset.fasta",
        seq=(
            "MYGKIIFVLLLSAIVSISASSTTGVAMHTSTSSSVTKSYISSQTNDTHKRDTYAATPRAH"
            "EVSEISVRTVYPPEEETGERVQLAHHFSEPEITLIIFGVMAGVIGTILLISYGIRRLIKK"
            "SPSDVKPLPSPDTDVPLSSVEIENPETSDQ"),
        refs=[
            dict(name="set160 (TMHMM benchmark, experimental topology; github.com/tjs23/python_ml_course set160.labels)",
                 rle=(
                    "o91 M23 i36")),
        ],
    ),
    "rho": dict(
        name="Rhodopsin, bovine (GPCR; 7 TM helices)", acc="P02699",
        source="set160.labels entry OPSD_BOVIN; identical to github.com/NENUBioCompute/DMCTOP Dataset.fasta",
        seq=(
            "MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLY"
            "VTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLG"
            "GEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIP"
            "EGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQES"
            "ATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAV"
            "YNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA"),
        refs=[
            dict(name="set160 (TMHMM benchmark, experimental topology; github.com/tjs23/python_ml_course set160.labels)",
                 rle=(
                    "o36 M25 i12 M26 o14 M20 i19 M24 o26 M28 i22 M24 o8 M25 i39")),
        ],
    ),
    "lacy": dict(
        name="Lactose permease LacY, E. coli (12 TM helices)", acc="P02920",
        source="set160.labels entry LACY_ECOLI; identical to the P02920 entry of DeepTMHMM.3line",
        seq=(
            "MYYLKNTNFWMFGLFFFFYFFIMGAYFPFFPIWLHDINHISKSDTGIIFAAISLFSLLFQ"
            "PLFGLLSDKLGLRKYLLWIITGMLVMFAPFFIFIFGPLLQYNILVGSIVGGIYLGFCFNA"
            "GAPAVEAFIEKVSRRSNFEFGRARMFGCVGWALCASIVGIMFTINNQFVFWLGSGCALIL"
            "AVLLFFAKTDAPSSATVANAVGANHSAFSLKLALELFRQPKLWFLSLYVIGVSCTYDVFD"
            "QQFANFFTSFFATGEQGTRVFGYVTTMGELLNASIMFFAPLIINRIGGKNALLLAGTIMS"
            "VRIIGSSFATSALEVVILKTLHMFEVPFLLVGCFKYITSQFEVRFSATIYLVCFCFFKQL"
            "AMIFMSVLAGNMYESIGFQGAYLVLGLVALGFTLISVFTLSGPGPLSLLRRQVNEVA"),
        refs=[
            dict(name="set160 (TMHMM benchmark, experimental topology; github.com/tjs23/python_ml_course set160.labels)",
                 rle=(
                    "i10 M23 o13 M21 i7 M25 o3 M23 i19 M19 o4 M20 i24 M23 o25 M22 i9 "
                    "M20 o4 M20 i12 M20 o13 M20 i18")),
            dict(name="DeepTMHMM training set (structure-derived; github.com/andregodinhodtu/TMHMM data/raw_data/DeepTMHMM.3line)",
                 rle=(
                    "I7 M25 O13 M22 I6 M18 O14 M25 I10 M21 O7 M20 I34 M23 O14 M20 I9 "
                    "M21 O4 M21 I15 M21 O11 M19 I17")),
        ],
    ),
    "aqp1": dict(
        name="Aquaporin-1, human (6 TM helices + 2 half-helices)", acc="P29972",
        source="identical in github.com/SBRG/Recon3D (P29972-1.fasta), github.com/korcsmarosgroup/HMIpipeline (P29972.fasta) and DMCTOP Dataset.fasta; no reference topology shipped",
        seq=(
            "MASEFKKKLFWRAVVAEFLATTLFVFISIGSALGFKYPVGNNQTAVQDNVKVSLAFGLSI"
            "ATLAQSVGHISGAHLNPAVTLGLLLSCQISIFRALMYIIAQCVGAIVATAILSGITSSLT"
            "GNSLGRNDLADGVNSGQGLGIEIIGTLQLVLCVLATTDRRRRDLGGSAPLAIGLSVALGH"
            "LLAIDYTGCGINPARSFGSAVITHNFSNHWIFWVGPFIGGALAVLIYDFILAPRSSDLTD"
            "RVKVWTSGQVEEYDLDADDINSRVEMKPK"),
    ),
    "b2ar": dict(
        name="Beta-2 adrenergic receptor, human (GPCR; 7 TM helices)", acc="P07550",
        source="identical in github.com/pablogainza/gprots (uniref/ADRB2/P07550.fasta), github.com/andrewcboardman/pyGEMME and ProteinGym (ai4protein/Pro-Prime); no reference topology shipped",
        seq=(
            "MGQPGNGSAFLLAPNGSHAPDHDVTQERDEVWVVGMGIVMSLIVLAIVFGNVLVITAIAK"
            "FERLQTVTNYFITSLACADLVMGLAVVPFGAAHILMKMWTFGNFWCEFWTSIDVLCVTAS"
            "IETLCVIAVDRYFAITSPFKYQSLLTKNKARVIILMVWIVSGLTSFLPIQMHWYRATHQE"
            "AINCYANETCCDFFTNQAYAIASSIVSFYVPLVIMVFVYSRVFQEAKRQLQKIDKSEGRF"
            "HVQNLSQVEQDGRTGHGLRRSSKFCLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHVIQD"
            "NLIRKEVYILLNWIGYVNSGFNPLIYCRSPDFRIAFQELLCLRRSSLKAYGNGYSSNGNT"
            "GEQSGYHVEQEKENKLLCEDLPGTEDFVGHQGTVPSDNIDSQGRNCSTNDSLL"),
    ),
    "ompa": dict(
        name="Outer membrane protein A, E. coli (precursor; 8-strand beta-barrel + periplasmic domain)", acc="P0A910",
        source="DeepTMHMM.3line entry P0A910; identical to UniProt FASTA in github.com/javieriserte/hack-a-ton-dp, AlessioDelConte/pytorch_disprot and YaoYinYing/FoldDock",
        seq=(
            "MKKTAIAIAVALAGFATVAQAAPKDNTWYTGAKLGWSQYHDTGFINNNGPTHENQLGAGA"
            "FGGYQVNPYVGFEMGYDWLGRMPYKGSVENGAYKAQGVQLTAKLGYPITDDLDIYTRLGG"
            "MVWRADTKSNVYGKNHDTGVSPVFAGGVEYAITPEIATRLEYQWTNNIGDAHTIGTRPDN"
            "GMLSLGVSYRFGQGEAAPVVAPAPAPAPEVQTKHFTLKSDVLFNFNKATLKPEGQAALDQ"
            "LYSQLSNLDPKDGSVVVLGYTDRIGSDAYNQGLSERRAQSVVDYLISKGIPADKISARGM"
            "GESNPVTGNTCDNVKQRAALIDCLAPDRRVEIEVKGIKDVVTQPQA"),
        refs=[
            dict(name="DeepTMHMM training set (structure-derived; github.com/andregodinhodtu/TMHMM data/raw_data/DeepTMHMM.3line)",
                 rle=(
                    "S21 P6 B8 O22 B8 P4 B7 O22 B8 P6 B8 O22 B8 P5 B8 O19 B8 P156")),
        ],
    ),
    "tpm1": dict(
        name="Tropomyosin alpha-1, human (continuous two-stranded coiled coil)", acc="P09493",
        source="identical UniProt FASTA in github.com/Thivby/Master_dissertation, dmx2/myocarditis and pixelatedbus/klasifikasi-protein",
        seq=(
            "MDAIKKKMQMLKLDKENALDRAEQAEADKKAAEDRSKQLEDELVSLQKKLKGTEDELDKY"
            "SEALKDAQEKLELAEKKATDAEADVASLNRRIQLVEEELDRAQERLATALQKLEEAEKAA"
            "DESERGMKVIESRAQKDEEKMEIQEIQLKEAKHIAEDADRKYEEVARKLVIIESDLERAE"
            "ERAELSEGKCAELEEELKTVTNNLKSLEAQAEKYSQKEDRYEEEIKVLSDKLKEAETRAE"
            "FAERSVTKLEKSIDDLEDELYAQKLKYKAISEELDHALNDMTSI"),
    ),
    "hras": dict(
        name="H-Ras, human (soluble negative control)", acc="P01112",
        source="translated from the four H-Ras coding exons in ../../triplets.txt (nt 6229-6339, 6607-6785, 6939-7098, 7796-7915)",
        seq=(
            "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAG"
            "QEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDL"
            "AARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPG"
            "CMSCKCVLS"),
    ),
}

# GCN4 leucine zipper and the Harbury et al. (1993) core variants.
# GCN4-p1 (PDB 2ZTA) in its heptad register (checked with `coils` below):
#
#   R M K Q L E D K V E E L L S K N Y H L E N E V A R L K K L V G E R
#   g a b c d e f g a b c d e f g a b c d e f g a b c d e f g a b c d
#
#   a positions: M2 V9 N16 V23 V30      d positions: L5 L12 L19 L26 R33
#
# The variants replace the core residues a = 9, 16, 23, 30 and d = 5, 12, 19, 26
# (M2 and R33 are unchanged). They are named p + residue at a + residue at d.
GCN4_P1 = "RMKQLEDKVEELLSKNYHLENEVARLKKLVGER"
GCN4_A = (9, 16, 23, 30)
GCN4_D = (5, 12, 19, 26)


def gcn4_variant(a_res, d_res):
    s = list(GCN4_P1)
    for i in GCN4_A:
        s[i - 1] = a_res
    for i in GCN4_D:
        s[i - 1] = d_res
    return "".join(s)


PRESETS.update({
    "gcn4-p1": dict(
        name="GCN4-p1 leucine zipper (wild type, dimer, PDB 2ZTA)", acc="P03069 res. 249-281",
        seq=GCN4_P1,
        source="Well-established 33-residue peptide; identical to the 2ZTA sequence in "
               "github.com/neeleshsoni21/COCONUT coconut/example/pipeline1/2zta/2zta.fasta"),
    "gcn4-pIL": dict(
        name="GCN4-pIL (Ile at a, Leu at d; dimer)", acc="variant",
        seq=gcn4_variant("I", "L"),
        source="Derived from GCN4-p1 by the substitution rule above (not independently checked)"),
    "gcn4-pII": dict(
        name="GCN4-pII (Ile at a, Ile at d; trimer, PDB 1GCM)", acc="variant",
        seq=gcn4_variant("I", "I"),
        source="Derived by substitution; identical to github.com/Indicator/RaptorX-SS8 "
               "examples/1gcm.seq"),
    "gcn4-pLI": dict(
        name="GCN4-pLI (Leu at a, Ile at d; tetramer, PDB 1GCL)", acc="variant",
        seq=gcn4_variant("L", "I"),
        source="Derived by substitution; residues 1-31 identical to SCOP/ASTRAL 2.06 domain "
               "d1gcla_ (github.com/wiwie/clustevalDockerRepository, astral_class_h)"),
})


def decode_rle(rle):
    return "".join(ch * int(n) for ch, n in re.findall(r"([A-Za-z])(\d+)", rle))


# --------------------------------------------------------------------------
# Shared helpers
# --------------------------------------------------------------------------

# Kyte & Doolittle (1982) hydropathy scale (same as structure_tutorial.py).
KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}


def wrap(s, width=60):
    return [s[i:i + width] for i in range(0, len(s), width)]


def sliding_mean(seq, scale, window):
    half = window // 2
    vals = [scale.get(a, 0.0) for a in seq]
    out = [None] * len(seq)
    for i in range(half, len(seq) - half):
        out[i] = sum(vals[i - half:i + half + 1]) / window
    return out


def runs(labels, chars):
    """1-based (start, end) of runs of any character in `chars`."""
    return [(m.start() + 1, m.end()) for m in re.finditer("[" + chars + "]+", labels)]


def print_tracks(rows, width=60):
    """rows = [(label, string)], all the same length; printed in blocks."""
    n = len(rows[0][1])
    for start in range(0, n, width):
        for k, (lab, s) in enumerate(rows):
            pos = f"{start + 1:>5}" if k == 0 else "     "
            print(f"{pos} {lab:<4} {s[start:start + width]}")
        print()


# --------------------------------------------------------------------------
# Part A1: transmembrane helices and the positive-inside rule
# --------------------------------------------------------------------------

def tm_segments(seq, window=19, cutoff=1.6):
    """
    Greedy peak picking on the windowed Kyte-Doolittle profile:
    take the highest window above `cutoff`, call its `window` residues a TM
    helix, then repeat with the next-highest window that does not overlap a
    helix already called. Returns [(start, end, peak_value)], 1-based.
    """
    hyd = sliding_mean(seq, KD, window)
    half = window // 2
    taken = [False] * len(seq)
    centres = sorted((i for i, v in enumerate(hyd) if v is not None and v > cutoff),
                     key=lambda i: -hyd[i])
    segs = []
    for c in centres:
        a, b = c - half, c + half
        if any(taken[a:b + 1]):
            continue
        for k in range(a, b + 1):
            taken[k] = True
        segs.append((a + 1, b + 1, hyd[c]))
    return sorted(segs), hyd


def loops_between(n, segs):
    out, prev = [], 0
    for a, b in segs:
        out.append((prev + 1, a - 1))
        prev = b
    out.append((prev + 1, n))
    return out


def positive_inside(seq, segs, flank=15):
    """
    von Heijne's positive-inside rule, simplified: count K + R within `flank`
    residues of each end of every TM segment (termini: only the end next to
    the membrane). Loops alternate sides, so the N-terminus is on the same side
    as loops 0, 2, 4, ... The side with more K + R is predicted to be inside.
    """
    loops = loops_between(len(seq), segs)
    counts = []
    for j, (a, b) in enumerate(loops):
        idx = set()
        if j > 0:                       # residues just after the previous helix
            idx |= set(range(a, min(b, a + flank - 1) + 1))
        if j < len(loops) - 1:          # residues just before the next helix
            idx |= set(range(max(a, b - flank + 1), b + 1))
        counts.append(sum(1 for i in idx if seq[i - 1] in "KR"))
    n_side = sum(counts[0::2])
    other = sum(counts[1::2])
    guess = "in" if n_side > other else "out" if other > n_side else "undecided"
    return loops, counts, n_side, other, guess


def match_segments(pred, ref, min_overlap=5):
    """A predicted segment matches a reference one if they share >= 5 residues."""
    used, hits = set(), 0
    for a, b in pred:
        for k, (c, d) in enumerate(ref):
            if k not in used and min(b, d) - max(a, c) + 1 >= min_overlap:
                used.add(k)
                hits += 1
                break
    return hits, len(pred) - hits, len(ref) - hits


def topology_string(n, segs, n_term):
    """Per-residue i/o/M string from TM segments and the N-terminal side."""
    out = []
    side = "i" if n_term == "in" else "o"
    prev = 0
    for a, b in segs:
        out.append(side * (a - 1 - prev))
        out.append("M" * (b - a + 1))
        side = "o" if side == "i" else "i"
        prev = b
    out.append(side * (n - prev))
    return "".join(out)


def ref_n_side(labels):
    """'in'/'out' for the first non-signal, non-membrane residue of a reference."""
    for ch in labels:
        if ch in "iIP":
            return "in"
        if ch in "oO":
            return "out"
    return "?"


def cmd_tm(seq, label, preset, window, cutoff, flank, use_reference):
    segs3, hyd = tm_segments(seq, window, cutoff)
    pred = [(a, b) for a, b, _ in segs3]
    refs = (preset or {}).get("refs", [])

    print(f"== {label} ==")
    print(f"Length {len(seq)} aa")
    vals = [v for v in hyd if v is not None]
    if vals:
        peak = max(vals)
        print(f"Kyte-Doolittle, {window}-residue window: max {peak:+.2f} at residue "
              f"{hyd.index(peak) + 1}; cutoff {cutoff}")
    print()

    if use_reference:
        if not refs or not runs(refs[0]["labels"], "M"):
            sys.exit("This preset has no reference helices; drop --use-reference.")
        segs = runs(refs[0]["labels"], "M")
        print(f"Using the REFERENCE helices from {refs[0]['name']} "
              "(tests the positive-inside rule on its own).")
    else:
        segs = pred
        if not segs:
            print(f"No TM helix called (no {window}-residue window with KD > {cutoff}).")
        else:
            print(f"Predicted TM helices ({len(segs)}):")
            print("   #   start   end   peak KD")
            for k, (a, b, p) in enumerate(segs3, 1):
                print(f"  {k:>2}   {a:>5} {b:>5}   {p:+.2f}")
            if segs[0][0] <= 20:
                print(f"  ! helix 1 starts at residue {segs[0][0]}: it could be a cleaved "
                      "signal peptide, not a TM helix. Hydropathy cannot tell; use Phobius, "
                      "SignalP 6.0 or DeepTMHMM.")
    print()

    if segs:
        loops, counts, n_side, other, guess = positive_inside(seq, segs, flank)
        print(f"Positive-inside rule: K + R within {flank} residues of each helix end")
        row1 = "  loop     " + " ".join(f"{('N' if j == 0 else 'C' if j == len(loops) - 1 else str(j)):>4}"
                                        for j in range(len(loops)))
        row2 = "  side     " + " ".join(f"{('A' if j % 2 == 0 else 'B'):>4}" for j in range(len(loops)))
        row3 = "  K+R      " + " ".join(f"{c:>4}" for c in counts)
        print(row1)
        print(row2)
        print(row3)
        print(f"  side A (N-terminal side) {n_side}   side B {other}"
              f"   ->  N-terminus predicted {guess.upper()}"
              + ("" if guess == "undecided" else
                 f", C-terminus {'in' if (guess == 'in') == (len(segs) % 2 == 0) else 'out'}"))
        print()

    for ref in refs:
        rseg = runs(ref["labels"], "M")
        rb = runs(ref["labels"], "B")
        print(f"Reference: {ref['name']}")
        if rseg:
            print(f"  {len(rseg)} TM {'helix' if len(rseg) == 1 else 'helices'}: "
                  + ", ".join(f"{a}-{b}" for a, b in rseg))
            print(f"  N-terminus {ref_n_side(ref['labels'])}")
            if not use_reference:
                hit, fp, fn = match_segments(pred, rseg)
                print(f"  prediction vs reference: {hit} matched, {fp} extra, {fn} missed "
                      "(match = overlap of 5 or more residues)")
        if rb:
            print(f"  {len(rb)} membrane beta-strands (no TM helices): "
                  + ", ".join(f"{a}-{b}" for a, b in rb))
        sp = runs(ref["labels"], "S")
        if sp:
            print(f"  signal peptide {sp[0][0]}-{sp[0][1]}")
        print()

    if segs:
        g = positive_inside(seq, segs, flank)[4]
        topo = topology_string(len(seq), segs, g if g != "undecided" else "in")
        rows = [("seq", seq), ("pred", topo)]
        for ref in refs:
            rows.append(("ref", ref["labels"]))
        print("Topology track: i inside, o outside, M membrane"
              + (" (reference letters as in its source)" if refs else ""))
        print_tracks(rows)


# --------------------------------------------------------------------------
# Part A2: beta-barrels and the hydrophobic moment
# --------------------------------------------------------------------------

def hydrophobic_moment(seq, window, degrees):
    """
    Eisenberg-style hydrophobic moment with the Kyte-Doolittle scale:
    |sum h_i * exp(i * delta * k)| / window. delta = 180 deg detects strict
    alternation (a beta-strand with one face in lipid), delta = 100 deg an
    amphipathic alpha-helix.
    """
    half = window // 2
    d = math.radians(degrees)
    out = [None] * len(seq)
    for i in range(half, len(seq) - half):
        c = s = 0.0
        for k in range(i - half, i + half + 1):
            h = KD.get(seq[k], 0.0)
            c += h * math.cos(d * k)
            s += h * math.sin(d * k)
        out[i] = math.hypot(c, s) / window
    return out


def cmd_barrel(seq, label, preset, window, threshold):
    mu = hydrophobic_moment(seq, window, 180)
    hyd = sliding_mean(seq, KD, 19)
    segs = tm_segments(seq)[0]
    print(f"== {label} ==")
    vals = [v for v in hyd if v is not None]
    print(f"TM helices called by hydropathy (19-window > 1.6): {len(segs)}"
          + (" at " + ", ".join(f"{a}-{b}" for a, b, _ in segs) if segs else ""))
    if vals:
        print(f"Max 19-window hydropathy: {max(vals):+.2f}")
    print(f"Hydrophobic moment at 180 deg, {window}-residue window "
          f"(high = hydrophobic and polar residues alternate)")
    track = "".join("." if v is None else "#" if v >= threshold else "-" for v in mu)
    mv = [v for v in mu if v is not None]
    print(f"  residues with moment >= {threshold}: {sum(1 for v in mv if v >= threshold)} of {len(mv)}")
    rows = [("seq", seq), ("mu", track)]
    refs = [r for r in (preset or {}).get("refs", []) if "B" in r["labels"]]
    if refs:
        lab = refs[0]["labels"]
        rows.append(("ref", lab))
        inb = [v for v, c in zip(mu, lab) if v is not None and c == "B"]
        outb = [v for v, c in zip(mu, lab) if v is not None and c != "B"]
        strands = runs(lab, "B")
        found = sum(1 for a, b in strands
                    if sum(1 for i in range(a - 1, b) if mu[i] is not None and mu[i] >= threshold) >= 3)
        print(f"  mean moment in reference strands  {sum(inb) / len(inb):.2f}")
        print(f"  mean moment everywhere else       {sum(outb) / len(outb):.2f}")
        print(f"  reference strands with 3 or more residues above {threshold}: "
              f"{found} of {len(strands)}")
        hs = [v for v, c in zip(hyd, lab) if v is not None and c == "B"]
        print(f"  mean 19-window hydropathy over strand residues {sum(hs) / len(hs):+.2f} "
              "(far below the 1.6 helix cutoff)")
    print(f"  '#' moment >= {threshold}   '-' below   '.' window does not fit")
    print()
    print_tracks(rows)


# --------------------------------------------------------------------------
# Part B: coiled coils
# --------------------------------------------------------------------------

HEPTAD = "abcdefg"
CORE = set("LIVMFY")          # hydrophobic residues that pack well at a and d
CHARGED = set("EKRD")         # salt-bridge formers at e and g
SMALL_POLAR = set("QNSTHA")   # acceptable at e and g
SURFACE = set("EKRDQNSTHAG")  # acceptable on the exposed b, c, f positions


def heptad_score(aa, pos):
    """
    Simplified, transparent scoring (a teaching method, not COILS/MTIDK):
      a, d : 2 if L I V M F Y, else 0
      e, g : 1 if E K R D, 0.5 if Q N S T H A, else 0
      b c f: 1 if polar/charged/small (E K R D Q N S T H A G), else 0
      Pro anywhere: -1 (breaks the helix)
    A perfect heptad scores 2+2+1+1+1+1+1 = 9.
    """
    if aa == "P":
        return -1.0
    if pos in "ad":
        return 2.0 if aa in CORE else 0.0
    if pos in "eg":
        return 1.0 if aa in CHARGED else 0.5 if aa in SMALL_POLAR else 0.0
    return 1.0 if aa in SURFACE else 0.0


def coils_scan(seq, window=28):
    """
    Slide a `window`-residue frame along the sequence; score each frame in all
    7 registers (normalised to 0-1, 1 = every heptad perfect). Every residue
    keeps the best score of any frame and register that covers it, and the
    heptad letter it has in that frame.
    """
    n = len(seq)
    w = min(window, n)
    best = [(0.0, "-")] * n
    for s in range(n - w + 1):
        for r in range(7):
            letters = [HEPTAD[(r + k) % 7] for k in range(w)]
            score = sum(heptad_score(seq[s + k], letters[k]) for k in range(w)) / (9.0 * w / 7.0)
            for k in range(w):
                if score > best[s + k][0]:
                    best[s + k] = (score, letters[k])
    return best


def oligomer_rule(seq, register, segs):
    """Harbury et al. (1993) rules of thumb from the a and d residues."""
    a = [seq[i] for s, e in segs for i in range(s - 1, e) if register[i] == "a"]
    d = [seq[i] for s, e in segs for i in range(s - 1, e) if register[i] == "d"]
    if not a or not d:
        return a, d, "no coiled-coil segment"
    beta_a = sum(1 for x in a if x in "IV") / len(a)
    leu_a = sum(1 for x in a if x == "L") / len(a)
    beta_d = sum(1 for x in d if x in "IV") / len(d)
    leu_d = sum(1 for x in d if x == "L") / len(d)
    if beta_a >= 0.5 and leu_d >= 0.5:
        verdict = "DIMER-like (beta-branched a, Leu d: like GCN4-p1 / pIL)"
    elif beta_a >= 0.5 and beta_d >= 0.5:
        verdict = "TRIMER-like (Ile at a and d: like GCN4-pII)"
    elif leu_a >= 0.5 and beta_d >= 0.5:
        verdict = "TETRAMER-like (Leu a, Ile d: like GCN4-pLI)"
    else:
        verdict = "no clear rule applies"
    if "N" in a:
        verdict += "; Asn at an a position favours a dimer"
    return a, d, verdict


def cmd_coils(seq, label, window, threshold):
    best = coils_scan(seq, window)
    scores = [s for s, _ in best]
    reg = "".join(r for _, r in best)
    segs = [(a, b) for a, b in runs("".join("C" if s >= threshold else "." for s in scores), "C")
            if b - a + 1 >= 14]
    print(f"== {label} ==")
    print(f"Simplified heptad scoring, {window}-residue window (teaching method, not COILS)")
    print(f"Score range {min(scores):.2f} to {max(scores):.2f}; "
          f"{sum(1 for s in scores if s >= threshold)} of {len(seq)} residues >= {threshold}")
    if segs:
        print("Coiled-coil segments (14 or more residues >= threshold): "
              + ", ".join(f"{a}-{b}" for a, b in segs))
    else:
        print("No coiled-coil segment called.")
    print()
    digits = "".join(str(min(9, int(s * 10))) if s > 0 else "0" for s in scores)
    shown = "".join(r if s >= threshold else "." for s, r in best)
    print("score: first decimal of the best score (8 = 0.80-0.89)   reg: heptad letter where score >= "
          f"{threshold}")
    print_tracks([("seq", seq), ("scr", digits), ("reg", shown)])
    a, d, verdict = oligomer_rule(seq, reg, segs)
    if segs:
        print(f"a positions: {''.join(a)}")
        print(f"d positions: {''.join(d)}")
        print(f"Oligomer rule of thumb: {verdict}")
        print("  (Only for short, GCN4-like zippers. Real oligomer state also depends on e/g "
              "residues, length and context: test 2, 3 and 4 copies in AlphaFold-Multimer.)")


def cmd_wheel(seq, label, start, end, register):
    best = coils_scan(seq)
    if end is None:
        end = min(len(seq), start + 17)
    if register is None:
        register = best[start - 1][1]
        if register == "-":
            register = "a"
        how = "best register from `coils`"
    else:
        how = "register chosen with --register"
    offset = HEPTAD.index(register)
    sub = [(i, seq[i - 1], HEPTAD[(offset + i - start) % 7]) for i in range(start, end + 1)]
    print(f"== {label}: residues {start}-{end}, residue {start} = '{register}' ({how}) ==")
    print()
    print("Helical wheel (100 degrees per residue, 3.6 residues per turn).")
    print("Spokes in the order they appear around the wheel; innermost residue first.")
    print()
    angle = {p: (HEPTAD.index(p) * 100) % 360 for p in HEPTAD}
    for p in sorted(HEPTAD, key=lambda p: angle[p]):
        res = [f"{aa}{i}" for i, aa, q in sub if q == p]
        tag = "core" if p in "ad" else "edge" if p in "eg" else "surface"
        print(f"  {p} {angle[p]:>3} deg  {tag:<7}  " + " ".join(res))
    print()
    print("Heptad table")
    print("   " + "  ".join(HEPTAD))
    row, rows = [" "] * 7, []
    for i, aa, p in sub:
        k = HEPTAD.index(p)
        if k == 0 and any(c != " " for c in row):
            rows.append(row)
            row = [" "] * 7
        row[k] = aa
    rows.append(row)
    for r in rows:
        print("   " + "  ".join(r))
    core = [aa for _, aa, p in sub if p in "ad"]
    hyd = sum(1 for x in core if x in CORE)
    print()
    print(f"Hydrophobic (L I V M F Y) at a/d: {hyd} of {len(core)}")
    ch = [aa for _, aa, p in sub if p in "eg"]
    print(f"Charged (E K R D) at e/g: {sum(1 for x in ch if x in CHARGED)} of {len(ch)}")


# --------------------------------------------------------------------------
# analyze: everything at once
# --------------------------------------------------------------------------

def cmd_analyze(seq, label, preset):
    segs3, hyd = tm_segments(seq)
    segs = [(a, b) for a, b, _ in segs3]
    vals = [v for v in hyd if v is not None]
    n = len(seq)
    print(f"== {label} ==")
    print(f"Length {n} aa   GRAVY {sum(KD.get(a, 0) for a in seq) / n:+.2f}")
    if vals:
        print(f"Max 19-window hydropathy {max(vals):+.2f}")
    print(f"TM helices (KD 19 > 1.6): {len(segs)}"
          + (" at " + ", ".join(f"{a}-{b}" for a, b in segs) if segs else ""))
    if segs:
        g = positive_inside(seq, segs)
        print(f"Positive-inside: K+R side A {g[2]}, side B {g[3]} -> N-terminus {g[4]}")
        if segs[0][0] <= 20:
            print("  first helix is near the N-terminus: possible signal peptide")
    mu = [v for v in hydrophobic_moment(seq, 9, 180) if v is not None]
    if mu:
        print(f"Residues with 180-deg hydrophobic moment >= 1.2 (9-window): "
              f"{sum(1 for v in mu if v >= 1.2)} of {len(mu)}")
    best = coils_scan(seq)
    scores = [s for s, _ in best]
    reg = "".join(r for _, r in best)
    cc = [(a, b) for a, b in runs("".join("C" if s >= 0.75 else "." for s in scores), "C")
          if b - a + 1 >= 14]
    print(f"Heptad score max {max(scores):.2f}; coiled-coil segments (>= 0.75): "
          + (", ".join(f"{a}-{b}" for a, b in cc) if cc else "none"))
    if cc:
        print(f"Oligomer rule of thumb: {oligomer_rule(seq, reg, cc)[2]}")
    for ref in (preset or {}).get("refs", []):
        m, b = runs(ref["labels"], "M"), runs(ref["labels"], "B")
        print(f"Reference ({ref['name'].split(' (')[0]}): {len(m)} TM helices, {len(b)} beta-strands, "
              f"N-terminus {ref_n_side(ref['labels'])}")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def get_seq(args):
    if args.seq:
        seq = re.sub(r"[^A-Za-z]", "", args.seq).upper()
        return seq, "user sequence", None
    key = args.preset
    if key not in PRESETS:
        sys.exit(f"Unknown preset {key!r}. Try: {', '.join(PRESETS)}")
    p = dict(PRESETS[key])
    p["refs"] = [dict(r, labels=decode_rle(r["rle"])) for r in p.get("refs", [])]
    for r in p["refs"]:
        assert len(r["labels"]) == len(p["seq"]), key
    return p["seq"], f"{key}: {p['name']}", p


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list", help="list the presets and their sources")

    def add(name, help_):
        p = sub.add_parser(name, help=help_)
        g = p.add_mutually_exclusive_group()
        g.add_argument("--preset", default="bR", help="built-in example (see `list`)")
        g.add_argument("--seq", help="analyse this protein sequence instead")
        return p

    add("fasta", "print the sequence as FASTA")
    p = add("tm", "hydropathy, TM helices and positive-inside topology")
    p.add_argument("--window", type=int, default=19)
    p.add_argument("--cutoff", type=float, default=1.6)
    p.add_argument("--flank", type=int, default=15)
    p.add_argument("--use-reference", action="store_true",
                   help="apply the positive-inside rule to the reference helices")
    p = add("barrel", "hydrophobic moment at 180 degrees (beta-barrels)")
    p.add_argument("--window", type=int, default=9)
    p.add_argument("--threshold", type=float, default=1.2)
    p = add("coils", "simplified heptad scoring")
    p.add_argument("--window", type=int, default=28)
    p.add_argument("--threshold", type=float, default=0.75)
    p = add("wheel", "text helical wheel and heptad table")
    p.add_argument("--start", type=int, default=1)
    p.add_argument("--end", type=int)
    p.add_argument("--register", choices=list(HEPTAD),
                   help="heptad letter of the first residue (default: best from coils)")
    add("analyze", "summary of all analyses")
    args = ap.parse_args()

    if args.cmd == "list":
        for k, p in PRESETS.items():
            print(f"{k:<9} {len(p['seq']):>4} aa  {p['name']}")
            print(f"{'':<16}{p['acc']}; {p['source']}")
            for r in p.get("refs", []):
                print(f"{'':<16}reference topology: {r['name']}")
        return

    seq, label, preset = get_seq(args)
    if args.cmd == "fasta":
        print(">" + label.replace(": ", " ", 1))
        print("\n".join(wrap(seq)))
    elif args.cmd == "tm":
        cmd_tm(seq, label, preset, args.window, args.cutoff, args.flank, args.use_reference)
    elif args.cmd == "barrel":
        cmd_barrel(seq, label, preset, args.window, args.threshold)
    elif args.cmd == "coils":
        cmd_coils(seq, label, args.window, args.threshold)
    elif args.cmd == "wheel":
        cmd_wheel(seq, label, args.start, args.end, args.register)
    elif args.cmd == "analyze":
        cmd_analyze(seq, label, preset)


if __name__ == "__main__":
    main()
