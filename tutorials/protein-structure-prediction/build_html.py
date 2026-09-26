#!/usr/bin/env python3
"""
Build protein-structure-tutorial.html from protein-structure-tutorial.src.html
by embedding the data files, so the page works offline as a single file.

    python3 build_html.py

Placeholders in the source page:
  __SS_TEST__ / __SS_TRAIN__   data/cb513.tsv, data/cb6133filtered.tsv (sequence<TAB>DSSP)
  __STRUCTURES_JSON__          CA coordinates (and haem irons) from data/structures/*.pdb
  __REFERENCE_SS_JSON__        data/reference_ss.tsv
  __EXCLUDED_JSON__            data/excluded_training.txt
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")


def tsv_body(name):
    with open(os.path.join(DATA, name)) as fh:
        return "".join(l for l in fh if not l.startswith("#")).strip()


def read_structure(name):
    chains, fe = {}, []
    with open(os.path.join(DATA, "structures", name)) as fh:
        for l in fh:
            if l.startswith("ATOM") and l[12:16].strip() == "CA":
                chains.setdefault(l[21], []).append(
                    [int(l[22:26]), l[17:20], round(float(l[30:38]), 2), round(float(l[38:46]), 2), round(float(l[46:54]), 2)])
            elif l.startswith("HETATM") and l[12:16].strip() == "FE":
                fe.append([l[21], round(float(l[30:38]), 2), round(float(l[38:46]), 2), round(float(l[46:54]), 2)])
    return {"chains": chains, "fe": fe}


def main():
    ref = {}
    with open(os.path.join(DATA, "reference_ss.tsv")) as fh:
        for l in fh:
            if l.startswith("#") or not l.strip():
                continue
            name, pdb, chain, first, seq, ss = l.rstrip("\n").split("\t")
            ref[name] = {"pdb": pdb, "chain": chain, "first": int(first), "seq": seq, "dssp": ss}
    with open(os.path.join(DATA, "excluded_training.txt")) as fh:
        excluded = [int(l.split("\t")[0]) for l in fh if l.strip() and not l.startswith("#")]
    structures = {"hras": read_structure("hras_1crr.pdb"), "hb": read_structure("hb_4hhb.pdb")}

    with open(os.path.join(HERE, "protein-structure-tutorial.src.html")) as fh:
        page = fh.read()
    for key, value in {
        "__SS_TEST__": tsv_body("cb513.tsv"),
        "__SS_TRAIN__": tsv_body("cb6133filtered.tsv"),
        "__STRUCTURES_JSON__": json.dumps(structures, separators=(",", ":")),
        "__REFERENCE_SS_JSON__": json.dumps(ref, separators=(",", ":")),
        "__EXCLUDED_JSON__": json.dumps(excluded),
    }.items():
        assert page.count(key) == 1, key
        page = page.replace(key, value)
    out = os.path.join(HERE, "protein-structure-tutorial.html")
    with open(out, "w") as fh:
        fh.write(page)
    print(f"wrote {out} ({len(page) / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
