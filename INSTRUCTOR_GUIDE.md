# Instructor Guide — AI-Assisted Bioinformatics Exercise

This exercise teaches **AI-assisted coding discipline** through bioinformatics.
The deliverable students hand in is a run notebook where every ✅ verify cell
passes. Grade the *verification*, not just the output.

## At a glance
- **Duration:** ~90 min (60 min core Sections 1–7, 30 min Section 8 challenge).
- **Group size:** works solo or in pairs (pair-prompting is great here).
- **Prereqs:** basic Python; DNA/codon/central-dogma concepts; an AI assistant.
- **Files:** student worksheet + pre-run solution in `notebooks/`.

## The one idea
> **Decompose → Prompt → Verify.** AI makes students faster; biology-based
> verification makes them right. An unverified plot is *worse* than no plot — it
> looks authoritative while being wrong.

## Suggested timing
| Time | Section | Focus |
|------|---------|-------|
| 0:00 | 0–1 | Framing + the 5 AI failure modes; import check. |
| 0:10 | 2 | Fetch from NCBI; the fallback pattern; verify the record. |
| 0:25 | 3–4 | GC content + first Matplotlib plot (verify vs library). |
| 0:45 | 5 | Six-frame ORF search — the hardest step; frame bugs. |
| 1:05 | 6–7 | Codon usage + amino-acid composition plots. |
| 1:20 | 8–9 | Open challenge + reflection. |

## Learning objectives → where they're assessed
1. Decompose problems → every 🧭 badge; Section 8 (no scaffold).
2. Prompt well → 💬 badges; reflection Q3 (rewrite a prompt).
3. Verify code → every ✅ cell; three techniques below.
4. Biopython (Entrez/SeqIO/Seq, ORFs, translation) → Sections 2, 5.
5. Matplotlib (GC landscape, codon, AA) → Sections 4, 6, 7.
6. Name the 5 failure modes → Section 0; reflection Q1.

## The three verification techniques to reinforce
- **Against a trusted library:** hand-rolled `gc_content` vs Biopython
  `gc_fraction` (Section 3).
- **Against known biology:** an ORF must start `M`, end `*`, have no internal
  stops; every codon must translate to the expected amino acid (Sections 5–6).
- **On tiny hand-checkable inputs + assertions:** the `ATGAAATAG → MK*` ORF test;
  base counts summing to length (Sections 4–7).

## AI pitfalls to seed on purpose
Encourage students to *let the AI fail* so they see verification catch it. Good
prompts to have them try (and then debug):
- Ask for GC content and many assistants suggest `Bio.SeqUtils.GC(seq)` — it's
  **deprecated/removed** in modern Biopython (→ `gc_fraction`). *Failure mode 5.*
- Ask for an ORF finder and a common bug is scanning only the **forward strand**
  or only **frame 0** — the six-frame requirement exposes it. *Failure mode 3.*
- Ask to "translate the sequence" and an assistant may translate a sequence whose
  length isn't a multiple of 3 (partial codon warning) or forget the strand.
  *Failure modes 2 & 4.*
- Ask for `seq.gc_content()` and you may get a confident **hallucinated method**.
  *Failure mode 1.*

## Expected results (offline fallback sequence)
Running the solution against `data/dna_sequence.fasta` yields:
- Length **10,309 bp**, overall GC **~65.9%**.
- Longest ORF: **330 aa**, **+** strand, translates cleanly (starts `M`, ends
  `*`, no internal stops).
- 331 codons counted (330 residues + stop); reverse strand is stop-riddled.

With internet and the default `NM_000518` (human *HBB*), students instead get a
clean **147-aa β-globin** CDS beginning `MVHLTPEEK…` — a great extra
verification target (compare the AI translation to the textbook protein).

## Grading rubric (100 pts)
| Criterion | Pts |
|-----------|-----|
| Sections 1–7 run top-to-bottom; all ✅ verify cells pass | 40 |
| Code is explained/understood (spot-check by asking a student to walk through one function) | 15 |
| Verifications are meaningful, not just "it ran" (e.g. compared to library/biology) | 20 |
| Section 8 challenge attempted, **with its own stated verification** | 15 |
| Reflection answers the 3 questions specifically (names a real AI error caught) | 10 |

## Common student stumbles
- **`Entrez.email` not set** → NCBI 400/blocked. The fallback still runs, but
  point out why the email matters for real work.
- **Editing the `.ipynb` and losing changes on regen** → tell students to edit
  cells, not `tools/build_notebooks.py`, unless they're maintaining the exercise.
- **"It ran, so it's right."** → the core misconception this whole exercise
  exists to break. Reward provable correctness.

## Adapting the exercise
- Swap `ACCESSION` for an organism relevant to your course.
- Harder version: remove the provided plotting cells in Sections 6–7 and make
  them TODOs too.
- Assessment variant: give a *broken* AI-written function and ask students to
  find the bug using only verification cells.
