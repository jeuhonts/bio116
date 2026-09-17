# BIO116 — Bioinformatics

## AI-Assisted Coding Exercise: Characterizing an Unknown DNA Sequence

A hands-on Jupyter exercise that teaches students **how to code with an AI
assistant** — using real bioinformatics tasks built on **Biopython** and
**Matplotlib** as the vehicle.

The biology (GC content, reading frames, codon usage, translation) is genuinely
useful, but the transferable skill is a disciplined way of working with an AI
coding assistant:

> ## 🧭 Decompose → 💬 Prompt → ✅ Verify
>
> Break the problem into a small step, prompt the assistant for *just that step*,
> then **verify the result against biology you already trust.** An AI assistant
> is a fast, confident colleague who is sometimes wrong — the student's job is to
> stay the expert who checks the work.

---

## What's in here

| Path | What it is |
|------|-----------|
| `notebooks/ai_bioinformatics_exercise.ipynb` | **Student worksheet** — 8 guided TODOs with prompts and verification checks. |
| `notebooks/ai_bioinformatics_solution.ipynb` | **Instructor solution** — fully worked and pre-executed, with plots. |
| `data/dna_sequence.fasta` | Offline fallback sequence (real ~10.3 kb fragment; see `data/README.md`). |
| `INSTRUCTOR_GUIDE.md` | Learning objectives, timing, grading rubric, and the AI pitfalls to seed. |
| `tools/build_notebooks.py` | Regenerates both notebooks from one shared definition. |
| `requirements.txt` | Python dependencies. |
| `triplets.txt` | Original raw data the fallback sequence was reconstructed from. |

---

## Learning objectives

By the end, a student can:

1. **Decompose** a bioinformatics task into steps an assistant can solve.
2. **Write focused prompts** that give context, constraints, and a definition of
   done.
3. **Verify AI-generated code** three ways — against a trusted library, against
   known biology (the genetic code, a translated protein), and with assertions
   on tiny hand-checkable inputs.
4. Use **Biopython** to fetch sequences from NCBI (`Bio.Entrez`/`Bio.SeqIO`),
   manipulate them (`Bio.Seq`), find ORFs, and translate.
5. Use **Matplotlib** to plot a GC landscape, codon usage, and amino-acid
   composition.
6. Recognize the **five common failure modes** of AI code in bioinformatics
   (hallucinated APIs, wrong genetic code, frame/off-by-one errors, type
   confusion, stale library usage).

**Audience:** students comfortable with basic Python and molecular-biology
concepts, new to Biopython and to AI-assisted coding.
**Time:** ~90 minutes. **Slides/prereqs:** none beyond this repo.

---

## Setup

```bash
git clone <this-repo>
cd bio116
python -m venv .venv && source .venv/bin/activate    # optional but recommended
pip install -r requirements.txt
jupyter notebook notebooks/ai_bioinformatics_exercise.ipynb
```

Then have your AI coding assistant open (Claude, Copilot, ChatGPT, or a local
model) and work through the notebook top to bottom.

### Online vs offline
- **With internet:** the notebook fetches a real gene live from NCBI (default:
  human β-globin `HBB`). Set your own email in the config cell — NCBI requires
  it.
- **Without internet (exam room, no wifi, locked-down lab):** the notebook
  automatically falls back to the bundled `data/dna_sequence.fasta`, so *every
  cell still runs.* Building that fallback is itself a lesson in robust code.

---

## How it's taught

Each task cell follows the same rhythm, marked with badges:

- 🧭 **Decompose** — the sub-problem, stated precisely.
- 💬 **Prompt** — a starter prompt to paste into the assistant (students are
  encouraged to improve it).
- ✍️ **Your code** — a `TODO` scaffold to fill using the assistant.
- ✅ **Verify** — a concrete check that must pass before trusting the output.

The verification is the point. A student who pastes AI code and gets a plot has
done nothing; a student who *proves* the plot is right has learned the job.

---

## Regenerating the notebooks

Both notebooks are generated from a single source of truth so they never drift:

```bash
python tools/build_notebooks.py
# validate the solution end to end (offline path):
cd notebooks && jupyter nbconvert --to notebook --execute --inplace ai_bioinformatics_solution.ipynb
```

Edit the content in `tools/build_notebooks.py`, not the `.ipynb` files directly.
