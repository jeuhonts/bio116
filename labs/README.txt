BIO 116 · HRAS bioinformatics labs (start page + 4 self-contained HTML files)

  index.html                   Start page: introduces the labs, links to each,
                               shows each student's progress, and can clear
                               saved progress (useful on shared computers)
  1-sequencing-assembly.html   Lab 1 · Sequencing & assembly
  2-gene-finding.html          Lab 2 · Gene finding
  3-genome-browser.html        Lab 3 · Genome browser (includes the graded exercise)
  4-genome-annotation.html     Lab 4 · Genome annotation

Keep all five files in the same folder and point students to index.html.
The "Lab n of 4" bar, the previous/next links and the "BIO 116 labs"
link back to the start page all use these file names.

Each file runs on its own in any modern browser. It loads nothing but
Google Fonts (and falls back to system fonts without them). Nothing a
student types is sent anywhere; practice progress is kept in that
student's own browser (localStorage).

Lab 3 graded exercise: "Submit and save file" downloads a JSON file for
students to upload to the course site. Grade the files with the separate
instructor grader (not included here; do not give it to students). The
answer keys are not in these student files, but the whole exercise runs
in the browser, so treat it as low-stakes: the submission "seal" catches
casual edits, not a determined student.
