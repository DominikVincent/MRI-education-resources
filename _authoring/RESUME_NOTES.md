# Resume notes

First draft is COMPLETE (2026-11-11): all 56 notebooks exist, execute cleanly, QA passes, book builds with 0 warnings.
Open work for the polishing pass: see POLISH_LIST.md.

## How to build
- One notebook (execute + link check):
  `/home/dominik/miniconda3/envs/GestureDetection/bin/python _authoring/build_nb.py _authoring/nbsrc/<folder>/<Name>.nbsrc <folder>/<Name>.ipynb --execute`
  (widget/interact cells are skipped during execution on purpose)
- QA over everything: `/home/dominik/miniconda3/envs/GestureDetection/bin/python _authoring/qa.py`
- Web book (Jupyter Book, env mri-book): `/home/dominik/miniconda3/envs/mri-book/bin/python _authoring/make_book.py`
  → `_build/html/index.html`. After adding/removing chapters run `jb clean .` first (stale pages keep dead links).
- Reference notebooks are generated from templates: `_authoring/reftpl/*.in` → `reftpl/expand.py` → nbsrc (edit the .in files).

## Notes
- Sources of truth: `_authoring/nbsrc/**.nbsrc` (the .ipynb are build outputs).
- AUTHORING_GUIDE.md still mentions the old /tmp scratchpad paths; read SCRATCH as `_authoring/` and SRC as a fresh clone of
  https://github.com/LarsonLab/MRI-education-resources (only needed for original images/text, already copied: images/, data/, txt/).
- Image files with spaces were duplicated without spaces (Birdcage_schematic.png etc.) because Sphinx can't load %20 paths.
