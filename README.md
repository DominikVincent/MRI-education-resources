# MRI — from intuition to engineering

A course on how MRI works, written as Jupyter notebooks: the physics, the mathematics, and the
hardware/engineering. It adapts and extends Peder Larson's open textbook
[*Principles of MRI*](https://larsonlab.github.io/MRI-education-resources) (MIT License) for readers
without a physics background. Every symbol is defined, every equation is explained in words, missing
derivations are added, there are new hardware chapters, and all code is in Python.

**Web version:** https://dominikvincent.github.io/MRI-education-resources/ (built by GitHub Actions from the `mri-course` branch).

**Open [`00_Start_Here.ipynb`](00_Start_Here.ipynb) first.**

## Layout

| Folder | Content |
|---|---|
| `00_Start_Here.ipynb` | Orientation, background self-check, course map, conventions |
| `pass1_big_picture/` | Pass 1: an intuitive tour of all of MRI in 5 notebooks |
| `pass2_deep_dive/` | Pass 2: the 30 original chapters, extended, plus 4 new hardware/safety chapters (B3–B6) |
| `refreshers/` | Optional math (RM1–RM6) and physics (RP1–RP4) background notebooks |
| `reference/` | Notation, glossary & acronyms, math cheat sheet, key equations, practice questions, software |
| `mri_tools.py` | Shared Python helpers (Bloch simulator, FFTs, phantoms, signal equations, coil fields) |
| `images/`, `data/` | Figures, animations and raw data from the original book |

## Running

The notebooks are stored with their outputs, so you can read them without running anything. To run them:

```bash
conda activate GestureDetection        # existing env with numpy/scipy/matplotlib/ipywidgets/jupyter
jupyter lab
```

Or create a fresh environment: `conda env create -f environment.yml && conda activate mri-course`.

## Building the web book locally

```bash
pip install -r requirements-book.txt          # or use the existing conda env "mri-book"
python _authoring/make_book.py                 # regenerates _toc.yml/_config.yml and runs jupyter-book
# -> _build/html/index.html
```

Notebook sources live in `_authoring/nbsrc/` (`.nbsrc` files). `_authoring/build_nb.py` turns them into executed notebooks; see `_authoring/RESUME_NOTES.md`.

## License / attribution

Text, figures, data and code ideas are adapted from P. Larson, *Principles of MRI*,
https://github.com/LarsonLab/MRI-education-resources, © 2021 Peder Larson, MIT License
(see `LICENSE_original_MIT.txt`).
