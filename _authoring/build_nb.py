#!/usr/bin/env python
"""Build (and optionally execute) a Jupyter notebook from a simple .nbsrc text file.

.nbsrc format: cells are separated by marker lines that contain ONLY
    <<<md>>>      -> the following lines form a markdown cell
    <<<code>>>    -> the following lines form a code cell
Everything is taken verbatim (no escaping needed for LaTeX, backslashes, quotes).
Leading/trailing blank lines of each cell are stripped.

Usage (always use the GestureDetection env python):
    PY=/home/dominik/miniconda3/envs/GestureDetection/bin/python
    $PY build_nb.py SRC.nbsrc OUT.ipynb            # build only
    $PY build_nb.py SRC.nbsrc OUT.ipynb --execute  # build + run all cells, store outputs
    $PY build_nb.py --check OUT.ipynb              # lint an existing notebook (links/images)

Execution runs with the notebook's own folder as working directory, a 600 s
per-cell timeout, and fails loudly (non-zero exit, traceback printed) if any
cell raises. After execution it also checks that every relative link and
image referenced in markdown exists on disk.
"""
import os
import re
import sys
import urllib.parse

import nbformat
from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook

MARK = re.compile(r"^<<<(md|code)>>>\s*$")


def parse(src_text):
    cells, kind, buf = [], None, []

    def flush():
        if kind is None:
            return
        text = "\n".join(buf).strip("\n")
        if text.strip():
            cells.append(new_markdown_cell(text) if kind == "md" else new_code_cell(text))

    for line in src_text.splitlines():
        m = MARK.match(line)
        if m:
            flush()
            kind, buf = m.group(1), []
        else:
            buf.append(line)
    flush()
    return cells


ROOT = "/home/dominik/Documents/other/mri"


def planned_notebooks():
    """Notebook paths (relative to ROOT) listed in the file map of AUTHORING_GUIDE.md."""
    guide = os.path.join(os.path.dirname(os.path.abspath(__file__)), "AUTHORING_GUIDE.md")
    planned, folder = set(), ""
    for line in open(guide, encoding="utf-8"):
        m = re.match(r"^([A-Za-z0-9_]+)/\s*$", line.strip())
        if m:
            folder = m.group(1)
            continue
        m = re.match(r"^\s*([A-Za-z0-9_.\-]+\.ipynb)", line)
        if m:
            name = m.group(1)
            planned.add(name if line[0] != " " else f"{folder}/{name}")
    return planned


def check_links(nb, nb_path):
    """Return list of problems with relative links / images in markdown cells.

    Links to notebooks that are planned in the guide's file map but not written yet are reported
    as 'pending' (not an error)."""
    base = os.path.dirname(os.path.abspath(nb_path))
    problems = []
    planned = planned_notebooks()
    pat = re.compile(r"!?\[[^\]]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)|<img[^>]+src=\"([^\"]+)\"")
    for i, c in enumerate(nb.cells):
        if c.cell_type != "markdown":
            continue
        for m in pat.finditer(c.source):
            target = m.group(1) or m.group(2)
            if re.match(r"^(https?:|mailto:|#)", target):
                continue
            path = urllib.parse.unquote(target.split("#")[0])
            if not path:
                continue
            full = os.path.normpath(os.path.join(base, path))
            rel = os.path.relpath(full, ROOT)
            if not os.path.exists(full) and rel in planned:
                print(f"  (pending: {rel} not written yet)")
            elif not os.path.exists(full):
                problems.append(f"cell {i}: missing target '{target}' -> {full}")
    return problems


WIDGET = re.compile(r"\binteract(ive)?\s*\(")


def execute(nb, nb_path):
    """Run all cells. Cells that create ipywidgets sliders (interact/interactive) are NOT executed:
    headless kernels sometimes hang on widget comms, and widget output is not stored anyway.
    They keep their source and a note, and run normally when the reader opens the notebook."""
    from nbclient import NotebookClient
    skipped = {}
    for i, c in enumerate(nb.cells):
        if c.cell_type == "code" and WIDGET.search(c.source):
            skipped[i] = c.source
            c.source = "pass"
    client = NotebookClient(nb, timeout=300, kernel_name="python3",
                            resources={"metadata": {"path": os.path.dirname(os.path.abspath(nb_path))}})
    try:
        client.execute()
    finally:
        for i, src in skipped.items():
            c = nb.cells[i]
            c.source = src
            c.outputs = [nbformat.v4.new_output("stream", name="stdout",
                         text="(Interactive widget cell: run this cell yourself in Jupyter to get the sliders.)\n")]
            c.execution_count = None
    return nb


def main(argv):
    if argv[0] == "--check":
        nb = nbformat.read(argv[1], as_version=4)
        probs = check_links(nb, argv[1])
        print("\n".join(probs) if probs else "links OK")
        return 1 if probs else 0
    src, out = argv[0], argv[1]
    do_exec = "--execute" in argv
    with open(src, encoding="utf-8") as f:
        cells = parse(f.read())
    nb = new_notebook(cells=cells)
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    nb.metadata["language_info"] = {"name": "python"}
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    rc = 0
    if do_exec:
        os.environ.setdefault("MPLBACKEND", "module://matplotlib_inline.backend_inline")
        try:
            execute(nb, out)
        except Exception as e:  # CellExecutionError carries the traceback text
            print("EXECUTION FAILED:\n", str(e)[-4000:])
            rc = 2
    nbformat.write(nb, out)
    probs = check_links(nb, out)
    if probs:
        print("LINK PROBLEMS:\n" + "\n".join(probs))
        rc = rc or 3
    n_md = sum(c.cell_type == "markdown" for c in cells)
    size_kb = os.path.getsize(out) / 1024
    print(f"wrote {out}: {len(cells)} cells ({n_md} md, {len(cells) - n_md} code), {size_kb:.0f} kB, "
          f"{'executed' if do_exec else 'not executed'}, rc={rc}")
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
