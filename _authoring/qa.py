"""QA over all notebooks in the course: structure, errors, links, size."""
import glob
import os
import re
import sys

import nbformat

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_nb  # noqa: E402

ROOT = build_nb.ROOT
planned = build_nb.planned_notebooks()
existing = {os.path.relpath(p, ROOT) for p in glob.glob(f"{ROOT}/**/*.ipynb", recursive=True)
            if ".ipynb_checkpoints" not in p and "/_build/" not in p}

print(f"planned {len(planned)}, existing {len(existing & planned)}")
missing = sorted(planned - existing)
if missing:
    print("MISSING:", *missing, sep="\n  ")
extra = sorted(existing - planned)
if extra:
    print("EXTRA (not in map):", *extra, sep="\n  ")

REQ = {
    "symbols table": re.compile(r"^##+ .*Symbols", re.M),
    "takeaways": re.compile(r"^##+ .*Key takeaways", re.M | re.I),
    "questions": re.compile(r"^##+ .*Check your understanding", re.M | re.I),
    "attribution": re.compile(r"Larson|New chapter written|adaptation", re.I),
}
rows = []
for rel in sorted(existing & planned):
    path = os.path.join(ROOT, rel)
    nb = nbformat.read(path, as_version=4)
    md = "\n".join(c.source for c in nb.cells if c.cell_type == "markdown")
    codes = [c for c in nb.cells if c.cell_type == "code"]
    errs = sum(1 for c in codes for o in c.get("outputs", []) if o.get("output_type") == "error")
    unexec = sum(1 for c in codes if c.get("execution_count") is None and c.source.strip() and not build_nb.WIDGET.search(c.source))
    imgs = sum(1 for c in codes for o in c.get("outputs", []) if "image/png" in o.get("data", {}))
    links = [p for p in build_nb.check_links(nb, path)]
    # links to planned-but-missing notebooks are printed by check_links as pending; recount here
    dead = []
    for c in nb.cells:
        if c.cell_type != "markdown":
            continue
        for m in re.finditer(r"\]\(([^)\s]+\.ipynb)[^)]*\)", c.source):
            t = os.path.normpath(os.path.join(os.path.dirname(path), m.group(1).split("#")[0]))
            if not os.path.exists(t):
                dead.append(os.path.relpath(t, ROOT))
    missing_req = [k for k, rx in REQ.items() if not rx.search(md)]
    if rel.startswith("reference/") or rel == "00_Start_Here.ipynb":
        missing_req = [k for k in missing_req if k == "attribution"]
    first = codes[0].source if codes else ""
    hdr_ok = (not codes) or ("mri_tools" in first and "sys.path" in first)
    size = os.path.getsize(path) / 1e6
    rows.append((rel, len(nb.cells), len(md) // 1000, len(codes), imgs, errs, unexec, size,
                 missing_req, links, sorted(set(dead)), hdr_ok))

print(f"\n{'notebook':55s} cells mdk code figs err unex   MB  issues")
for r in rows:
    rel, n, mdk, nc, imgs, errs, unexec, size, mreq, links, dead, hdr = r
    issues = []
    if errs: issues.append(f"{errs} ERROR outputs")
    if unexec: issues.append(f"{unexec} unexecuted")
    if mreq: issues.append("missing " + "/".join(mreq))
    if links: issues.append(f"{len(links)} bad links/images")
    if dead: issues.append(f"dead nb links: {dead}")
    if not hdr: issues.append("non-standard first code cell")
    if size > 5: issues.append("LARGE")
    print(f"{rel:55s} {n:5d} {mdk:3d} {nc:4d} {imgs:4d} {errs:3d} {unexec:4d} {size:5.1f}  {'; '.join(issues)}")
    for l in links:
        print("      ", l)
tot = sum(r[7] for r in rows)
print(f"\ntotal size {tot:.1f} MB, total markdown {sum(r[2] for r in rows)}k chars")
