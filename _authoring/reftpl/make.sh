#!/bin/bash
# usage: make.sh Name   (expands reftpl/Name.in -> nbsrc/reference/Name.nbsrc, builds+executes into ROOT/reference)
S=/tmp/claude-1000/-home-dominik-Documents-other-mri/cfd9f6c6-3c4e-45cd-a910-a339740d0872/scratchpad
PY=/home/dominik/miniconda3/envs/GestureDetection/bin/python
$PY $S/reftpl/expand.py $S/reftpl/$1.in $S/nbsrc/reference/$1.nbsrc || exit 1
$PY $S/build_nb.py $S/nbsrc/reference/$1.nbsrc /home/dominik/Documents/other/mri/reference/$1.ipynb --execute 2>&1 | grep -v pending
echo "pending links: $($PY $S/build_nb.py --check /home/dominik/Documents/other/mri/reference/$1.ipynb 2>&1 | grep -c pending)"
