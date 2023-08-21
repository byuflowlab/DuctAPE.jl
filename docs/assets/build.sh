#!/bin/bash

python3 gen_xdsm.py

pdflatex dtxdsm.tex
makeindex dtxdsm.nlo -s nomencl.ist -o dtxdsm.nls
pdflatex dtxdsm.tex
