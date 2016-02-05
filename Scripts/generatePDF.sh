#!/bin/bash

echo "" > output.tex

python generating_Jesse.py -n "${1}" -num "${2}" -t "${3}"


pdflatex output.tex

rm output.log
rm output.aux
rm output.tex
