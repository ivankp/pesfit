#!/bin/bash

for s in 68 90
do
  ./bin/confnoteplot2 -c plot2.cfg -o plot2_${s}_.pdf --sigma-frac=0.$s
  pdftk A=plot2_${s}_.pdf cat A1 output plot2_${s}.pdf
  ./scripts/fix_pdf plot2_${s}.pdf
done
