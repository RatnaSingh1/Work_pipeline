#!/bin/bash
for f in *.pdbqt; do  b=`basename $f .pdbqt`; echo Processing ligand $b; mkdir -p data02; vina --config conf.txt --ligand $f --out data02/$f.pdbqt --log data02/$f.txt; done
cd data02
grep " 1 " *.txt | cut -c1-13,20-45 >>result2
grep " 2 " *.txt | cut -c1-13,20-45 >>result2
grep " 3 " *.txt | cut -c1-13,20-45 >>result2
grep " 4 " *.txt | cut -c1-13,20-45 >>result2
grep " 5 " *.txt | cut -c1-13,20-45 >>result2
