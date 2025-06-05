#!/bin/bash
 
salmon index -t Nta.fa.gz -i Nta_index

for fn in data/R*;
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i Nta_index -l A \
         -1 ${fn}/${samp}_1.fq.gz \
         -2 ${fn}/${samp}_2.fq.gz \
         -p 8 --validateMappings -o quants/${samp}_quant
done 
