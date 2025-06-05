#!/bin/bash
 


for fn in data/R*;
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"

         
         STAR --genomeDir indexes --readFilesIn ${fn}/${samp}_1.fq.gz ${fn}/${samp}_2.fq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts  --outFileNamePrefix alignments1//${samp}_quant --runThreadN 14
done 


featureCounts -T 8 -s 2 -p -t exon -g gene_id -a genes.gtf -o counts.txt alignments1/*.bam
