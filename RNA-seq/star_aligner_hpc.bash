#!/bin/bash

 
#SBATCH --nodes=1                 
#SBATCH --ntasks-per-node=72                         
#SBATCH --partition=normal       
#SBATCH --time=24:00:00            
 
#SBATCH --job-name=star         
#SBATCH --output=output.dat        
#SBATCH --mail-type=ALL            
#SBATCH --mail-user=


module load  palma/2019a  GCC/8.2.0-2.31.1

module load STAR/2.7.3a
STAR-Fusion/1.8.0-Perl-5.28.1-Python-3.7.2
 
 STAR --runMode genomeGenerate --genomeDir index --genomeFastaFiles Nbv5_transcriptome.fa --sjdbGTFfile Nbv5tr_Nbv05.gtf --sjdbOverhang 50 --outFileNamePrefix index --limitGenomeGenerateRAM 99158965632


for fn in data/R*;
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"

         
         STAR --genomeDir index --readFilesIn ${fn}/${samp}_1.fq.gz ${fn}/${samp}_2.fq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts  --outFileNamePrefix alignments//${samp}_quant --runThreadN 14
done 


featureCounts -T 8 -s 2 -p -t exon -g gene_id -a genes.gtf -o counts.txt alignments/*.bam
