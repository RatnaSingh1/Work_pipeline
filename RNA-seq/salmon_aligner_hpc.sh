#!/bin/bash

#SBATCH --nodes=1                 
#SBATCH --ntasks-per-node=72                         
#SBATCH --partition=normal       
#SBATCH --time=24:00:00            
#SBATCH --job-name=star         
#SBATCH --output=output.dat        
#SBATCH --mail-type=ALL            
#SBATCH --mail-user=

module load palma/2020b  GCC/10.2.0  OpenMPI/4.0.5
module load Salmon/1.4.0


 
salmon index -t transcriptome.fa -i index

for fn in data/R*;
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i index -l A \
         -1 ${fn}/${samp}_1.fq.gz \
         -2 ${fn}/${samp}_2.fq.gz \
         -p 8 --validateMappings -o quants/${samp}_quant
done 
