#!/bin/sh 
#PBS -o output.dat
#PBS -l walltime=48:00:00,nodes=1:westmere:ppn=12
#PBS -A r0bpmo

#PBS -m ae
#PBS -q default
#PBS -N DP6
#PBS -j oe
cd $PBS_O_WORKDIR
module load palma/2019a

module load FastQC/0.11.9-Java-11
module add bowtie2
module add palma/2019b  GCC/8.3.0  OpenMPI/3.1.4
modue add HISAT2/2.2.1

module add palma/2019a  GCC/8.2.0-2.31.1
module add SAMtools/1.9
module load Subread/2.0.0




# Step 1: Extract splice sites and exon locations from the GTF file
extract_splice_sites.py chrX_data/genes/genes.gtf > chrX.ss
extract_exons.py chrX_data/genes/genes.gtf > chrX.exon

# Step 2: Build HISAT2 index with exon and splice site information
hisat2-build -p 14 --ss chrX.ss --exon chrX.exon genome.fa chrX_tran

# Step 3: Create directory for mapped data
mkdir map

# Step 4: Map RNA-seq reads to the indexed genome using HISAT2
# Mapping samples from different conditions and time points
for sample in Sam-CS-2h-1 Sam-CS-2h-2 Sam-CS-2h-3 Sam-CS-5h-1 Sam-CS-5h-2 Sam-CS-5h-3 \
              Sam-H2O-2h-1 Sam-H2O-2h-2 Sam-H2O-2h-3 Sam-H2O-5h-1 Sam-H2O-5h-2 Sam-H2O-5h-3
do
  hisat2 -p 14 -x chrX_tran -q samples/${sample}.fastq -S map/${sample}M.sam
done

# Step 5: Convert SAM files to sorted BAM files using samtools
for sample in Sam-CS-2h-1 Sam-CS-2h-2 Sam-CS-2h-3 Sam-CS-5h-1 Sam-CS-5h-2 Sam-CS-5h-3 \
              Sam-H2O-2h-1 Sam-H2O-2h-2 Sam-H2O-2h-3 Sam-H2O-5h-1 Sam-H2O-5h-2 Sam-H2O-5h-3
do
  samtools sort -@ 8 -o map/${sample}B.bam map/${sample}M.sam
done

#step 6 : Assemble Transcripts Using StringTie

mkdir assembly
for sample in Sam-CS-2h-1 Sam-CS-2h-2 Sam-CS-2h-3 Sam-CS-5h-1 Sam-CS-5h-2 Sam-CS-5h-3 \
              Sam-H2O-2h-1 Sam-H2O-2h-2 Sam-H2O-2h-3 Sam-H2O-5h-1 Sam-H2O-5h-2 Sam-H2O-5h-3
do
  stringtie map/${sample}.bam -l ${sample}S -p 14 -G genes.gtf -o assembly/${sample}.gtf
done

#step 7 Merge Transcript Assemblies

stringtie --merge -p 8 -G genes.gtf -o stringtie_merged.gtf mergelist.txt

# Step 6: Count reads  genes using featureCounts
featureCounts -T 8 -s 2 -p -t exon -g gene_id -a stringtie_merged.gtf -o counts.txt map/*.bam

