#!/bin/bash

mkdir fastqc_output_all_2
mkdir fastq_trimmed_all_2



for f in $(ls Fastq/*_R1_001.fastq.gz | sed 's/Fastq\///' | sed 's/_R1_001.fastq.gz//'); do
  echo "trim_galore -j 8 --fastqc  --fastqc_args \"-o fastqc_output_all_2/\" --paired --max_n 5 -o fastq_trimmed_all_2/ Fastq/${f}_R1_001.fastq.gz Fastq/${f}_R2_001.fastq.gz"
  
  trim_galore -j 8 -a CTGTCTCTTATACACATCT  --fastqc  --quality 25   --fastqc_args "-o fastqc_output_all_2/" --paired --max_n 5 -o fastq_trimmed_all_2/ Fastq/${f}_R1_001.fastq.gz Fastq/${f}_R2_001.fastq.gz 
done


