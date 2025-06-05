#!/bin/bash

# set the number of nodes
#SBATCH --nodes=1

# set the number of CPU cores per node
#SBATCH --ntasks-per-node 30

# set a partition
#SBATCH --partition normal

# set max wallclock time
#SBATCH --time=20:00:00

# set name of job
#SBATCH --job-name=test123

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# set an output file
#SBATCH --output output.dat

# send mail to this address
#SBATCH --mail-user=singhr@uni-muenster.de


 module load  palma/2019a
module load QIIME2/2019.7



# Set the path to the directory containing the paired-end FASTQ files
#SEQ_DIR=$fastq

# Import the paired-end FASTQ files into QIIME 2
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path rename_trim_fastqc \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path paired-end-demux.qza

# Perform quality control on the paired-end FASTQ files
qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization paired-end-demux-summary.qzv

# Denoise the paired-end sequences to generate feature tables and representative sequences
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --p-trim-left-f 13 \
  --p-trim-left-r 13 \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 225 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats denoising-stats.qza

# Tabulate the representative sequences
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv  

#de novo chimera checking
qiime vsearch uchime-denovo --i-sequences rep-seqs.qza --i-table table.qza --o-chimeras chimera-denova.qza --o-nonchimeras non-chimera-denova.qza --o-stats chimera-denova-stats.qza

#filterout chimeras
qiime feature-table filter-features --m-metadata-file non-chimera-denova.qza --i-table table.qza --o-filtered-table table-nonchimeric.qza
qiime feature-table filter-seqs --i-data rep-seqs.qza --m-metadata-file non-chimera-denova.qza --o-filtered-data rep-seqs-nonchimeric.qza

#clustering
#E2. Clustering of FeatureTable[Frequency] and FeatureData[Sequence] (final denovo filtering)
qiime vsearch cluster-features-de-novo --i-table table-nonchimeric.qza --i-sequences rep-seqs-nonchimeric.qza --o-clustered-table table-dn-97.qza --o-clustered-sequences rep-seqs-dn-97.qza --p-perc-identity 0.97 


# Assign taxonomy to the representative sequences using the SILVA database
qiime feature-classifier classify-sklearn \
  --i-classifier ./Classifier/silva-classifier-138-99.qza \
  --i-reads rep-seqs-dn-97.qza \
  --o-classification taxonomy-dn-97_138.qza

qiime metadata tabulate --m-input-file taxonomy-dn-97_138.qza --o-visualization taxonomy-dn-97_138.qzv

#G1.5 Taxonomy-based filtering of tables 

qiime taxa filter-table --i-table table-dn-97.qza --i-taxonomy taxonomy-dn-97_138.qza --p-exclude D_4__Mitochondria,D_3__Chloroplast --o-filtered-table table-final.qza 

#G1.5 Taxonomy-based filtering of sequences
qiime taxa filter-seqs --i-sequences rep-seqs-dn-97.qza --i-taxonomy taxonomy-dn-97_138.qza --p-exclude D_4__Mitochondria,D_3__Chloroplast --o-filtered-sequences rep-seqs-final_138.qza 

# Tabulate the representative sequences
qiime feature-table tabulate-seqs \
  --i-data rep-seqs-final_138.qza  \
  --o-visualization rep-seqs-final_138.qzv

## Tabulate the table
qiime feature-table summarize --i-table table-final.qza --o-visualization table-final.qzv
# Generate a taxonomic bar plot
qiime taxa barplot --i-table table-final.qza  --i-taxonomy taxonomy-dn-97_138.qza --m-metadata-file metadata.tsv --o-visualization taxa-bar-plots_138.qzv

