#!/bin/bash
 
# set the number of nodes
#SBATCH --nodes=1

# set the number of CPU cores per node
#SBATCH --ntasks-per-node 4

# How much memory is needed (per node). Possible units: K, G, M, T
#SBATCH --mem=100G

# set a partition
#SBATCH --partition gpuv100

# set numer of gpus
#SBATCH --gres=gpu:1

# set max wallclock time
#SBATCH --time=1-00:00:00

# set name of job
#SBATCH --job-name=Alphafold

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL
#SBATCH --output output1.dat
# run the application



ml palma/2020b fosscuda/2020b
ml AlphaFold/2.1.1


export ALPHAFOLD_DATA_DIR=/Applic.HPC/data/alphafold

cd /scratch/tmp/singhr/alphafold 

alphafold \
--fasta_paths MD.fasta \
--model_preset monomer \
--output_dir results1/results \
--max_template_date 2021-07-31 \
--data_dir /Applic.HPC/data/alphafold/ \

#ENV

