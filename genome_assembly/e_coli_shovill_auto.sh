#!/bin/bash

mkdir shovill_assembly

# Loop through each _R1 file in the strepto directory
for fn1 in *_R1_001_val_1.fq; do
    
    # Extract the common part of the file name
    samp=$(basename "${fn1}" "_R1_001_val_1.fq")
    
    # Check if the corresponding _R2 file exists
    fn2="${samp}_R2_001_val_2.fq"
    if [ -f "$fn2" ]; then
        echo "Processing sample ${samp}"
        
        # Run megahit with both _R1 and _R2 input files
       singularity exec /home/ratna/shovill.sif shovill --R1 "${fn1}" --R2 "${fn2}" --outdir "shovill_assembly/${samp}_assembly" --cpus 14 

    else
        echo "Error: Corresponding _R2 file not found for ${fn1}"
    fi
done
