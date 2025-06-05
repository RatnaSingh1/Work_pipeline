#!/bin/bash



# Loop through each _R1 file in the strepto directory
for fn1 in *_R1_001_val_1.fq; do
    
    # Extract the common part of the file name
    samp=$(basename "${fn1}" "_R1_001_val_1.fq")
    
    # Check if the corresponding _R2 file exists
    fn2="${samp}_R2_001_val_2.fq"
    if [ -f "$fn2" ]; then
        echo "Processing sample ${samp}"
        
        # Run megahit with both _R1 and _R2 input files
       spades.py  --meta -1 "${fn1}" -2 "${fn2}" -o "spades_assembly/${samp}_assembly" -t 14 
    else
        echo "Error: Corresponding _R2 file not found for ${fn1}"
    fi
done


