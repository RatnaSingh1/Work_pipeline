#this script collect the information from snp.txt from each folder and concatenate in snp_summries.txt

#!/bin/bash

set -e


echo "Saving summary..."
> ../snp_summaries.txt  # Clear or create snp_summaries.txt

for folder in *; do
    if [ -f "$folder/snps.txt" ]; then
        cat "$folder/snps.txt" >> snp_summaries.txt
        echo -e "\n" >> ../snp_summaries.txt  # Add a complete blank line after each file
    else
        echo "Warning: $folder/snps.txt not found"
    fi
done
