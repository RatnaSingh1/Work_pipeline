

#!/bin/bash

set -e

# Create the fastqc_report folder if it doesn't exist
mkdir -p fastqc_report

echo "Saving summary..."
> fastqc_report/fastqc_data_summaries.txt  # Clear or create fastqc_data_summaries.txt

for folder in fastqc_output_all_2/*_fastqc; do
    if [ -f "$folder/fastqc_data.txt" ]; then
        # Extract lines from >>Basic Statistics to >>END_MODULE and append to summary
        awk '/#Measure/,/%GC/' "$folder/fastqc_data.txt" >> fastqc_report/fastqc_data_summaries.txt
        echo -e "\n" >> fastqc_report/fastqc_data_summaries.txt  # Add a blank line after each file
    else
        echo "Warning: $folder/fastqc_data.txt not found"
    fi
done