# this script take the snp_summries.txt and generate the variant_summary.csv, which is usefil for generating plot using R

#!/bin/bash

# Define output CSV file
output_csv="variant_summary.csv"

# Add header to the CSV file
echo "Sample,Variant-DEL,Variant-INS,Variant-SNP,Variant-COMPLEX,Variant-MNP,VariantTotal" > "$output_csv"

# Process each entry in the snp_summaries.txt file
# Read the snp_summaries.txt file line by line
while read -r line; do
    # Extract relevant information from each line
    case "$line" in
        DateTime*)
            # Reset variables for a new sample
            sample_name=""
            variant_del=0
            variant_ins=0
            variant_snp=0
            variant_complex=0
            variant_mnp=0
            variant_total=0
            ;;
        ReadFiles*)
            # Extract base name from the ReadFiles path
            # Assuming the base name is the first part before "_R1" or "_R2"
            read_file_path=$(echo "$line" | cut -f2 -d$'\t')
            sample_name=$(basename "$read_file_path" | sed 's/_R[12].*//')
            ;;
        Variant-DEL*)
            variant_del=$(echo "$line" | cut -f2 -d$'\t')
            ;;
        Variant-INS*)
            variant_ins=$(echo "$line" | cut -f2 -d$'\t')
            ;;
        Variant-SNP*)
            variant_snp=$(echo "$line" | cut -f2 -d$'\t')
            ;;
        Variant-COMPLEX*)
            variant_complex=$(echo "$line" | cut -f2 -d$'\t')
            ;;
        Variant-MNP*)
            variant_mnp=$(echo "$line" | cut -f2 -d$'\t')
            ;;
        VariantTotal*)
            variant_total=$(echo "$line" | cut -f2 -d$'\t')
            # Write the current sample's data to the CSV file
            if [ -n "$sample_name" ]; then
                echo "$sample_name,$variant_del,$variant_ins,$variant_snp,$variant_complex,$variant_mnp,$variant_total" >> "$output_csv"
            fi
            ;;
    esac
done < snp_summaries.txt

echo "CSV file generated: $output_csv"
