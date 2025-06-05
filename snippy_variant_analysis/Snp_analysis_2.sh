#!/bin/bash

# this script use to compare the snps.tab from each folder to get list of genes that are unique and common in samples – output snp_comparison_summary
#First Count (Count in the otput file): represents the total number of mutations (or occurrences of a particular Product) in all the samples (folders). Each instance of a mutation #across samples is counted, irrespective of whether it appears in the same or different folders.
#Second Count (Sample_Count in the output file): This represents the number of unique samples (folders) in which a specific mutation (Product) is present. It shows how widespread a #mutation is across the samples.


# Directory containing folders with SNP files
base_dir="$(pwd)"

# List of main folders to process
folders=("e_coli" "plasmid1" "plasmid2" "plasmid3" "plasmid4" "plasmid5")

# Loop through each main folder
for folder_name in "${folders[@]}"; do
    folder="$base_dir/$folder_name"

    if [ -d "$folder" ]; then
        echo "Processing main folder: $folder_name"

        # Define the output file with the base folder name included
        output_file="$folder/${folder_name}_snp_comparison_summary.txt"
        > "$output_file"  # Clear or create the output file

        # Initialize the output file with headers
        echo -e "Product\tCount\tFolders\tEffect\tGene\tSample_Count" > "$output_file"

        # Temporary file to hold Product entries with folder info
        temp_file="$folder/temp_snp_data.txt"
        > "$temp_file"

        # Process each subfolder within the current main folder
        for subfolder in "$folder"/*; do
            if [ -d "$subfolder" ]; then
                # Extract the subfolder name
                subfolder_name=$(basename "$subfolder")

                # Define the SNP file path in the subfolder
                snp_file="$subfolder/snps.tab"

                if [ -f "$snp_file" ]; then
                    echo "Processing SNP file: $snp_file"

                    # Extract Product (Column14), Effect (Column11), and Gene (Column13), skipping empty fields
                    awk -v folder="$subfolder_name" \
                        'BEGIN {FS="\t"; OFS="\t"} NR > 1 {if ($14 != "" && $11 != "" && $13 != "") print $14, folder, $11, $13}' \
                        "$snp_file" >> "$temp_file"
                else
                    echo "Warning: No $snp_file file found in $subfolder"
                fi
            fi
        done

        # Step 2: Sort the data to ensure identical Product entries are grouped together
        sorted_temp_file="$folder/sorted_temp_snp_data.txt"
        sort -k1,1 "$temp_file" > "$sorted_temp_file"

        # Step 3: Process the sorted data to count occurrences and gather folder names
        awk -F "\t" '{
            # Extract product, folder, effect, and gene
            product = $1
            folder = $2
            effect = $3
            gene = $4

            # Increment the count for each unique product
            product_count[product]++

            # Append the folder name to the list for each product if not already added
            if (!seen[product FS folder]++) {
                folders[product] = (folders[product] ? folders[product] ", " : "") folder
                Sample_Count[product]++  # Increment the folder count for the product
            }

            # Store the effect for each product (assumes effects are the same for a product)
            effects[product] = effect

            # Store the gene for each product (assumes genes are the same for a product)
            genes[product] = gene
        } END {
            # Output each product with its count, folders, effect, gene, and folder count
            for (prod in product_count) {
                print prod "\t" product_count[prod] "\t" folders[prod] "\t" effects[prod] "\t" genes[prod] "\t" Sample_Count[prod]
            }
        }' "$sorted_temp_file" >> "$output_file"

        echo "Results saved in $output_file"

        # Optional: Clean up temporary files
        rm "$temp_file" "$sorted_temp_file"
    else
        echo "Warning: Main folder $folder does not exist"
    fi

done

echo "Processing complete."
