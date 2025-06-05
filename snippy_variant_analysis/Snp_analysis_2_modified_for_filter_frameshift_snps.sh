#!/bin/bash
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
        output_file="$folder/${folder_name}_snp_comparison_summary_filter.txt"
        > "$output_file"  # Clear or create the output file

        # Initialize the output file with headers
        echo -e "Product\tCount\tFolders\tEffect" > "$output_file"

        # Temporary file to hold entries with folder info
        temp_file="$folder/temp_snp_data.txt"
        > "$temp_file"

        # Process each subfolder within the current main folder
        for subfolder in "$folder"/*; do
            if [ -d "$subfolder" ]; then
                # Extract the subfolder name
                subfolder_name=$(basename "$subfolder")

                # Define the SNP file path in the subfolder
                snp_file="$subfolder/filter_leftover_snps.tab"

                if [ -f "$snp_file" ]; then
                    echo "Processing SNP file: $snp_file"

                    # Extract Product (Column14) and Effect (Column11), skipping empty fields
                    awk -v folder="$subfolder_name" \
                        'BEGIN {FS="\t"; OFS="\t"} NR > 1 {if ($14 != "" && $11 != "") print $14, folder, $11}' \
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
            # Extract product, folder, and annotation
            product = $1
            folder = $2
            annotation = $3

            # Increment the count for each unique product
            product_count[product]++

            # Append the folder name to the list for each product if not already added
            if (!seen[product FS folder]++) {
                folders[product] = (folders[product] ? folders[product] ", " : "") folder
            }

            # Store the annotation for each product (assumes annotations are the same for a product)
            annotations[product] = annotation
        } END {
            # Output each product with its count, folders, and annotation
            for (prod in product_count) {
                print prod "\t" product_count[prod] "\t" folders[prod] "\t" annotations[prod]
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
