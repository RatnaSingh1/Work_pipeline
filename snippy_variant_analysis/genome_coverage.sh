#!/bin/bash
# geneome coverage from Bam output file
# Directory containing folders with SNP files (current working directory)
base_dir="$(pwd)"

# List of folders to process (specify main folder names here)
folders=("plasmid1" "plasmid2" ..........)

# Create a file to store the summary of coverage results
summary_file="$base_dir/genome_coverage_summary_1.txt"
echo -e "Subfolder\tRName\tStartPos\tEndPos\tNumReads\tCovBases\tCoverage\tMeanDepth\tMeanBaseQ\tMeanMapQ" > "$summary_file"

# Loop through each folder in the base directory
for folder in "${folders[@]}"; do
    # Define the path to the main folder
    main_folder="$base_dir/$folder"
    
    # Check if the main folder exists
    if [[ -d "$main_folder" ]]; then
        # Loop through each subfolder in the current main folder
        for subfolder in "$main_folder"/*/; do
            # Check if the subfolder contains a snps.bam file
            bam_file="$subfolder/snps.bam"
            
            if [[ -f "$bam_file" ]]; then
                # Run samtools coverage and capture the output
                coverage_output=$(samtools coverage "$bam_file")
                
                # Extract the relative path (e.g., plasmid1/folder1) from the full path
                relative_path="${subfolder#$base_dir/}"
                
                # Append the relative path and coverage output to the summary file
                echo "$coverage_output" | awk -v folder="$relative_path" 'NR > 1 {print folder "\t" $0}' >> "$summary_file"
            else
                echo "No snps.bam found in $subfolder" >&2
            fi
        done
    else
        echo "Directory $main_folder does not exist!" >&2
    fi
done

# Final summary file has been saved in the base directory
echo "Genome coverage summary has been saved to $summary_file"
