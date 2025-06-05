#!/bin/bash

# Define the input directory containing the FNA files
input_directory="split_sequences"

# Define the output directory
output_parent_directory="output_data_pangenome"

# Define the directory to store GFF files
gff_folder="gff_folder"

# Create the GFF folder if it doesn't exist
mkdir -p "$gff_folder"

# Check if prokka is available
if ! command -v prokka &> /dev/null; then
    echo "Error: prokka not found. Make sure it's installed and accessible in your PATH."
    exit 1
fi

echo "Starting the processing..."

# Loop through each FNA file in the input directory
for fna_file in "$input_directory"/*.fasta; do
    if [ -f "$fna_file" ]; then
        echo "Processing file: $fna_file"
        # Extract the filename without extension
        filename=$(basename -- "$fna_file" .fna)

        # Create a directory for the output
        output_directory="$output_parent_directory/$filename"
        mkdir -p "$output_directory"

        # Run Prokka on the FNA file and output to the corresponding directory
        echo "Running Prokka on $fna_file"
        prokka --kingdom Viruses --species Bacteriophage --centre X --compliant --addmrna   --rfam  --cpus 15 --outdir "$output_directory" --prefix "$filename" --force "$fna_file"   
prokka_exit_code=$?

        if [ $prokka_exit_code -ne 0 ]; then
            echo "Error: Prokka failed for file '$fna_file'."
            exit 1
        fi

        # Find the GFF file in the Prokka output directory
        gff_file=$(find "$output_directory" -name "*.gff" -type f)

        # Move the GFF file to the GFF folder
        echo "Moving GFF file to $gff_folder/${filename}.gff"
       cp "$gff_file" "$gff_folder/${filename}.gff"
        echo "Finished processing $fna_file"
    fi
done

echo "Processing completed successfully."

cd *_genomic/gff_folder

roary --mafft -e -n -qc  -r -y -p 15 *.gff -o clustered_proteins -f output_rorary

