

#!/bin/bash




# Define the input directory containing the FNA files
input_directory="shovill_assembly"

# Define the output directory
output_parent_directory="output_data_shovill_e_coli_prokka_scaffolds"

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
for fna_file in "$input_directory"/*/spades.fasta; do
    if [ -f "$fna_file" ]; then
        echo "Processing file: $fna_file"
        # Extract the foldername without extension
       foldername=$(basename -- "$(dirname -- "$fna_file")")

        # Create a directory for the output
        output_directory="$output_parent_directory/$foldername"
        mkdir -p "$output_directory"

        # Run Prokka on the FNA file and output to the corresponding directory
        echo "Running Prokka on $fna_file"
        prokka --kingdom Viruses --species Bacteriophage --centre X   --cpus 15 --outdir "$output_directory" --prefix "$foldername" --force "$fna_file"   
prokka_exit_code=$?

        if [ "$prokka_exit_code" -ne 0 ]; then
            echo "Error: Prokka failed for file '$fna_file'."
            exit 1
        fi
        # Find the GFF file in the Prokka output directory
        gff_file=$(find "$output_directory" -name "*.gff" -type f)

        # Move the GFF file to the GFF folder
        echo "Copying GFF file to $gff_folder/${foldername}.gff"
        cp "$gff_file" "$gff_folder/${foldername}.gff"
        echo "Finished processing $fna_file"
    fi
done

echo "Processing completed successfully."


# Create the annotation directory outside the loop
#mkdir prokka_annotation

# Loop through each _L001_assembly.tsv file in the output_data_shovill_e_coli directory
#for fn1 in output_data_shovill_e_coli/*/*_L001_assembly.tsv; do
 #   
  #  # Extract the common part of the file name
   # samp=$(basename "${fn1}" "_L001_assembly.tsv")
    #
    # Apply awk to filter the lines and save the result to the annotation directory
    #awk '$3 != "hypothetical" && $3 != "hypothetical protein"' "$fn1" > "annotation/${samp}_assembly.txt"
    
#done
