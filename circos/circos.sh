
#!/bin/bash

# Iterate through the directories
for dir_path in *_assembly/; do
    # Ensure it's a directory
    if [[ -d "$dir_path" ]]; then
        for fasta_file in "$dir_path"/*.fasta; do
            # Ensure there is at least one matching fasta file
            if [[ -e "$fasta_file" ]]; then
                # Get the base name of the fasta file (without the directory path)
                base_name=$(basename "$fasta_file" .fasta)
                # Create a directory with this base name in the circos folder
                output_dir="circos/$base_name"
              
                # Run the mummer2circos command and save the output in the created directory
                apptainer exec mummer2circos.simg mummer2circos -r "$fasta_file" -q "$dir_path"/*assembly.fasta -f -a promer -o "$output_dir"
            fi
        done
    fi
done

