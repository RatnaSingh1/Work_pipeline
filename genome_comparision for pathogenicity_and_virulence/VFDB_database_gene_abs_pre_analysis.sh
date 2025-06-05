#!/bin/bash

#from pan genome analysis gene list gene_presence_absence.csv obtained

#*.fna run again the PHI database (phi-base_current.fas) to find the putative pathogenic genes with a cutoff 30%identity and 50% query cover output_parameters_50.txt

#steps to prepare header

 echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen" > output_parameters_70.txt

#for the output result list

blastp -query PROKKA_10022024.faa -subject VFDB_setA_pro.fas  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' | awk '$3 >= 70 && ($4/$14)*100 >= 50' >> output_parameters_70.txt


# Extract the header from gene_presence_absence.csv and store it in a temporary file
head -n 1 gene_presence_absence.csv >> temp_matched_columns.txt

# Append the common rows from gene_presence_absence.csv to the temporary file
grep -Ff <(cut -f1 output_parameters_70.txt) gene_presence_absence.csv >> temp_matched_columns.txt

# Rename the temporary file to matched_columns1.txt
mv temp_matched_columns.txt matched_columns.txt




#in this way gene present in orthologous and there function reterived

#In next step fasta file reterived for the related gene



#grepping first gene list from output_parameters_70.txt

 awk '{print $1}' output_parameters_70.txt > matched_gene_list.txt

seqtk subseq PROKKA_10022024.faa matched_gene_list.txt > filter_sequence.faa


#to remove sequences from filter_sequence.faa that have a length of less than 100 amino acids, 

awk '/^>/ {if (length(seq) >= 100) {print header ORS seq}; header = $0; seq = ""; next} {seq = seq $0} END {if (length(seq) >= 100) {print header ORS seq}}' filter_sequence.faa > filtered_sequences.faa

# copy the header of filtered_sequences.faa which is in this format 
#than  take only the text just after the sign > , like this NKFGOEFH_00004 and exclude everything befor and after it and #then comapre the #extracted with a another file named output.txt and match and copy the txt from there to a new file


# Extract gene names from filtered_sequences.faa
awk '/^>/ {gsub(/[()]/, "", $0); split($0, arr, " "); print substr(arr[1], 2)}' filtered_sequences.faa > gene_names.txt

# Compare with gene_names.txt and copy corresponding text
awk 'NR==FNR{genes[$0]; next} $1 in genes {print $0}' gene_names.txt output_parameters_70.txt > matched_output_final.txt




#to include the header in the matched_columns_final.txt

# Extract the header from gene_presence_absence.csv and store it in a temporary file
head -n 1 gene_presence_absence.csv > temp_matched_columns1.txt

# Append the common rows from gene_presence_absence.csv to the temporary file
grep -Ff <(cut -f1 matched_output_final.txt) gene_presence_absence.csv >> temp_matched_columns1.txt

# Rename the temporary file to matched_columns1.txt
mv temp_matched_columns1.txt matched_columns_final.txt




