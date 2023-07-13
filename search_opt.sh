#!/bin/bash
### This code creates a list of Samples with hits of viral taxa when searched against a list of accession numbers from RefSeq


# Set the file containing the list of strings to search for
string_list_file="sequences.acc"

# Set the directory containing the files to parse
data_dir="/xdisk/uschuch/corykeith/VICAT_COTTON_FINAL/out_P002/all_blast"

# Set the output file for the matching strings
output_file="alphalist.txt"
PLATE_NUM="P002"
# Read the search strings into an array
readarray -t search_strings < "$string_list_file"

# Loop through all files in the data directory
for file in "$data_dir"/*; do
    # Skip directories
    if [[ -d "$file" ]]; then
        continue
    fi
    # Extract the filename and suffix from the file path using sed
    filename=$(echo "$file" | sed 's/.*all_blast\(.*\)_blastn.*/\1/')
    filename="${filename#/}"
    suffix="${filename}_filtered_scaffolds.fasta"
    echo "Processing $file: suffix=$suffix, filename=$filename"
    # Process each file
    while IFS=$'\t' read -r -a myArray; do
     Accession="${myArray[1]#*|}"
     Accession="${Accession%|}"
     # Search for each search string in the line
        for search_string in "${search_strings[@]}"; do
            # If the search string is found in the line, output the first column
            if awk -v search="$search_string" -v filename="$suffix" -v accession="$search_string" '$0 ~ search {print filename"\t"$1"\t"accession}' "$file" >> "$output_file"; then
                echo "$search_string found in $file"
            fi
        done
    done < "$file"
done

sort -u $output_file > ${P002}_sorted${output_file}
