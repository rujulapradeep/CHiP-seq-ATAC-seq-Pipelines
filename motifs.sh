#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <input_file> <output_dir>"
  exit 1
fi

# Assign the positional parameters to variables
input_file="$1"
output_dir="$2"

# Check if the input file exists
if [ ! -f "$input_file" ]; then
  echo "Input file '$input_file' does not exist."
  exit 1
fi

# Check if the output directory exists or create it
if [ ! -d "$output_dir" ]; then
  mkdir -p "$output_dir"
fi

# Function to create directories and run Homer analysis for each group
find_motifs() {
  
  # Create directory for the group within the output directory
  group_dir="${output_dir}/${group}_motifs"
  mkdir -p "$group_dir"
  
  # Run Homer analysis for female peaks
  echo "Running Homer analysis for female peaks in group: $group"
  findMotifsGenome.pl "$female_peaks" hg38 "$group_dir/homer_results_female/" -size given
  
  # Run Homer analysis for male peaks
  echo "Running Homer analysis for male peaks in group: $group"
  findMotifsGenome.pl "$male_peaks" hg38 "$group_dir/homer_results_male/" -size given
}

# Read the input file line by line and process each line
while read -r line; do
  # Split the line into fields using space as delimiter
  fields=($line)
  
  # Extract individual fields
  female_peaks="${fields[0]}"
  male_peaks="${fields[1]}"
  group="${fields[2]}"
  
  echo "Finding Motifs: $group"
  
  # Call the function to create directories and run Homer analysis
  find_motifs "$group" "$female_peaks" "$male_peaks"
done < "$input_file"

echo "Script execution complete."
