#!/bin/bash

# Function to check if a command exists
command_exists() {
  command -v "$1" >/dev/null 2>&1
}

# Function to install a package using pip
install_with_pip() {
  pip install "$1"
}

# Function to install a package using conda
install_with_conda() {
  conda install -y "$1"
}

# Define input file and output directory
input_file=$1
output_dir=$2

# Check if bedtools is installed
if ! command_exists bedtools; then
  echo "bedtools is not installed. Attempting to install with pip..."
  if command_exists pip; then
    install_with_pip bedtools
  else
    echo "pip is not installed. Attempting to install with conda..."
    if command_exists conda; then
      install_with_conda bedtools
    else
      echo "Neither pip nor conda is installed. Please install bedtools manually."
      exit 1
    fi
  fi
fi

# Read input file line by line
while IFS=$'\t' read -r group_name && IFS=$'\t' read -r -a female_files && IFS=$'\t' read -r -a male_files; do
    # Create directory for the group
    group_dir="${output_dir}/${group_name}"
    mkdir -p "$group_dir"

    # Process female bed files
    for female_bed in "${female_files[@]}"; do
        # Sort female bed file
        bedtools sort -i "$female_bed" > "${group_dir}/$(basename "$female_bed")_sortedf.bed"
    done

    # Process male bed files
    for male_bed in "${male_files[@]}"; do
        # Sort male bed file
        bedtools sort -i "$male_bed" > "${group_dir}/$(basename "$male_bed")_sortedm.bed"
    done

    # Count the number of female and male files
    num_female=$(ls -1 "${group_dir}"/*_sortedf.bed | wc -l)
    num_male=$(ls -1 "${group_dir}"/*_sortedm.bed | wc -l)

    # Find common peak regions in female bed files
    bedtools multiinter -i "${group_dir}"/*_sortedf.bed | awk -v num="$num_female" '$4 == num' | cut -f 1-4 > "${group_dir}/female_common_peaks.bed"

    # Find common peak regions in male bed files
    bedtools multiinter -i "${group_dir}"/*_sortedm.bed | awk -v num="$num_male" '$4 == num' | cut -f 1-4 > "${group_dir}/male_common_peaks.bed"

    # Intersect original bed files with common peaks for females
    for file in "${group_dir}"/*_sortedf.bed; do
        filename=$(basename "$file" _sortedf.bed)
        bedtools intersect -a "$file" -b "${group_dir}/female_common_peaks.bed" -wa -u -f 0.50 > "${group_dir}/${filename}_specific_females.bed"
    done

    # Intersect original bed files with common peaks for males
    for file in "${group_dir}"/*_sortedm.bed; do
        filename=$(basename "$file" _sortedm.bed)
        bedtools intersect -a "$file" -b "${group_dir}/male_common_peaks.bed" -wa -u -f 0.50 > "${group_dir}/${filename}_specific_males.bed"
    done

    # Combine all female peaks into one file
    cat "${group_dir}"/*_specific_females.bed > "${group_dir}/female_peaks.bed"

    # Combine all male peaks into one file
    cat "${group_dir}"/*_specific_males.bed > "${group_dir}/male_peaks.bed"

    # Find female specific peaks by subtracting male peaks from female peaks
    bedtools subtract -A -a "${group_dir}/female_peaks.bed" -b "${group_dir}/male_peaks.bed" -f 0.50 > "${group_dir}/female_specific_peaks.bed"

    # Find male specific peaks by subtracting female peaks from male peaks
    bedtools subtract -A -a "${group_dir}/male_peaks.bed" -b "${group_dir}/female_peaks.bed" -f 0.50 > "${group_dir}/male_specific_peaks.bed"

    # Clean up intermediate files
    rm -f "${group_dir}"/*_sorted*.bed "${group_dir}"/*_common_peaks*.bed "${group_dir}"/*.bed_specific_*.bed
done < "$input_file"

echo "Script execution complete."