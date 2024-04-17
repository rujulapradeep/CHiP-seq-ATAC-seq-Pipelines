options(repos = c(CRAN = "https://cran.r-project.org"))

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 2) {
  stop("Usage: Rscript script_name.R <input_file> <output_directory>")
}

# Input file path
input_file <- args[1]

# Output directory
output_dir <- args[2]

# Read input file
input_data <- read.table(input_file, header = FALSE, stringsAsFactors = TRUE, sep = "\t")

# Create output directory if it doesn't exist
if (!file.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


#function to install package
install_required_packages <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      BiocManager::install(package)
    }
  }
}

# Define required packages
required_packages <- c(
  "BiocParallel", 
  "ChIPseeker", 
  "TxDb.Hsapiens.UCSC.hg38.knownGene", 
  "rtracklayer", 
  "org.Hs.eg.db", 
  "readr"
)

# Call the function to install required packages
install_required_packages(required_packages)

# Load necessary libraries
print("Loading libraries...")
library(BiocParallel)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(rtracklayer)
library(org.Hs.eg.db)
library(parallel)
library(readr)

# Function to annotate peaks
annotate_peaks <- function(peaks_file, group, gender, output_folder) {
  cat("Annotating peaks for group:", group, "\n")
  
  # Load TxDb
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
 # Determine the number of available cores
num_cores <- parallel::detectCores()


# Use all available cores for parallel processing
param <- MulticoreParam(num_cores)
print(paste("Using", num_cores, "cores..."))

register(param)

   #make peaks file into correct format
    peaks_file <- toString(peaks_file)
   peaks_table <- read.table(peaks_file, header = FALSE, stringsAsFactors = TRUE, sep = "\t")
   peaks_table_modified <- peaks_table[c(1:6, 10)]

  if (gender == "female") {
  write_delim(peaks_table_modified, file.path(group_folder, "female_specific_peaks_modified.bed"), delim = "\t", col_names = FALSE)
  regions_to_annotate <- import(file.path(group_folder, "female_specific_peaks_modified.bed"), format = "BED")
} else {
  write_delim(peaks_table_modified, file.path(group_folder, "male_specific_peaks_modified.bed"), delim = "\t", col_names = FALSE)
  regions_to_annotate <- import(file.path(group_folder, "male_specific_peaks_modified.bed"), format = "BED")
}


  # Using CHiP-seeker's annotatePeak function
  annotate_chunk <- function(chunk) {
    annotatePeak(chunk, tssRegion = c(-3000, 3000),
                 TxDb = txdb, level = "gene", annoDb = "org.Hs.eg.db",
                 sameStrand = FALSE, ignoreOverlap = FALSE, overlap = "TSS")
  }
 
  # Adjust the size of each chunk as needed
  chunk_size <- 100
  chunks <- split(regions_to_annotate, ceiling(seq_along(regions_to_annotate) / chunk_size))
  
 
  # Annotate peaks in parallel
  annotations <- bplapply(chunks, FUN = annotate_chunk)
  
  
  # Create a connection to the output file
  if (gender == "female") {
    output_file <- file.path(output_folder, paste0(group, "_female_annotations.txt"))
  } else {
    output_file <- file.path(output_folder, paste0(group, "_male_annotations.txt"))
  }
  conn <- file(output_file, open = "w")
  
  # Combine annotations from each chunk to one output file
  print("Writing annotations to output file...")
  for (i in seq_along(annotations)) {
    write.table(annotations[[i]], conn, append = TRUE, sep = "\t", quote = FALSE)
  }
  
  print("Finished...")
  
  # Close the connection
  close(conn)
}

# Loop over each line in input data
for (i in 1:nrow(input_data)) {
  # Get input parameters
  female_specific_peaks <- input_data[i, "V1"]
  print(female_specific_peaks)
  male_specific_peaks <- input_data[i, "V2"]
  group <- input_data[i, "V3"]
  
  # Create group-specific folder in output directory
  group_folder <- file.path(output_dir, paste0(group, "_annotations"))
  dir.create(group_folder, showWarnings = FALSE)
  
  # Run annotation for female peaks
  annotate_peaks(female_specific_peaks, group, "female", group_folder)
  
  # Run annotation for male peaks
  annotate_peaks(male_specific_peaks, group, "male", group_folder)
}

print("Script execution complete.")
