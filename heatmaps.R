options(repos = c(CRAN = "https://cran.r-project.org"))

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) < 4) {
  stop("Usage: Rscript script_name.R <input_file> <output_directory> <genome_file> <chromosome_size>")
}

# Input file path
input_file <- args[1]
cat(input_file)

# Output directory
output_dir <- args[2]

# Genome file
genome_file <- args[3]

# Chromosome size file
chromosome_size <- args[4]

# Read input file
input_data <- read.table(input_file, header = FALSE, stringsAsFactors = TRUE, sep = "\t")
print(input_data)

# Create output directory if it doesn't exist
if (!file.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to install R packages if not already installed
install_r_packages <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg)
    }
  }
}

# Function to install external tools if not already installed
install_external_tools <- function(tools) {
  for (tool in tools) {
    if (system(paste("which", tool), intern = TRUE) == "") {
      if (system("conda --version", intern = TRUE) != "") {
        system(paste("conda install -c bioconda", tool))
      } else {
        print(paste("ERROR:", tool, "is not installed and conda is not available. Please install it manually."))
        quit(save = "no", status = 1)
      }
    }
  }
}

# Load necessary R libraries
print("Loading libraries...")
required_packages <- c("circlize", "rtracklayer", "GenomicRanges", "EnrichedHeatmap", "data.table", "readr", "ComplexHeatmap", "ggplot2")
install_r_packages(required_packages)

# Install external tools
tools <- c("bedtools", "samtools", "ucsc-bedgraphtobigwig")
install_external_tools(tools)


# Function to create heatmaps for a group
create_heatmaps <- function(female_peaks_bed, male_peaks_bed, group, output_dir, genome_file, chromosome_size) {
  cat("Creating heatmaps for group:", group, "\n")
  
  # Create directory for the group
  group_dir <- file.path(output_dir, group)
  if (!file.exists(group_dir)) {
    dir.create(group_dir, showWarnings = FALSE)
  }
  
 
  # Run the data preprocessing commands for female peaks
  female_cmds <- c(
    paste0("bedtools bedtobam -i ", female_peaks_bed, " -g ", genome_file, " > ", file.path(group_dir, "female_peaks.bam")),
    paste0("samtools sort ", file.path(group_dir, "female_peaks.bam"), " -o ", file.path(group_dir, "female_peaks_sorted.bam")),
    paste0("samtools index ", file.path(group_dir, "female_peaks_sorted.bam")),
    paste0("bedtools genomecov -ibam ", file.path(group_dir, "female_peaks_sorted.bam"), " -bg > ", file.path(group_dir, "female_peaks.bedgraph")),
    paste0("bedGraphToBigWig ", file.path(group_dir, "female_peaks.bedgraph"), " ", chromosome_size, " ", file.path(group_dir, "female_peaks.bw"))
  )
  
  # Run the data preprocessing commands for male peaks
  male_cmds <- c(
    paste0("bedtools bedtobam -i ", male_peaks_bed, " -g ", genome_file, " > ", file.path(group_dir, "male_peaks.bam")),
    paste0("samtools sort ", file.path(group_dir, "male_peaks.bam"), " -o ", file.path(group_dir, "male_peaks_sorted.bam")),
    paste0("samtools index ", file.path(group_dir, "male_peaks_sorted.bam")),
    paste0("bedtools genomecov -ibam ", file.path(group_dir, "male_peaks_sorted.bam"), " -bg > ", file.path(group_dir, "male_peaks.bedgraph")),
    paste0("bedGraphToBigWig ", file.path(group_dir, "male_peaks.bedgraph"), " ", chromosome_size, " ", file.path(group_dir, "male_peaks.bw"))
  )
  
  # Execute the commands
  for (cmd in c(female_cmds, male_cmds)) {
    system(cmd)
  }
    
  print(female_peaks_bed)
   female_peaks <- paste0(female_peaks_bed)
   male_peaks <- paste0(male_peaks_bed)
    print(female_peaks)
  
  # Load data and create heatmaps
    bed_data_female <- read.table(female_peaks, header = FALSE, stringsAsFactors = TRUE, sep = "\t")
    bed_data_male <- read.table(male_peaks, header = FALSE, stringsAsFactors = TRUE, sep = "\t")
  female_bw <- import(file.path(group_dir, "female_peaks.bw"), format = "BigWig")
  male_bw <- import(file.path(group_dir, "male_peaks.bw"), format = "BigWig")
  
  # code from diffBind package
  bed_data_female <- bed_data_female[c(1:6, 10)]
  bed_data_male <- bed_data_male[c(1:6, 10)]
  # Write modified BED files
  write_delim(bed_data_female, file.path(group_dir, "female_specific_peaks_modified.bed"), delim = "\t", col_names = FALSE)
  write_delim(bed_data_male, file.path(group_dir, "male_specific_peaks_modified.bed"), delim = "\t", col_names = FALSE)

    
female_bed <- import(file.path(group_dir, "female_specific_peaks_modified.bed"), format = "BED")
male_bed <- import(file.path(group_dir, "male_specific_peaks_modified.bed"), format = "BED")


#code for female specific peaks heatmap
  female.10kb <- resize(female_bed, width = 10000, fix = "center")
  female.10kb.center<- resize(female.10kb, width =1, fix = "center")
  female.mat<- normalizeToMatrix(female_bw, female.10kb.center, value_column = "score",mean_mode="w0", w=100, extend = 5000)
  male.mat<- normalizeToMatrix(male_bw, female.10kb.center, value_column = "score",mean_mode="w0", w=100, extend = 5000)
  quantile(female.mat, probs = c(0.005, 0.5,0.90))
  quantile(male.mat, probs = c(0.005, 0.5,0.90))
  library(circlize)
  col_fun_female<- circlize::colorRamp2(c(0, 1), c("white", "red"))
  col_fun_male<- circlize::colorRamp2(c(0, 1), c("white", "red"))
  png(file=file.path(group_dir, "female_heatmap.png"))
  female_heatmap <- EnrichedHeatmap(female.mat, axis_name_rot = 0, name = "Female", column_title = "Female", use_raster = TRUE, col = col_fun_female, top_annotation = HeatmapAnnotation(lines = anno_enriched())) + EnrichedHeatmap(male.mat, axis_name_rot = 0, name = "Male", column_title = "Male", use_raster = TRUE, col = col_fun_male, top_annotation = HeatmapAnnotation(lines = anno_enriched()))
    draw(female_heatmap)
    dev.off()
  
#code for male specific peaks heatmap
  male.10kb<- resize(male_bed, width = 10000, fix = "center")
  male.10kb.center<- resize(male.10kb, width =1, fix = "center")
  male.mat<- normalizeToMatrix(male_bw, male.10kb.center, value_column = "score",mean_mode="w0", w=100, extend = 5000)
  female.mat<- normalizeToMatrix(female_bw, male.10kb.center, value_column = "score",mean_mode="w0", w=100, extend = 5000)
  quantile(male.mat, probs = c(0.005, 0.5,0.90))
  quantile(female.mat, probs = c(0.005, 0.5,0.90))
  col_fun_female<- circlize::colorRamp2(c(0, 1), c("white", "red"))
  col_fun_male<- circlize::colorRamp2(c(0, 1), c("white", "red"))
   png(file=file.path(group_dir, "male_heatmap.png"))
  male_heatmap <- EnrichedHeatmap(male.mat, axis_name_rot = 0, name = "Male",column_title = "Male", use_raster = TRUE, col = col_fun_male, top_annotation = HeatmapAnnotation(lines = anno_enriched())) + EnrichedHeatmap(female.mat, axis_name_rot = 0, name = "Female", column_title = "Female", use_raster = TRUE, col = col_fun_female,top_annotation = HeatmapAnnotation(lines = anno_enriched()))
    draw(male_heatmap)
    dev.off()

}

# Loop over each line in input data
for (i in 1:nrow(input_data)) {
  # Get input parameters
  female_peaks_bed <- input_data[i, "V1"]
  male_peaks_bed <- input_data[i, "V2"]
  group <- input_data[i, "V3"]
  
  # Create heatmaps for the group
  create_heatmaps(female_peaks_bed, male_peaks_bed, group, output_dir, genome_file, chromosome_size)
}

print("Script execution complete.")
