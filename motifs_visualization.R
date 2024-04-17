#!/usr/bin/env Rscript
options(repos = c(CRAN = "https://cran.r-project.org"))



# Define required packages
required_packages <- c(
  "htmlTable", "htmlwidgets"
)

# Function to install R packages if not already installed
install_r_packages <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, dependencies = TRUE)
    }
  }
}

# Load necessary R libraries
print("Loading libraries...")
install_r_packages(required_packages)

if (!requireNamespace("xfun", quietly = TRUE)) {
  install.packages("xfun", dependencies = TRUE)
}

if (!requireNamespace("htmlTable", quietly = TRUE)) {
  install.packages("htmlTable", dependencies = TRUE, force = TRUE)
}
library(htmlTable)


if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
  install.packages("htmlwidgets", dependencies = TRUE, force = TRUE)
}
library(htmlwidgets)


# Check if correct number of command-line arguments are provided
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript motifs_visualization.R <input_file> <output_directory>")
}

# Assign command-line arguments to variables
input_file <- args[1]
output_dir <- args[2]

# Check if output directory exists, create if not
if (!file.exists(output_dir)) {
  dir.create(output_dir)
  cat("Output directory created:", output_dir, "\n")
} 

# Read input file
input_data <- read.table(input_file, header = FALSE, stringsAsFactors = FALSE)

# Process data for each group
for (i in 1:nrow(input_data)) {
  
  group <- input_data[i, 1]
  female_file <- input_data[i, 2]
  male_file <- input_data[i, 3]
  
  group_dir <- file.path(output_dir, paste0(group, "_motifs"))
  
  
  # Create group directory if it doesn't exist
  if (!file.exists(group_dir)) {
    dir.create(group_dir)
    cat("Group directory created:", group_dir, "\n")
  }
  
  cat("Visualizing motifs for:", group, "\n")
  
 #function to create circle for female motifs
  create_circle <- function(log_p_value) {
    if (is.finite(log_p_value)) {
      scaled_value <- 1 - pmax(0, pmin(1, (log_p_value - min(female_motif_data$V4)) / (max(female_motif_data$V4) - min(female_motif_data$V4))))
      min_size <- 5
      size <- min_size + scaled_value * (30 - min_size)
      circle_html <- sprintf('<div style="width: %fpx; height: %fpx; border-radius: 50%%; background-color: #3498db;"></div>', size, size)
      return(circle_html)
    } 
  }
  #function to create circle for male motifs
  create_circle_male <- function(log_p_value) {
    if (is.finite(log_p_value)) {
      scaled_value <- 1 - pmax(0, pmin(1, (log_p_value - min(male_motif_data$V4)) / (max(male_motif_data$V4) - min(male_motif_data$V4))))
      min_size <- 5
      size <- min_size + scaled_value * (30 - min_size)
      circle_html <- sprintf('<div style="width: %fpx; height: %fpx; border-radius: 50%%; background-color: #e74c3c;"></div>', size, size)
      return(circle_html)
    } 
  }
  
  # Read data for female motifs and create circles
  female_motif_data <- read.table(female_file, header = FALSE, sep = "\t")
  female_motif_data$Motif <- sub("/.*", "", female_motif_data$V1)
  female_motif_data$Circle <- sapply(female_motif_data$V4, create_circle)
  
  # Read data for male motifs and create circles
  male_motif_data <- read.table(male_file, header = FALSE, sep = "\t")
  male_motif_data$Motif <- sub("/.*", "", male_motif_data$V1)
  male_motif_data$Circle <- sapply(male_motif_data$V4, create_circle_male)
  
  # Combine unique motifs
  unique_motifs <- unique(c(female_motif_data$Motif, male_motif_data$Motif))
  
  # Create combined data frame
  combined_data <- data.frame(Motif = unique_motifs)
  combined_data$Female_Circles <- sapply(unique_motifs, function(motif) paste(female_motif_data$Circle[female_motif_data$Motif == motif], collapse = ""))
  combined_data$Male_Circles <- sapply(unique_motifs, function(motif) paste(male_motif_data$Circle[male_motif_data$Motif == motif], collapse = ""))
  combined_data <- combined_data[order(nchar(combined_data$Female_Circles) > 0 & nchar(combined_data$Male_Circles) > 0, nchar(combined_data$Male_Circles) > 0, decreasing = TRUE), ]
  combined_data$Motif_Num <- seq_len(nrow(combined_data))
  combined_data <- combined_data[, c("Motif_Num", "Motif", "Female_Circles", "Male_Circles")]
  
  # Create HTML table
  combined_table_html <- htmlTable(combined_data, header = c("#", "Motif", "Female Circles", "Male Circles"), align = c("c", "l", "c", "c"))
  
  # Print or display the combined HTML table
  print(combined_table_html)
  
  # Save HTML table to group directory
  output_file <- file.path(group_dir, "motifs_table.html")
  cat(as.character(combined_table_html), file = output_file)
  cat("Table saved to:", output_file, "\n")
}

print("Script execution complete.")
