# Load necessary libraries
library(tidyverse)
library(factoextra)

# Define the directory containing the folders
base_dir <- "Desktop/ShinyApp/data/PRJNA894167_4_processed/fivepseq_counts"

# List all sample folders 
sample_folders <- list.dirs(base_dir, recursive = FALSE)

# Initialize empty list to store the samples
data_list <- list()

# Loop through each sample folder
for (folder in sample_folders) {
  # Construct the file path to the amino_acid_pauses.txt file
  file_path <- file.path(folder, "protein_coding", "amino_acid_pauses.txt")
  
  # Check if the file exists
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }
  
  # Read the file
  data <- read.table(file_path, header = TRUE, sep = "\t", check.names = FALSE)
  
  # Extract all amino acid vectors (excluding the first column, which contains amino acid names)
  amino_acid_data <- data[, -1]  # Remove the first column (amino acid names)
  rownames(amino_acid_data) <- data[, 1]  # Set amino acid names as row names
  
  # Store the data in the list with the sample name as the key
  sample_name <- basename(folder)
  data_list[[sample_name]] <- amino_acid_data
}

# Combine all amino acid vectors into a single data frame
# Each row will represent a sample, and each column will represent an amino acid at a specific position
combined_data <- do.call(rbind, lapply(names(data_list), function(sample) {
  # Flatten the amino acid data for each sample into a single row
  flattened_data <- as.vector(t(data_list[[sample]]))
  names(flattened_data) <- paste(rownames(data_list[[sample]]), colnames(data_list[[sample]]), sep = "_")
  flattened_data
}))

# Add sample names as row names
rownames(combined_data) <- names(data_list)

# Perform PCA on the combined data
pca_result <- prcomp(combined_data, scale. = TRUE)

# Visualize the PCA results
fviz_pca_ind(pca_result, 
             geom.ind = "point", # Show points only
             col.ind = rownames(combined_data), # Color by sample
             palette = "jco", # Use a nice color palette
             legend.title = "Samples"
) + 
  ggtitle("PCA of Amino Acid Pauses (All Amino Acids)") +
  theme_minimal()
