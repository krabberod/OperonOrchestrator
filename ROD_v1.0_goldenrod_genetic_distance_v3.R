library(tidyverse)
library(ggplot2)
library(ape)
library(msa)
source("02_scripts/calculate_weighted_distances.R")


ROD_v1.0 <- read_rds("./../04_github/ROD/ROD_v1.0_operon_variants.rds") %>% as_tibble()
ROD_v1.0_genome_stats <- read_rds("./../04_github/ROD/ROD_v1.0_genome_statistics.rds") %>% as_tibble()
# saveRDS(ROD_v1.0_genome_stats,"./../04_github/ROD/ROD_v1.0_genome_stats.rds")
head(ROD_v1.0$assembly_id)

# setwd("/Users/anderkkr/Dropbox/Projects/00_Master_projects/01_Active_Projects/02_OMG/17_OperonOrchestrator/05_goldenROD/")

# # Find the one where a msa is needed
# multicopy <- ROD_v1.0_genome_stats %>% filter(rDNA_variants > 1)
# assembly_id_list <- multicopy$assembly_id
# ROD_v1.0 %>% filter(assembly_id %in% assembly_id_list)
# genomes <- df$assembly_id %>% unique()
# Break down by genus? 
# df <- ROD_v1.0
# genomes <- unique(multicopy$assembly_id )[1:10]



### GENETIC DISTANCE PR. ACCESSION
# Read alignments: 
all_alignments <- sort(list.files("/Users/anderkkr/Dropbox/Projects/00_Master_projects/01_Active_Projects/02_OMG/17_OperonOrchestrator/07_goldenROD_alignment/02_pr_genome/01_mafft/", full.names = TRUE))
all_alignments <- sort(list.files("/Users/anderkkr/Dropbox/Projects/00_Master_projects/01_Active_Projects/02_OMG/17_OperonOrchestrator/test_aling/", full.names = TRUE))
all_alignments <- sort(list.files("/Users/anderkkr/Dropbox/Projects/00_Master_projects/01_Active_Projects/02_OMG/17_OperonOrchestrator/05_goldenROD/pr_subdivision/mafft_big/", full.names = TRUE))

# TESTING one by one: 
alignment <- readDNAMultipleAlignment(all_alignments[2])
basename(all_alignments[7]) %>% str_remove(".mafft.fasta") %>% str_remove("ROD_v1.0_")
alignment_ape <- as.DNAbin(alignment)
distances <- dist.dna(alignment_ape, model = "TN93", as.matrix = T)
distances
# If idetity is wanted:
1 - as.matrix(distances)
# distances <- 1 - as.matrix(distances)
# weighted_distance <- calculate_weighted_distances(distances)
# sd(distances)
# distances <- distances[lower.tri(distances, diag = FALSE)]
# sd(distances)

weighted_distance <- calculate_weighted_distances(distances)
1-weighted_distance
distances <-weighted_distance
# distances <- distances[lower.tri(distances, diag = FALSE)]
# distances <- distances[!is.nan(distances)]
which(is.nan(distances))
mean(1-distances)
median(distances)
max(distances)
min(distances)
max(distances) - min(distances)
var(distances)
max(distances)

# Loop through each accession, calculate summary stats, and store them in the list
# need 

distance_matrices <- list()
distance_matrices_weighted <- list()
summary_stats_list <- list()
source("02_scripts/calculate_weighted_distances.R")

# Set up progress bar
total <- length(all_alignments[1])
pb <- txtProgressBar(min = 0, max = total, style = 3)

# Initialize loop counter
counter <- 0

for(alignment in all_alignments[2]) {
  alignmentMsa <- readDNAMultipleAlignment(alignment)
  alignment_ape <- as.DNAbin(alignmentMsa)
  
  distances <- dist.dna(alignment_ape, model = "TN93", as.matrix = T)
  weighted_distance <- calculate_weighted_distances(distances)
  alignment_name <- basename(alignment) %>% str_remove(".mafft.fasta") %>% str_remove("ROD_v1.0_")
  distance_matrices_weighted[[alignment_name]] <- weighted_distance
  
  # distances <- distances[lower.tri(distances, diag = FALSE)]
  distances <- distances[lower.tri(weighted_distance, diag = FALSE)]
  
  # Calculate summary statistics
  mean_distance <- mean(distances)
  median_distance <- median(distances)
  distance_range <- max(distances) - min(distances)
  std_dev <- sd(distances)
  min_distance <- min(distances)
  max_distance <- max(distances)

  # Create a tibble for the current alignment with the summary stats
  alignment_summary <- tibble(
    assembly_id = alignment_name,
    mean_distance = mean_distance,
    min_distance = min_distance, 
    max_distance = max_distance,
    median_distance = median_distance,
    distance_range = distance_range,
    std_dev = std_dev
  
  )
  
  # Append the species summary to the list
  summary_stats_list[[basename(alignment)]] <- alignment_summary
  # Update counter and progress bar
  counter <- counter + 1
  setTxtProgressBar(pb, counter)
}

distance_matrices_weighted[2]

# Combine all species summaries into a single tibble
summary_stats <- bind_rows(summary_stats_list)
summary_stats
colnames(summary_stats)[1] <- "subdivision"

write.table(summary_stats, "./01_ROD_results/class_distance_summary_stats_2.tab", row.names = F, quote = F, sep = "\t")
# saveRDS(distance_matrices_weighted, "./../01_results/class_distance_matrices_weighted_1.rds")
saveRDS(distance_matrices_weighted, "./../01_results/class_distance_matrices_weighted_2.rds")

### FINDING THE REPRESENTATIV ROD


# Assuming 'summary_stats' is a data frame with species names and their corresponding summary statistics
# summary_stats <- data.frame(species = c(...), mean_distance = c(...), median_distance = c(...), range = c(...), variance = c(...))

# Boxplot for visualizing the distribution of mean distances across species
ggplot(summary_stats, aes(x = subdivision, y = mean_distance)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Mean Genetic Distance Distribution by genome", x = "genome", y = "Mean Genetic Distance")

# Heatmap for one species (example)
# Assuming 'distance_matrix_species1' is your distance matrix for a specific species
library(pheatmap)
pheatmap(distance_matrices_weighted[[3]], color = colorRampPalette(c("blue", "white", "red"))(50), 
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
         main = "Genetic Distance Heatmap for Species 1")

# Bar chart for comparing mean distances across species
ggplot(summary_stats, aes(x = subdivision, y = mean_distance, fill = alignment)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, ),legend.position = "none") +
  labs(title = "Mean Genetic Distance Comparison Across genome", x = "genome", y = "Mean Genetic Distance")


# Scatter plot (if you have an additional variable, e.g., habitat diversity)
# Assuming 'habitat_diversity' is part of your summary_stats data frame
ggplot(summary_stats, aes(x = distance_range, y = mean_distance)) +
  geom_point(aes(color = subdivision)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, ),legend.position = "none") +
  labs(title = "Genetic Distance vs.", x = "distance_range", y = "Mean Genetic Distance")


#### CALCULATING THE GENTIC DISTANCE PR.  

