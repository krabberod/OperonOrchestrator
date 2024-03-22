library(tidyverse)
library(ggplot2)
library(ape)
library(msa)

ROD_v0.4 <- read_rds("./../04_github/ROD/ROD_v0.4.rds") %>% as_tibble()
ROD_v0.4_genome_stats <- read_rds("./../04_github/ROD/ROD_v0.4_genome_stats.rds") %>% as_tibble()
head(ROD_v0.4$assembly_id)

# # Find the one where a msa is needed
# multicopy <- ROD_v0.4_genome_stats %>% filter(rDNA_variants > 1)
# assembly_id_list <- multicopy$assembly_id
# ROD_v0.4 %>% filter(assembly_id %in% assembly_id_list)
# genomes <- df$assembly_id %>% unique()
# Break down by genus? 
# df <- ROD_v0.4
# 
# # set the output directory 
# for (genome in genomes){
#   df_sub <- df %>% dplyr::filter(assembly_id==genome) 
#   # header <- paste0(df_subf$seqid,";","size=",df_sub$size,"|",df_sub$lineage , " taxid=",df_sub$taxid, ";")
#   header <- paste0(df_sub$seqid,";","size=",df_sub$size)
#   Xfasta <- character(nrow(df_sub) * 2)
#   Xfasta[c(TRUE, FALSE)] <- paste0(">", header)
#   Xfasta[c(FALSE, TRUE)] <- df_sub$sequence
#   writeLines(Xfasta, paste0("./../05_goldenROD/pr_genome_multicopy/ROD_v0.4_", genome ,".fasta"))
#   getwd()
#   #print(header)
#   print(paste0("./../05_goldenROD/pr_genome_multicopy/ROD_v0.4_", genome ,".fasta"))
# }

genomes <- unique(multicopy$assembly_id )[1:10]
distance_matrices <- list()

for(genome in genomes){
  temp <- ROD_v0.4 %>% filter(assembly_id == genome) 
  dna_sequences <- DNAStringSet(temp$sequence)
  alignment <- msa(dna_sequences)
  alignment_ape <- as.DNAbin(alignment)
  distances <- dist.dna(alignment_ape, model = "JC") 
  distance_matrices[[genome]] <- distances
}


summary_stats_list <- list()

# Loop through each species, calculate summary stats, and store them in the list
for(genome in names(distance_matrices)) {
  # Extract the distance matrix for the current species
  distances <- as.vector(as.matrix(distance_matrices[[genome]]))
  
  # Removing diagonal and upper triangle as they are redundant or zero
  distances <- distances[lower.tri(distance_matrices[[genome]])]
  
  # Calculate summary statistics
  mean_distance <- mean(distances)
  median_distance <- median(distances)
  distance_range <- max(distances) - min(distances)
  variance_distance <- var(distances)
  
  # Create a tibble for the current species with the summary stats
  genome_summary <- tibble(
    genome = genome,
    mean_distance = mean_distance,
    median_distance = median_distance,
    distance_range = distance_range,
    variance_distance = variance_distance
  )
  
  # Append the species summary to the list
  summary_stats_list[[genome]] <- genome_summary
}

# Combine all species summaries into a single tibble
summary_stats <- bind_rows(summary_stats_list)

# Print the summary statistics tibble
print(summary_stats)




# Assuming 'summary_stats' is a data frame with species names and their corresponding summary statistics
# summary_stats <- data.frame(species = c(...), mean_distance = c(...), median_distance = c(...), range = c(...), variance = c(...))

# Boxplot for visualizing the distribution of mean distances across species
ggplot(summary_stats, aes(x = genome, y = mean_distance)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Mean Genetic Distance Distribution by genome", x = "genome", y = "Mean Genetic Distance")

# Heatmap for one species (example)
# Assuming 'distance_matrix_species1' is your distance matrix for a specific species
library(pheatmap)
pheatmap(distance_matrices[[10]], color = colorRampPalette(c("blue", "white", "red"))(50), 
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
         main = "Genetic Distance Heatmap for Species 1")

# Bar chart for comparing mean distances across species
ggplot(summary_stats, aes(x = genome, y = mean_distance, fill = genome)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Mean Genetic Distance Comparison Across genome", x = "genome", y = "Mean Genetic Distance")

# Scatter plot (if you have an additional variable, e.g., habitat diversity)
# Assuming 'habitat_diversity' is part of your summary_stats data frame
# ggplot(summary_stats, aes(x = habitat_diversity, y = mean_distance)) +
#   geom_point(aes(color = genome)) +
#   labs(title = "Genetic Distance vs. Habitat Diversity", x = "Habitat Diversity", y = "Mean Genetic Distance")


