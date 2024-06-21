### TESTING METACODER FOR ROD
library(tidyverse)
library(wesanderson)
library(ggplot2)
library(metacoder)

# Read the current version: 
# ROD_v1.1 <- read_rds("./../04_github/ROD/ROD_v1.1.rds") %>% as_tibble()
ROD_v1.1_genome_stats <- read_rds("./../04_github/ROD/ROD_v1.1_genome_stats.rds") %>% as_tibble()
ROD_v1.1_genome_stats$mean_distance
ROD_v1.1_genome_stats$median_distance[is.na(ROD_v1.1_genome_stats$median_distance)] <- 0


#### Metacoder ####
# ROD_metacoder <- parse_tax_data(ROD_v1.1_genome_stats, class_cols = 3:10 , named_by_rank = TRUE)
# ROD_metacoder_family <- parse_tax_data(ROD_v1.1_genome_stats, class_cols = 3:9 , named_by_rank = TRUE)
# ROD_metacoder_order <- parse_tax_data(ROD_v1.1_genome_stats, class_cols = 3:8 , named_by_rank = TRUE)
# ROD_metacoder_class <- parse_tax_data(ROD_v1.1_genome_stats, class_cols = 3:7 , named_by_rank = TRUE)
ROD_metacoder_subdivision <- parse_tax_data(ROD_v1.1_genome_stats, class_cols = 3:6 , named_by_rank = TRUE)
# ROD_metacoder_division <- parse_tax_data(ROD_v1.1_genome_stats, class_cols = 3:5 , named_by_rank = TRUE)

obj <- ROD_metacoder_subdivision
obj$data$rDNA_copies <- calc_taxon_abund(obj, "tax_data", cols = "rDNA_copies")
obj$data$rDNA_variants <- calc_taxon_abund(obj, "tax_data", cols = "rDNA_variants")
obj$data$mean_distance <- calc_taxon_abund(obj, "tax_data", cols = "mean_distance")
obj$data$mean_distance$mean_distance/obj$n_obs()

print(obj)
obj$data$tax_data


# Testing some palettes
wes_palette("Royal1")
wes_palette("Zissou1")
wes_palette("Rushmore1")


pal <- wes_palette(15, name = "GrandBudapest1", type = "continuous")#[1:4]
pal[1:8]

pal <- wes_palette("GrandBudapest1")#[1:4]
# pal2 <- wes_palette(12, name = "Rushmore1", type = "continuous")
pal2 <- wes_palette("Rushmore1")
pal3 <- wes_palette("Zissou1")


set.seed(666) # This makes the plot appear the same each time it is run 
heat_tree(obj, 
          node_label = taxon_names,
          edge_label = obj$n_obs(),
          # node_size = obj$data$rDNA_variants$rDNA_variants/obj$n_obs(), 
          node_color = obj$data$rDNA_variants$rDNA_variants/obj$n_obs(), 
          node_size_axis_label = "rDNA variants (size)",
          # node_color = obj$data$mean_distance$mean_distance/obj$n_obs(),
          edge_color = obj$data$mean_distance$mean_distance/obj$n_obs(),
          node_color_axis_label = "mean genetic distance",
          edge_size = obj$n_obs(),
          edge_size_axis_label = "number of genomes",
          # edge_color = obj$data$rDNA_copies$rDNA_copies/obj$n_obs(),
          # node_color = obj$data$rDNA_copies$rDNA_copies/obj$n_obs(),
          node_size = obj$data$rDNA_copies$rDNA_copies/obj$n_obs(),
          edge_color_axis_label = "rDNA copies (col)",
          node_color_range = c(pal[1],pal[3]),
          edge_color_range = c(pal3[2],pal2[4]), 
          edge_label_max = 1500,
          node_label_max = 1500,
          edge_label_color = "white",
          node_size_trans = "area",
          # layout ="davidson-harel", # The primary layout algorithm. 
          # Setting layot to "reingold-tilford" and removing the initial layout results in a cricel
          layout ="reingold-tilford")
          #initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations
# ggsave("metacoder_2_subdivision_GenDist_v5.pdf") 

