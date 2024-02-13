### TESTING METACODER FOR ROD
library(tidyverse)
library(wesanderson)
library(ggplot2)
library(metacoder)
library(ape)


ROD_v0.4 <- read_rds("./../04_github/ROD/ROD_v0.4.rds") %>% as_tibble()
ROD_v0.4_genome_stats <- read_rds("./../04_github/ROD/ROD_v0.4_genome_stats.rds") %>% as_tibble()
# ROD_v0.4_genome_stats <- 

#colnames(ROD_v0.4)
# ROD_v0.4_genome_stats <- ROD_v0.4 %>% select(-c(seqid, length, size, sequence)) %>% unique() %>% as_tibble()
#sum(ROD_v0.4$size)


#total_rdna_size <- ROD_v0.4 %>% 
#group_by(assembly_id) %>%
# summarise(rDNA_copies = sum(size), .groups = 'drop')
#ROD_v0.4_genome_stats <- left_join(ROD_v0.4_genome_stats,total_rdna_size) 
(dim(ROD_v0.4_genome_stats %>% filter(rDNA_copies == 2)))


ROD_metacoder <- parse_tax_data(ROD_v0.4, class_cols = 7:13 , named_by_rank = TRUE)
ROD_metacoder <- parse_tax_data(ROD_v0.4_genome_stats, class_cols = 3:10 , named_by_rank = TRUE)
ROD_metacoder <- parse_tax_data(ROD_v0.4_genome_stats, class_cols = 3:8 , named_by_rank = TRUE)
ROD_metacoder <- parse_tax_data(ROD_v0.4_genome_stats, class_cols = 3:7 , named_by_rank = TRUE)

obj <- ROD_metacoder
obj$data$type_abund <- calc_taxon_abund(obj, "tax_data", cols = "rDNA_copies")

print(obj)
obj$data$tax_data


pal <- wes_palette(15, name = "GrandBudapest1", type = "continuous")#[1:4]
pal[1:8]
# pal2 <- wes_palette(15, name = "Zissou1", type = "continuous")

pal2 <- wes_palette(12, name = "Rushmore1", type = "continuous")
pal2



set.seed(666) # This makes the plot appear the same each time it is run 
heat_tree(obj, 
          node_label = taxon_names,
          node_size = obj$data$type_abund$rDNA_copies/obj$n_obs(),
          node_color = obj$data$type_abund$rDNA_copies/obj$n_obs(), 
          node_color_range = pal[1:8], 
          edge_color_range = pal2[5:8], 
          # edge_label_color_range = quantative_palette(),
          edge_size = obj$n_obs(),
          edge_color = obj$n_obs(),
          edge_label = obj$n_obs(),
          node_size_axis_label = "rDNA copies",
          node_color_axis_label = "rDNA copies pr. genome",
          edge_size_axis_label = "number of genomes",
          edge_color_axis_label = "number of genomes",
          edge_label_max = 1500,
          node_label_max = 1500,
          edge_label_color = "white",
          node_size_trans = "area",
          # layout ="davidson-harel", # The primary layout algorithm. 
          # Setting layot to "reingold-tilford" and removing the initial layout results in a cricel
          layout ="reingold-tilford")
          # initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations
# ggsave("metacoder_class_genome_pr.pdf") 
ggsave("metacoder_class_genome_pr_wheel.pdf") 

wes_palette("Royal1")
wes_palette("Zissou1",2)
wes_palette("Rushmore1")

pal2 <- wes_palette(12, name = "Rushmore1", type = "continuous")
pal2


# BottleRocket1, BottleRocket2, Rushmore1, Royal1, Royal2, Zissou1, Darjeeling1, Darjeeling2, Chevalier1 , 
# FantasticFox1 , Moonrise1, Moonrise2, Moonrise3, Cavalcanti1, GrandBudapest1, GrandBudapest2, 
# IsleofDogs1, IsleofDogs2, FrenchDispatch, AsteroidCity2, AsteroidCity2, AsteroidCity3


