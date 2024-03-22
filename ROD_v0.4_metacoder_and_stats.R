### TESTING METACODER FOR ROD
library(tidyverse)
library(wesanderson)
library(ggplot2)
library(metacoder)
library(ape)
library(ggridges)

# Read the current version: 
ROD_v0.4 <- read_rds("./../04_github/ROD/ROD_v0.4.rds") %>% as_tibble()
ROD_v0.4_genome_stats <- read_rds("./../04_github/ROD/ROD_v0.4_genome_stats.rds") %>% as_tibble()


#### General Stats ####
# Total number of rDNA copies
sum(ROD_v0.4_genome_stats$rDNA_copies)
sum(ROD_v0.4$size)

# Counting copies: 
# Pr.genome: 
# copy pr. genome, 1
ROD_v0.4_genome_stats %>% filter(rDNA_copies == 1) %>% 
  summarize(total_rDNA_copies = sum(rDNA_copies),average_rDNA_copies = mean(rDNA_copies))

# copy pr. genome, 2-10
ROD_v0.4_genome_stats %>% filter(rDNA_copies >= 2) %>% filter(rDNA_copies <= 10) %>% 
  summarize(total_rDNA_copies = sum(rDNA_copies),average_rDNA_copies = mean(rDNA_copies))

# copy pr. genome, 11-100
ROD_v0.4_genome_stats %>%  filter(rDNA_copies >= 11) %>% filter(rDNA_copies <= 100) %>%
  summarize(total_rDNA_copies = sum(rDNA_copies),average_rDNA_copies = mean(rDNA_copies))

# copy pr. genome, 101-1000
ROD_v0.4_genome_stats %>%  filter(rDNA_copies >= 101) %>% filter(rDNA_copies <= 1000) %>%
  summarize(total_rDNA_copies = sum(rDNA_copies),average_rDNA_copies = mean(rDNA_copies))

ROD_v0.4_genome_stats %>%  filter(rDNA_copies >= 1001) %>% #filter(rDNA_copies <= 1000) %>%
  summarize(total_rDNA_copies = sum(rDNA_copies), average_rDNA_copies = mean(rDNA_copies))

## Lenght pr. genome
## THIS IS NOW ADDED TO THE GENOME STATS DATA
# group_of_interest <- "assembly_id"
# df <- ROD_v0.4 %>%
#   #arrange(.data[[group_of_interest]], desc(size)) %>%
#   group_by(.data[[group_of_interest]]) %>%
#   reframe(rDNA_variants = n(),
#           length_min = min(length),
#           length_max= max(length),
#           total_size = sum(size),
#           size_max = max(size),
#           size_min = min(size),
#           size_second_max = ifelse(n_distinct(size) > 1, sort(size, decreasing = TRUE)[2], NA),
#           size_max_second_diff = size_max - ifelse(n_distinct(size) > 1, sort(size, decreasing = TRUE)[2], NA),
#           size_max_prop =  size_max / total_size,
#           size_second_prop =  size_second_max / total_size,
#           remaining_size = total_size - size_max - size_second_max,
#           all_sizes = paste(size, collapse = ";"),
#           all_lengths = paste(length, collapse = ";")) %>%
# ungroup()
# ROD_v0.4_genome_stats <- left_join(ROD_v0.4_genome_stats, df)
# saveRDS(ROD_v0.4_genome_stats, "./../04_github/ROD/ROD_v0.4_genome_stats.rds")

# Max copynumber in genome
(ROD_v0.4_genome_stats %>% filter(total_size == max(total_size)))$species
max(ROD_v0.4_genome_stats$length_max)
min(ROD_v0.4_genome_stats$length_min)
sum(ROD_v0.4_genome_stats$rDNA_copies>1)/length(ROD_v0.4_genome_stats$rDNA_copies)

# Lenght pr. other groups
group_of_interest <- "species"
df <- ROD_v0.4 %>% 
  #arrange(.data[[group_of_interest]], desc(size)) %>%
  group_by(.data[[group_of_interest]]) %>%
  reframe(rDNA_variants = n(),
          genomes = length(unique(assembly_id)),
            total_size = sum(size),
            size_max = max(size),
            size_min = min(size), 
            size_second_max = ifelse(n_distinct(size) > 1, sort(size, decreasing = TRUE)[2], NA), 
            size_max_second_diff = size_max - ifelse(n_distinct(size) > 1, sort(size, decreasing = TRUE)[2], NA),
            size_max_prop =  size_max / total_size,
            size_second_prop =  size_second_max / total_size, 
            remaining_size = total_size - size_max - size_second_max, 
            all_sizes = paste(size, collapse = ";"),
            all_lengths = paste(length, collapse = ";"),
            supergroup = paste(unique(supergroup),collapse = ";"),
            division = paste(unique(division),collapse = ";"),
            subdivision = paste(unique(subdivision),collapse = ";"),
            class = paste(unique(class),collapse = ";"),
            order = paste(unique(order),collapse = ";"),
            family = paste(unique(family),collapse = ";"),
            genus = paste(unique(genus),collapse = ";"),
            species = paste(unique(species),collapse = ";")
            ) %>% ungroup() 

View(df[842,])
head(df$species)
length(df$species)
sum(df$genomes)

ROD_v0.4_genome_stats %>% filter(species =="Brassica juncea")

# Order data according to some variable
ordered_data <- df %>% 
#  filter(rDNA_variants ==1 ) %>% 
  arrange(desc(total_size))

ordered_data
#df$group_of_interest <- factor(df$group_of_interest, levels = df$group_of_interest)
df[[group_of_interest]] <- factor(df[[group_of_interest]], levels = unique(df[[group_of_interest]]))


# Make data long, for plotting
long_data <- ordered_data %>%
  pivot_longer(cols = c(size_max, size_second_max, remaining_size), names_to = "category", values_to = "value")

# Plotting (only makes sense if group_of_interest = "assembly_id")
ggplot(long_data, aes(x = .data[[group_of_interest]], y = value, fill = category)) +
  geom_bar(stat = "identity") +  # Stacked bar plot
  theme_minimal() +
  labs(title = "Stacked Bar Plot of Sizes by Assembly ID",
       x = "Assembly ID",
       y = "Total Size",
       fill = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Adjusting x labels for readability

# long_data <- ordered_data %>%
#   pivot_longer(cols = c(size_max, size_min), names_to = "variable", values_to = "value")
# 
# # Plotting
# ggplot(long_data, aes(x = assembly_id, y = value, color = variable, group = variable)) +
#   geom_line() +  # Using geom_line() to connect points; you can change to geom_point() if you prefer dots
#   theme_minimal() +
#   labs(title = "Comparison of Maximum Size and Difference to Second Max",
#        x = "Assembly ID",
#        y = "Value",
#        color = "Variable") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Adjusting x labels for readability


# RIDGEPLOT
lengths_separated <- df %>% 
  separate_rows(all_lengths, sep = ";") %>%
  mutate(all_lengths = as.numeric(all_lengths))

lengths_separated

lengths_separated$lineage <- paste(lengths_separated$supergroup, lengths_separated$division,
                                   lengths_separated$subdivision, 
                                   lengths_separated$class,
                                   sep = " - ") %>% as.factor()


ggplot(lengths_separated, aes(x = all_lengths, y = lineage, fill = supergroup)) +
  geom_density_ridges2(rel_min_height = 0.005) +
  theme_ridges(line_size = 0.5, font_size = 6, grid = TRUE) +  # Use a suitable theme for ridge plots
  labs(title = "Ridge Plot of Lengths by class",
       x = "Length",
       y = group_of_interest) +
  scale_fill_manual(values = wes_palette("FantasticFox1",6,type = "continuous")) +
  scale_y_discrete(limits = rev(levels(lengths_separated$lineage)))
  # theme(legend.position = "bottom")
ggsave("ridge_class.pdf")



# Genomes and rDNA copy numbers pr. taxonomic divisions: 
# colnames(ROD_v0.4_genome_stats)
# 
# ROD_v0.4 %>%
#   group_by(division) %>%
#   summarise(
#     count = length(unique(assembly_id)), 
#     rDNA_copies_sum = sum(size), 
#     mean_rDNA_copies = mean(size)
#   ) %>%
#   ungroup()


#### Metacoder ####

# ROD_metacoder <- parse_tax_data(ROD_v0.4, class_cols = 7:13 , named_by_rank = TRUE)
ROD_metacoder <- parse_tax_data(ROD_v0.4_genome_stats, class_cols = 3:10 , named_by_rank = TRUE)
ROD_metacoder <- parse_tax_data(ROD_v0.4_genome_stats, class_cols = 3:8 , named_by_rank = TRUE)
ROD_metacoder <- parse_tax_data(ROD_v0.4_genome_stats, class_cols = 3:7 , named_by_rank = TRUE)

obj <- ROD_metacoder
obj$data$type_abund <- calc_taxon_abund(obj, "tax_data", cols = "rDNA_copies")

print(obj)
obj$data$tax_data


# Testing some palettes
wes_palette("Royal1")
wes_palette("Zissou1")
wes_palette("Rushmore1")

pal2 <- wes_palette(12, name = "Rushmore1", type = "continuous")
pal2

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



# BottleRocket1, BottleRocket2, Rushmore1, Royal1, Royal2, Zissou1, Darjeeling1, Darjeeling2, Chevalier1 , 
# FantasticFox1 , Moonrise1, Moonrise2, Moonrise3, Cavalcanti1, GrandBudapest1, GrandBudapest2, 
# IsleofDogs1, IsleofDogs2, FrenchDispatch, AsteroidCity2, AsteroidCity2, AsteroidCity3


