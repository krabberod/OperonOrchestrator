
### RDIGEPLOT FOR ROD
library(tidyverse)
library(wesanderson)
library(ggplot2)
library(metacoder)
library(ggridges)

# Read the current version: 
# ROD_v1.1 <- read_rds("./../04_github/ROD/ROD_v1.1.rds") %>% as_tibble()
# ROD_v1.1_genome_stats <- read_rds("./../04_github/ROD/ROD_v1.1_genome_stats.rds") %>% as_tibble()



group_of_interest <- "subdivision"
df <- ROD_v1.1 %>% 
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

lengths_separated <- df %>% 
  separate_rows(all_lengths, sep = ";") %>%
  mutate(all_lengths = as.numeric(all_lengths))

lengths_separated

lengths_separated$lineage <- paste(lengths_separated$supergroup, lengths_separated$division,
                                   lengths_separated$subdivision, 
                                   sep = " - ") %>% as.factor()


#### RIDGEPLOT SUBDIVISION
# Assuming lengths_separated is your original data frame
# Calculate the number of all_lengths entries per lineage
counts <- lengths_separated %>%
  group_by(lineage) %>%
  summarize(count = n())

# Join the counts back to the original data frame
lengths_separated_with_counts <- lengths_separated %>%
  left_join(counts, by = "lineage") %>%
  mutate(lineage_with_count = paste(lineage, ", (n=", count, ")", sep = ""))

lengths_separated_with_counts$lineage_with_count <- factor(lengths_separated_with_counts$lineage_with_count, levels = rev(unique(lengths_separated_with_counts$lineage_with_count)))


combined_palette <- c(wes_palette("AsteroidCity1"), wes_palette("Cavalcanti1"))
# Create a data frame where each row represents a color
palette_df <- data.frame(color = combined_palette, index = 1:length(combined_palette))

# Plotting the combined palette
ggplot(palette_df, aes(x = factor(index), y = 1, fill = color)) + 
  geom_bar(stat = "identity") + 
  scale_fill_identity() + 
  theme_minimal() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(title = "Combined Color Palette")


# Update your ggplot code to use lineage_with_count for the y axis
p <- ggplot(lengths_separated_with_counts, aes(x = all_lengths, y = lineage_with_count, fill = supergroup)) +
  geom_density_ridges2(rel_min_height = 0.005) +
  theme_ridges(line_size = 0.5, font_size = 6, grid = TRUE) +
  labs(title = "Ridge Plot of Lengths by subdivision",
       x = "Length",
       y = group_of_interest) +
  scale_fill_manual(values = combined_palette) +
  scale_y_discrete(limits = rev(sort(levels(lengths_separated_with_counts$lineage_with_count))))

p <- p + geom_vline(xintercept = c(5779), linetype = "dashed", color = "red", linewidth = 0.5)

# Add points for single-entry lineages
single_entries_with_counts <- lengths_separated_with_counts %>%
  filter(count < 3)
p + geom_point(data = single_entries_with_counts, aes(x = all_lengths, y = lineage_with_count), color = "black", size = 2)


# ggsave("01_ROD_results/ridge_subdiv_v1.1.pdf")

