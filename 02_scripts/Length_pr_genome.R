library(tidyverse)


# Lenght pr. genome
# THIS IS NOW ADDED TO THE GENOME STATS DATA
group_of_interest <- "assembly_id"
df <- ROD_v0.4 %>%
  #arrange(.data[[group_of_interest]], desc(size)) %>%
  group_by(.data[[group_of_interest]]) %>%
  reframe(rDNA_variants = n(),
          length_min = min(length),
          length_max= max(length),
          total_size = sum(size),
          size_max = max(size),
          size_min = min(size),
          size_second_max = ifelse(n_distinct(size) > 1, sort(size, decreasing = TRUE)[2], NA),
          size_max_second_diff = size_max - ifelse(n_distinct(size) > 1, sort(size, decreasing = TRUE)[2], NA),
          size_max_prop =  size_max / total_size,
          size_second_prop =  size_second_max / total_size,
          remaining_size = total_size - size_max - size_second_max,
          all_sizes = paste(size, collapse = ";"),
          all_lengths = paste(length, collapse = ";")) %>%
ungroup()
ROD_v0.4_genome_stats <- left_join(ROD_v0.4_genome_stats, df)
saveRDS(ROD_v0.4_genome_stats, "./../04_github/ROD/ROD_v0.4_genome_stats.rds")