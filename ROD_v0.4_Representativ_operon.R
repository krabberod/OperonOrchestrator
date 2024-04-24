### FINDING THE REPRESENTATIV OPERON IN  ROD
library(tidyverse)
library(ggplot2)

ROD_v0.4 <- read_rds("./../04_github/ROD/ROD_v0.4.rds") %>% as_tibble()
ROD_v0.4_genome_stats <- read_rds("./../04_github/ROD/ROD_v0.4_genome_stats.rds") %>% as_tibble()


# write.table(ROD_v0.4_genome_stats, "ROD_v0.4_genome_stats.tab", quote = F, sep = "\t", row.names = F)
# write.table(ROD_v0.4, "ROD_v0.4.tab", quote = F, sep = "\t", row.names = F)


# Criteria: the representative need to be "more common" than the other operons
# One thing to check easily: If the size_max > 75%
# ROD_v0.4_genome_stats %>% select(size_max,size_min,size_second_max,size_max_prop)
# ROD_v0.4_genome_stats %>% filter(size_max_prop < 0.5)  %>% select(size_max,size_min,size_second_max,size_max_prop,total_size)
# ROD_v0.4_genome_stats %>% filter(size_max_prop > 0.6) %>% select(taxid) %>% unique()
# ROD_v0.4_genome_stats %>% filter(size_max_prop < 0.6) %>% filter(size_max < size_second_max) %>% select(size_max,size_min,size_second_max,size_max_prop)
ROD_v0.4_genome_stats %>% filter(rDNA_copies != size_max)
similar_size <- ROD_v0.4_genome_stats %>% filter(size_max > 1) %>% filter(size_max == size_second_max) # %>% .$rDNA_variants %>% sum()

similar_size_ass <- similar_size$assembly_id
# Pick the representative sequence = max_size
# Print the rest


max_size_tibble_unique <- ROD_v0.4 %>% filter(!assembly_id %in% similar_size_ass) %>%
  group_by(assembly_id) %>%
  filter(size == max(size)) %>%
  dplyr::slice(1) %>%
  ungroup()

non_max_size_tibble <- anti_join(ROD_v0.4, max_size_tibble_unique, by = c("assembly_id"))
write.table(similar_size, "./../08_rexRODney_manually_curated/for_manuall_curation.tab", quote = F, sep = "\t", row.names = F)
ROD_v0.4_genome_stats %>% filter(assembly_id == "GCA_001541825") %>% .$max_distance

colnames(ROD_v0.4_genome_stats)
ROD_v0.4_genome_stats <- ROD_v0.4_genome_stats %>% select(-c(size_second_max,size_max_second_diff,size_second_prop,remaining_size),)
write.table(max_size_tibble_unique, "./../04_github/ROD/ROD_v1.0_reference_sequences.tab", quote = F, sep = "\t", row.names = F)
write.table(ROD_v0.4_genome_stats, "./../04_github/ROD/ROD_v1.0_genome_statistics.tab", quote = F, sep = "\t", row.names = F)
write.table(ROD_v0.4, "./../04_github/ROD/ROD_v1.0_operon_variants.tab", quote = F, sep = "\t", row.names = F)
saveRDS(ROD_v0.4_genome_stats, "./../04_github/ROD/ROD_v1.0_genome_statistics.rds")
saveRDS(max_size_tibble_unique, "./../04_github/ROD/ROD_v1.0_reference_sequences.rds")
