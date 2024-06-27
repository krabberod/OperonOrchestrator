### FINDING THE REPRESENTATIV OPERON IN  ROD
library(tidyverse)
library(ggplot2)

ROD_v1.2 <- read_rds("./../04_github/ROD/ROD_v1.2.rds") %>% as_tibble()
ROD_v1.2_genome_stats <- read_rds("./../04_github/ROD/ROD_v1.2_genome_stats.rds") %>% as_tibble()
# write.table(ROD_v1.2_genome_stats, "ROD_v1.2_genome_stats.tab", quote = F, sep = "\t", row.names = F)
# write.table(ROD_v1.2, "ROD_v1.2.tab", quote = F, sep = "\t", row.names = F)


# Criteria: the representative need to be "more common" than the other operons
# One thing to check easily: If the size_max > 75%
second_max) %>% select(size_max,size_min,size_second_max,size_max_prop)
ROD_v1.2_genome_stats %>% filter(rDNA_copies != size_max)
similar_size <- ROD_v1.2_genome_stats %>% filter(size_max > 1) %>% filter(size_max == size_second_max) # %>% .$rDNA_variants %>% sum()

similar_size_ass <- similar_size$assembly_id
# Pick the representative sequence = max_size
# Print the rest

ROD_v1.2_reference_sequences <- ROD_v1.2 %>% #filter(!assembly_id %in% similar_size_ass) %>%
  group_by(assembly_id) %>%
  filter(size == max(size)) %>%
  dplyr::slice(1) %>%
  ungroup()

# non_max_size_tibble <- anti_join(ROD_v1.2, max_size_tibble_unique, by = c("assembly_id"))
# write.table(similar_size, "./../08_rexRODney_manually_curated/for_manuall_curation.tab", quote = F, sep = "\t", row.names = F)
# ROD_v1.2_genome_stats %>% filter(assembly_id == "GCA_001541825") %>% .$max_distance

colnames(ROD_v1.2_genome_stats)
ROD_v1.2_genome_stats <- ROD_v1.2_genome_stats %>% select(-c(size_second_max,size_max_second_diff,size_second_prop,remaining_size),)
write.table(ROD_v1.2_reference_sequences, gzfile("./../04_github/ROD/ROD_v1.2_reference_sequences.tab.gz"), quote = FALSE, sep = "\t", row.names = FALSE)


# Ensure the seqid columns are character type for the join operation
ROD_v1.2_reference_sequences <- ROD_v1.2_reference_sequences %>%
  mutate(seqid = as.character(seqid))

ROD_v1.2 <- ROD_v1.2 %>%
  mutate(seqid = as.character(seqid))

# Perform the left join and create the reference column
ROD_v1.2 <- ROD_v1.2 %>%
  left_join(ROD_v1.2_reference_sequences %>%
              select(seqid) %>%
              mutate(reference = "Yes"), by = "seqid") %>%
  mutate(reference = ifelse(is.na(reference), "No", reference)) %>%
  select(seqid, reference, everything())

